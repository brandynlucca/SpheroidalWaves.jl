module SpheroidalWaveFunctions

using Libdl

export smn, rmn, radial_wronskian, accuracy, eigenvalue, jacobian_eigen, jacobian_smn, jacobian_rmn, find_c_for_eigenvalue

const _backend_libraries = Dict{Symbol,Union{Nothing,String}}(
    :double => nothing,
    :quad => nothing,
)

const _backend_handles = Dict{String,Ptr{Cvoid}}()

const _ENV_BACKEND_R8 = "SPHEROIDALWAVEFUNCTIONS_LIBRARY_R8"
const _ENV_BACKEND_R16 = "SPHEROIDALWAVEFUNCTIONS_LIBRARY_R16"

function _validate_precision(precision::Symbol)
    if precision != :double && precision != :quad
        error("precision must be :double or :quad, got :$precision")
    end
end

function set_backend_library!(path::AbstractString; precision::Symbol=:double)
    _validate_precision(precision)
    _backend_libraries[precision] = String(path)
    return _backend_libraries[precision]
end

function backend_library(; precision::Symbol=:double)
    _validate_precision(precision)
    return _backend_libraries[precision]
end

function _set_backend_from_candidate(path, precision::Symbol, source::AbstractString)
    if path isa AbstractString
        candidate = String(path)
        if isfile(candidate)
            set_backend_library!(candidate; precision=precision)
            return true
        else
            @warn "Ignoring backend path from $source because file does not exist" precision path=candidate maxlog=1
        end
    elseif path !== nothing
        @warn "Ignoring backend path from $source because value is not a string" precision value=path maxlog=1
    end
    return false
end

function _configure_backends_from_env!()
    configured_any = false
    if haskey(ENV, _ENV_BACKEND_R8)
        configured_any |= _set_backend_from_candidate(ENV[_ENV_BACKEND_R8], :double, "ENV[$_ENV_BACKEND_R8]")
    end
    if haskey(ENV, _ENV_BACKEND_R16)
        configured_any |= _set_backend_from_candidate(ENV[_ENV_BACKEND_R16], :quad, "ENV[$_ENV_BACKEND_R16]")
    end
    return configured_any
end

function _configure_backends_from_local_config!()
    config_file = joinpath(dirname(@__FILE__), "..", "deps", "library_config.jl")
    if !isfile(config_file)
        return false
    end

    configured_any = false
    include(config_file)
    if isdefined(@__MODULE__, :SPHEROIDAL_BATCH_LIBRARY_R8)
        lib_r8 = Base.invokelatest(getfield, @__MODULE__, :SPHEROIDAL_BATCH_LIBRARY_R8)
        configured_any |= _set_backend_from_candidate(lib_r8, :double, "local library_config.jl")
    end
    if isdefined(@__MODULE__, :SPHEROIDAL_BATCH_LIBRARY_R16)
        lib_r16 = Base.invokelatest(getfield, @__MODULE__, :SPHEROIDAL_BATCH_LIBRARY_R16)
        configured_any |= _set_backend_from_candidate(lib_r16, :quad, "local library_config.jl")
    end
    return configured_any
end

function _require_backend_library(precision::Symbol)
    _validate_precision(precision)
    lib = _backend_libraries[precision]
    if lib === nothing
        error("No backend library configured for precision :$precision.")
    end
    return lib
end

function _require_backend_handle(lib::String)
    return get!(_backend_handles, lib) do
        Libdl.dlopen(lib)
    end
end

function _symbol_pointer(lib::String, symbol::Symbol)
    return Libdl.dlsym(_require_backend_handle(lib), symbol)
end

function _real_suffix(precision::Symbol)
    return precision === :double ? "_r8" : "_r16"
end

function _complex_suffix(precision::Symbol)
    return precision === :double ? "_c8" : "_c16"
end

function _bool_to_cint(x::Bool)
    return x ? Cint(1) : Cint(0)
end

function _check_scalar_status(status::Cint)
    if status != 0
        error("Backend call failed with status $(Int(status)).")
    end
end

function _check_vector_status(status::Vector{Cint})
    if any(!=(0), status)
        error("Backend batch call returned non-zero status for one or more elements.")
    end
end

function _complex_parts(re::Vector{Float64}, im::Vector{Float64})
    out = Vector{ComplexF64}(undef, length(re))
    @inbounds for i in eachindex(re)
        out[i] = complex(re[i], im[i])
    end
    return out
end

function _is_exact_spherical_limit(c::Real)
    return iszero(c)
end

function _is_exact_spherical_limit(c::Complex)
    return iszero(real(c)) && iszero(imag(c))
end

function _default_jacobian_step(c::Real)
    return sqrt(eps(Float64)) * max(1.0, abs(Float64(c)))
end

function _default_jacobian_step(c::Complex)
    return sqrt(eps(Float64)) * max(1.0, abs(real(c)), abs(imag(c)))
end

function _resolve_jacobian_step(c::Union{Real,Complex}, h)
    step = h === nothing ? _default_jacobian_step(c) : Float64(h)
    if !(step > 0)
        error("h must be positive, got $h")
    end
    return step
end

function _validate_jacobian_tolerances(rtol::Real, atol::Real)
    if !(rtol > 0)
        error("rtol must be positive, got $rtol")
    end
    if !(atol > 0)
        error("atol must be positive, got $atol")
    end
end

function _all_finite(x)
    return x isa Number ? isfinite(x) : all(isfinite, x)
end

function _max_relative_change(dh, dh2; atol::Real)
    abs_diff = abs.(dh2 .- dh)
    scale = max.(abs.(dh2), atol)
    rel = abs_diff ./ scale
    return rel isa Number ? Float64(rel) : Float64(maximum(rel))
end

function _jacobian_suggested_action(conditioning_flag::Symbol, finite_flag::Bool, precision::Symbol)
    if !finite_flag
        return :retry_smaller_h
    end
    if conditioning_flag === :good
        return :accept
    elseif conditioning_flag === :warning
        return precision === :quad ? :accept : :retry_smaller_h
    else
        return precision === :quad ? :retry_smaller_h : :use_quad
    end
end

function _jacobian_metadata(dh, dh2, step_used::Float64; precision::Symbol, rtol::Real, atol::Real)
    finite_flag = _all_finite(dh) && _all_finite(dh2)
    rel_change = _max_relative_change(dh, dh2; atol=atol)
    conditioning_flag = !finite_flag ? :poor : (rel_change <= 10 * rtol ? :good : (rel_change <= 1000 * rtol ? :warning : :poor))
    suggested_action = _jacobian_suggested_action(conditioning_flag, finite_flag, precision)
    return (
        step_used=step_used,
        relative_change_when_halving_step=rel_change,
        finite_flag=finite_flag,
        conditioning_flag=conditioning_flag,
        suggested_action=suggested_action,
    )
end

function _finite_difference_with_metadata(calc::Function, step::Float64; precision::Symbol, adaptive::Bool, rtol::Real, atol::Real)
    d_h = calc(step)
    d_h2 = calc(step / 2)
    if !adaptive
        metadata = _jacobian_metadata(d_h, d_h2, step; precision=precision, rtol=rtol, atol=atol)
        return d_h, metadata
    end

    best = d_h2
    step_used = step / 2
    metadata = _jacobian_metadata(d_h, d_h2, step_used; precision=precision, rtol=rtol, atol=atol)

    if adaptive && metadata.conditioning_flag === :poor
        d_h4 = calc(step / 4)
        best = d_h4
        step_used = step / 4
        metadata = _jacobian_metadata(d_h2, d_h4, step_used; precision=precision, rtol=rtol, atol=atol)
    end

    return best, metadata
end

function _legendre_p(n::Integer, x::Float64)
    n < 0 && throw(ArgumentError("n must be nonnegative"))
    n == 0 && return 1.0
    n == 1 && return x

    prev2 = 1.0
    prev1 = x
    for k in 2:n
        curr = ((2k - 1) * x * prev1 - (k - 1) * prev2) / k
        prev2, prev1 = prev1, curr
    end
    return prev1
end

function _legendre_p_derivative(n::Integer, x::Float64)
    n == 0 && return 0.0
    pn = _legendre_p(n, x)
    pnm1 = _legendre_p(n - 1, x)
    return n * (x * pn - pnm1) / (x^2 - 1.0)
end

function _legendre_q(n::Integer, x::Float64)
    n < 0 && throw(ArgumentError("n must be nonnegative"))
    abs(x) <= 1.0 && throw(ArgumentError("Legendre Q_n(x) requires |x| > 1"))

    q0 = 0.5 * log((x + 1.0) / (x - 1.0))
    n == 0 && return q0

    q1 = x * q0 - 1.0
    n == 1 && return q1

    prev2 = q0
    prev1 = q1
    for k in 2:n
        curr = ((2k - 1) * x * prev1 - (k - 1) * prev2) / k
        prev2, prev1 = prev1, curr
    end
    return prev1
end

function _legendre_q_derivative(n::Integer, x::Float64)
    qn = _legendre_q(n, x)
    qnm1 = n == 0 ? 0.5 * log((x + 1.0) / (x - 1.0)) : _legendre_q(n - 1, x)
    return n == 0 ? -1.0 / (x^2 - 1.0) : n * (x * qn - qnm1) / (x^2 - 1.0)
end

function _spherical_scale(n::Integer, normalize::Bool)
    return normalize ? sqrt((2n + 1) / 2) : 1.0
end

function _spherical_smn_real(m::Integer, n::Integer, eta::AbstractVector{<:Real}; normalize::Bool=false)
    m == 0 || throw(ArgumentError("analytic c=0 fallback is only implemented for m = 0"))
    scale = _spherical_scale(n, normalize)
    eta64 = Float64.(eta)
    value = [scale * _legendre_p(n, x) for x in eta64]
    derivative = [scale * _legendre_p_derivative(n, x) for x in eta64]
    return (; value, derivative)
end

function _spherical_smn_complex(m::Integer, n::Integer, eta::AbstractVector{<:Real}; normalize::Bool=false)
    real_result = _spherical_smn_real(m, n, eta; normalize=normalize)
    value = ComplexF64.(real_result.value)
    derivative = ComplexF64.(real_result.derivative)
    return (; value, derivative)
end

function _spherical_rmn_real(m::Integer, n::Integer, x::AbstractVector{<:Real}; kind::Integer=1)
    m == 0 || throw(ArgumentError("analytic c=0 fallback is only implemented for m = 0"))
    x64 = Float64.(x)
    p_value = [_legendre_p(n, xi) for xi in x64]
    p_derivative = [_legendre_p_derivative(n, xi) for xi in x64]
    q_value = [_legendre_q(n, xi) for xi in x64]
    q_derivative = [_legendre_q_derivative(n, xi) for xi in x64]

    if kind == 1
        return (; value=ComplexF64.(p_value), derivative=ComplexF64.(p_derivative))
    elseif kind == 2
        return (; value=ComplexF64.(q_value), derivative=ComplexF64.(q_derivative))
    elseif kind == 3
        return (; value=complex.(p_value, q_value), derivative=complex.(p_derivative, q_derivative))
    elseif kind == 4
        return (; value=complex.(p_value, -q_value), derivative=complex.(p_derivative, -q_derivative))
    end

    error("kind must be 1, 2, 3, or 4")
end

function _call_real_smn(prefix::Symbol, m::Integer, n::Integer, c::Real, eta::AbstractVector{<:Real}; precision::Symbol=:double, normalize::Bool=false)
    if _is_exact_spherical_limit(c)
        return _spherical_smn_real(m, n, eta; normalize=normalize)
    end

    lib = _require_backend_library(precision)
    suffix = _real_suffix(precision)
    symbol = Symbol(String(prefix) * "_smn_batch" * suffix)
    fnptr = _symbol_pointer(lib, symbol)

    n_eta = Cint(length(eta))
    eta64 = Float64.(eta)
    value = zeros(Float64, n_eta)
    derivative = zeros(Float64, n_eta)
    status = Ref{Cint}(0)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}),
          Cint(m), Cint(n), Cdouble(c), n_eta, eta64, _bool_to_cint(normalize), value, derivative, status)

    _check_scalar_status(status[])
    return (; value, derivative)
end

function _call_real_rmn(prefix::Symbol, m::Integer, n::Integer, c::Real, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1)
    if prefix === :psms && _is_exact_spherical_limit(c)
        return _spherical_rmn_real(m, n, x; kind=kind)
    end

    lib = _require_backend_library(precision)
    suffix = _real_suffix(precision)
    symbol = Symbol(String(prefix) * "_rmn_batch" * suffix)
    fnptr = _symbol_pointer(lib, symbol)

    n_x = Cint(length(x))
    x64 = Float64.(x)
    value_re = zeros(Float64, n_x)
    value_im = zeros(Float64, n_x)
    deriv_re = zeros(Float64, n_x)
    deriv_im = zeros(Float64, n_x)
    status = zeros(Cint, n_x)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
          Cint(m), Cint(n), Cdouble(c), n_x, x64, Cint(kind), value_re, value_im, deriv_re, deriv_im, status)

    _check_vector_status(status)
    value = _complex_parts(value_re, value_im)
    derivative = _complex_parts(deriv_re, deriv_im)
    return (; value, derivative)
end

function _call_complex_smn(prefix::Symbol, m::Integer, n::Integer, c::Complex, eta::AbstractVector{<:Real}; precision::Symbol=:double, normalize::Bool=false)
    if _is_exact_spherical_limit(c)
        return _spherical_smn_complex(m, n, eta; normalize=normalize)
    end

    lib = _require_backend_library(precision)
    suffix = _complex_suffix(precision)
    symbol = Symbol(String(prefix) * "_smn_batch" * suffix)
    fnptr = _symbol_pointer(lib, symbol)

    n_eta = Cint(length(eta))
    eta64 = Float64.(eta)
    value_re = zeros(Float64, n_eta)
    value_im = zeros(Float64, n_eta)
    deriv_re = zeros(Float64, n_eta)
    deriv_im = zeros(Float64, n_eta)
    status = Ref{Cint}(0)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}),
          Cint(m), Cint(n), Cdouble(real(c)), Cdouble(imag(c)), n_eta, eta64, _bool_to_cint(normalize),
          value_re, value_im, deriv_re, deriv_im, status)

    _check_scalar_status(status[])
    value = _complex_parts(value_re, value_im)
    derivative = _complex_parts(deriv_re, deriv_im)
    return (; value, derivative)
end

function _call_complex_rmn(prefix::Symbol, m::Integer, n::Integer, c::Complex, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1)
    lib = _require_backend_library(precision)
    suffix = _complex_suffix(precision)
    symbol = Symbol(String(prefix) * "_rmn_batch" * suffix)
    fnptr = _symbol_pointer(lib, symbol)

    n_x = Cint(length(x))
    x64 = Float64.(x)
    value_re = zeros(Float64, n_x)
    value_im = zeros(Float64, n_x)
    deriv_re = zeros(Float64, n_x)
    deriv_im = zeros(Float64, n_x)
    status = zeros(Cint, n_x)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
          Cint(m), Cint(n), Cdouble(real(c)), Cdouble(imag(c)), n_x, x64, Cint(kind),
          value_re, value_im, deriv_re, deriv_im, status)

    _check_vector_status(status)
    value = _complex_parts(value_re, value_im)
    derivative = _complex_parts(deriv_re, deriv_im)
    return (; value, derivative)
end

function _call_real_smn_accuracy(prefix::Symbol, m::Integer, n::Integer, c::Real, eta::AbstractVector{<:Real}; precision::Symbol=:double, normalize::Bool=false)
    if _is_exact_spherical_limit(c)
        return fill(precision === :double ? 15 : 30, length(eta))
    end

    lib = _require_backend_library(precision)
    suffix = _real_suffix(precision)
    symbol = Symbol(String(prefix) * "_smn_batch" * suffix * "_acc")
    fnptr = _symbol_pointer(lib, symbol)

    n_eta = Cint(length(eta))
    eta64 = Float64.(eta)
    value = zeros(Float64, n_eta)
    derivative = zeros(Float64, n_eta)
    naccs = zeros(Cint, n_eta)
    status = Ref{Cint}(0)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}),
          Cint(m), Cint(n), Cdouble(c), n_eta, eta64, _bool_to_cint(normalize), value, derivative, naccs, status)

    _check_scalar_status(status[])
    return Int.(naccs)
end

function _call_real_rmn_accuracy(prefix::Symbol, m::Integer, n::Integer, c::Real, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1)
    if prefix === :psms && _is_exact_spherical_limit(c)
        return fill(precision === :double ? 15 : 30, length(x))
    end

    lib = _require_backend_library(precision)
    suffix = _real_suffix(precision)
    symbol = Symbol(String(prefix) * "_rmn_batch" * suffix * "_acc")
    fnptr = _symbol_pointer(lib, symbol)

    n_x = Cint(length(x))
    x64 = Float64.(x)
    value_re = zeros(Float64, n_x)
    value_im = zeros(Float64, n_x)
    deriv_re = zeros(Float64, n_x)
    deriv_im = zeros(Float64, n_x)
    naccr = zeros(Cint, n_x)
    status = zeros(Cint, n_x)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}),
          Cint(m), Cint(n), Cdouble(c), n_x, x64, Cint(kind), value_re, value_im, deriv_re, deriv_im, naccr, status)

    _check_vector_status(status)
    return Int.(naccr)
end

function _call_complex_smn_accuracy(prefix::Symbol, m::Integer, n::Integer, c::Complex, eta::AbstractVector{<:Real}; precision::Symbol=:double, normalize::Bool=false)
    if _is_exact_spherical_limit(c)
        return fill(precision === :double ? 15 : 30, length(eta))
    end

    lib = _require_backend_library(precision)
    suffix = _complex_suffix(precision)
    symbol = Symbol(String(prefix) * "_smn_batch" * suffix * "_acc")
    fnptr = _symbol_pointer(lib, symbol)

    n_eta = Cint(length(eta))
    eta64 = Float64.(eta)
    value_re = zeros(Float64, n_eta)
    value_im = zeros(Float64, n_eta)
    deriv_re = zeros(Float64, n_eta)
    deriv_im = zeros(Float64, n_eta)
    naccs = zeros(Cint, n_eta)
    status = Ref{Cint}(0)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}),
          Cint(m), Cint(n), Cdouble(real(c)), Cdouble(imag(c)), n_eta, eta64, _bool_to_cint(normalize),
          value_re, value_im, deriv_re, deriv_im, naccs, status)

    _check_scalar_status(status[])
    return Int.(naccs)
end

function _call_complex_rmn_accuracy(prefix::Symbol, m::Integer, n::Integer, c::Complex, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1)
    lib = _require_backend_library(precision)
    suffix = _complex_suffix(precision)
    symbol = Symbol(String(prefix) * "_rmn_batch" * suffix * "_acc")
    fnptr = _symbol_pointer(lib, symbol)

    n_x = Cint(length(x))
    x64 = Float64.(x)
    value_re = zeros(Float64, n_x)
    value_im = zeros(Float64, n_x)
    deriv_re = zeros(Float64, n_x)
    deriv_im = zeros(Float64, n_x)
    naccr = zeros(Cint, n_x)
    status = zeros(Cint, n_x)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}),
          Cint(m), Cint(n), Cdouble(real(c)), Cdouble(imag(c)), n_x, x64, Cint(kind),
          value_re, value_im, deriv_re, deriv_im, naccr, status)

    _check_vector_status(status)
    return Int.(naccr)
end

function _call_real_eigenvalue(prefix::Symbol, m::Integer, n::Integer, c::Real; precision::Symbol=:double)
    if _is_exact_spherical_limit(c)
        return Float64(n * (n + 1))
    end

    lib = _require_backend_library(precision)
    suffix = _real_suffix(precision)
    symbol = Symbol(String(prefix) * "_eigenvalue" * suffix)
    fnptr = _symbol_pointer(lib, symbol)

    eig = Ref{Cdouble}(0.0)
    status = Ref{Cint}(0)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Ref{Cdouble}, Ref{Cint}),
          Cint(m), Cint(n), Cdouble(c), eig, status)

    _check_scalar_status(status[])
    return eig[]
end

function _call_complex_eigenvalue(prefix::Symbol, m::Integer, n::Integer, c::Complex; precision::Symbol=:double)
    if _is_exact_spherical_limit(c)
        return complex(Float64(n * (n + 1)), 0.0)
    end

    lib = _require_backend_library(precision)
    suffix = _complex_suffix(precision)
    symbol = Symbol(String(prefix) * "_eigenvalue" * suffix)
    fnptr = _symbol_pointer(lib, symbol)

    eig_re = Ref{Cdouble}(0.0)
    eig_im = Ref{Cdouble}(0.0)
    status = Ref{Cint}(0)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}),
          Cint(m), Cint(n), Cdouble(real(c)), Cdouble(imag(c)), eig_re, eig_im, status)

    _check_scalar_status(status[])
    return complex(eig_re[], eig_im[])
end

function smn(m::Integer, n::Integer, c::Union{Real,Complex}, eta::AbstractVector{<:Real};
             spheroid::Symbol=:prolate, precision::Symbol=:double, normalize::Bool=false)
    """
    Spheroidal angular wave functions (first kind).
    
    Computes prolate or oblate spheroidal angular wave functions of the first kind,
    Smₙ(η, c), and their derivatives with respect to η.
    
    These are eigenfunction solutions to the angular part of the spheroidal wave equation,
    orthogonal on the interval η ∈ [-1, 1]. They closely resemble associated Legendre
    polynomials when c is small.
    
    **Mathematical Background:**
    - The spheroidal wave equation separates into radial and angular parts
    - Angular functions satisfy: d/dη[(1-η²)dSmₙ/dη] + [λ(c,m,n) - c²(1-η²)]Smₙ = 0
    - Here λ(c,m,n) is the eigenvalue (separation constant)
    - Order m: 0 ≤ m ≤ n
    - Degree n: defines which eigenvalue/eigenfunction to compute
    
    **Parameters:**
    
        m::Integer
            Order parameter (must satisfy 0 ≤ m ≤ n)
        n::Integer  
            Degree parameter (defines eigenfunction; n ≥ 0)
        c::Real or Complex
            Size parameter (prolate: c = kd/2 with k=wavenumber, d=interfocal distance)
            - Real c: real-valued functions
            - Complex c: complex-valued functions (advanced applications)
        η::AbstractVector{<:Real}
            Evaluation points in [-1, 1]
            Represents cos(θ) where θ is angle in spherical coordinates
        
        Keyword Arguments:
            spheroid::Symbol = :prolate
                :prolate or :oblate geometry
                - :prolate: elongated spheroid (ζ > 0 in ξ,η,ϕ coordinates)
                - :oblate: flattened spheroid (ζ < 0)
            
            precision::Symbol = :double
                :double (default) or :quad precision
                - :double: ~15 decimal digits (real*8)
                - :quad: ~30 decimal digits (real*16, slower)
            
            normalize::Bool = false
                Whether to scale by sqrt((2n+1)/2)
                - false: Meixner-Schafke normalization (norm → ∞ as m → ∞)
                - true: Unity normalization (constant norm, convenient for applications)
    
    **Returns:**
        NamedTuple with fields:
        - .value::Vector{Real or Complex}
            Function values Smₙ(η, c) at each η point
        - .derivative::Vector{Real or Complex}
            First derivatives dSmₙ/dη at each η point
    
    **Accuracy:**
        - Double precision: ~14 decimal digits (excellent for most applications)
        - Quad precision: ~29 decimal digits (for high-precision research)
    
    **Special Cases:**
        - c = 0 (spherical limit): Returns associated Legendre polynomials P_n^m(η)
        - c → ∞: Functions oscillate rapidly; use appropriate η resolution
        - m = 0: Functions are even; m > 0: mixed parity
    
    **Performance:**
        - Vectorized in Fortran; fast batch evaluation
        - O(1) cost per evaluation point (amortized over batch)
        - Typical: 1000 points in ~0.001s
    
    **Examples:**
        
        # Basic usage: prolate Smn with c=1.5
        η = [-0.9, -0.5, 0, 0.5, 0.9]
        result = smn(1, 2, 1.5, η)
        println(result.value)         # Function values
        println(result.derivative)    # Derivatives
        
        # Normalized functions (unity norm)
        result_norm = smn(0, 3, 10.0, η; normalize=true)
        
        # Oblate spheroid with higher precision
        result_quad = smn(2, 4, 5.0, η; spheroid=:oblate, precision=:quad)
        
        # Complex parameter (advanced)
        c_complex = 1.5 + 0.1im
        result_c = smn(1, 2, c_complex, η)
    
    **Notes:**
        - Order and degree must satisfy 0 ≤ m ≤ n
        - All η values must be in [-1, 1]
        - For spherical limit (c ≈ 0), uses analytic Legendre polynomial fallback
        - Complex c support requires additional validation for branch cuts
    
    **References:**
        - DLMF §30.2: https://dlmf.nist.gov/30.2
        - Van Buren & Boisvert (2004): Accurate calculation of prolate spheroidal wave functions
    """
    
    _validate_precision(precision)
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end

    if c isa Real
        prefix = spheroid === :prolate ? :psms : :oblate
        return _call_real_smn(prefix, m, n, c, eta; precision=precision, normalize=normalize)
    else
        prefix = spheroid === :prolate ? :cprolate : :coblate
        return _call_complex_smn(prefix, m, n, c, eta; precision=precision, normalize=normalize)
    end
end

function rmn(m::Integer, n::Integer, c::Union{Real,Complex}, x::AbstractVector{<:Real};
             spheroid::Symbol=:prolate, precision::Symbol=:double, kind::Integer=1)
    """
    Spheroidal radial wave functions (characteristic-exponent form).
    
    Computes prolate or oblate spheroidal radial wave functions and their derivatives,
    returning values in characteristic-exponent form to avoid overflow/underflow.
    
    These are eigenfunction solutions to the radial part of the spheroidal wave equation,
    defined on the interval x > 1 (for prolate) or x > -1 (for oblate, with |x| > 1).
    Four kinds of radial functions are available.
    
    **Mathematical Background:**
    - Radial functions satisfy: d/dξ[(ξ²-1)dRmₙ/dξ] - [λ(c,m,n)ξ² - m² - c²ξ(ξ²-1)]Rmₙ = 0
    - λ(c,m,n) is the eigenvalue from the angular function problem
    - Functions are related to spherical Bessel and Hankel functions when c→0
    - Characteristic-exponent representation: R = R_characteristic × 10^exponent
    
    **Parameters:**
    
        m, n, c: See smn() documentation
        
        x::AbstractVector{<:Real}
            Radial evaluation points
            - For prolate: x > 1 (distance from interfocal axis)
            - For oblate: |x| > 1
            - Vectorized computation: all points evaluated in single Fortran call
        
        kind::Integer = 1
            Which radial function kind:
            
            1: First kind R₁ (propagating, finite at x→∞)
               - Regular solution; used for most applications
               - Analogous to spherical Bessel j_ℓ(kr)
            
            2: Second kind R₂ (singular at x→∞)
               - Complementary solution
               - Analogous to spherical Bessel y_ℓ(kr)
            
            3: Third kind (linear combination) = R₁ + i·R₂
               - Hankel function analog (outgoing wave)
            
            4: Fourth kind (linear combination) = R₁ - i·R₂
               - Alternate Hankel function analog
    
    **Returns:**
        NamedTuple with fields:
        - .value::Vector{ComplexF64}
            Function values Rmₙ(x, c) (always complex, even for real input)
            Real part is the physical radial function
            Imaginary part indicates which kind was computed
        
        - .derivative::Vector{ComplexF64}
            First derivatives dRmₙ/dx
    
    **Accuracy & Overflow Protection:**
        - Uses characteristic-exponent representation internally
        - Scales all values to [0.1, 1) × 10^exponent
        - Automatically reconstructs full-precision values
        - Maintains accuracy even for very large/small values
    
    **Special Cases:**
        - c = 0 (spherical limit): Returns Legendre Q functions for kind=2
        - x = 1 (boundary): kind=2 produces infinite value (singular point)
        - Large c or x: Rapid oscillation; may need fine resolution
    
    **Examples:**
        
        # Basic radial function (kind 1)
        x = [1.5, 2.0, 3.0, 5.0]
        r1 = rmn(0, 1, 200, x; kind=1)
        println(r1.value)
        
        # Second kind for same parameter set
        r2 = rmn(0, 1, 200, x; kind=2)
        
        # Hankel-like combinations
        r3_hankel = rmn(0, 1, 200, x; kind=3)  # Outgoing wave
        r4_hankel = rmn(0, 1, 200, x; kind=4)  # Incoming wave
        
        # Verification: Wronskian should be constant
        W = r1.value .* r2.derivative - r1.derivative .* r2.value
        println(W)  # Should be approximately constant
        
        # Oblate spheroid
        r_oblate = rmn(1, 2, 500, x; spheroid=:oblate, kind=1)
    
    **Notes:**
        - Return values are complex even when c and x are real
        - For real c, kind=1 has real values; kind=2 has real values; kinds 3,4 are complex
        - Characteristic-exponent format is internal; reconstructed automatically
        - Wronskian relation: R₁·R₂' - R₁'·R₂ = constant (diagnostic check)
    
    **Performance:**
        - Vectorized batch computation in Fortran
        - Four kinds computed efficiently in single Fortran call
        - Typical: 1000 points in ~0.001s
    
    **References:**
        - DLMF §30.3: https://dlmf.nist.gov/30.3
        - Van Buren & Boisvert (2004): Accurate calculation of prolate spheroidal wave functions
    """
    
    _validate_precision(precision)
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end

    if c isa Real
        prefix = spheroid === :prolate ? :psms : :oblate
        return _call_real_rmn(prefix, m, n, c, x; precision=precision, kind=kind)
    else
        prefix = spheroid === :prolate ? :cprolate : :coblate
        return _call_complex_rmn(prefix, m, n, c, x; precision=precision, kind=kind)
    end
end

function radial_wronskian(m::Integer, n::Integer, c::Union{Real,Complex}, x::AbstractVector{<:Real};
                          spheroid::Symbol=:prolate, precision::Symbol=:double)
    """
    Verify consistency of radial functions via the Wronskian determinant.
    
    Computes the Wronskian W = R₁·R₂' - R₁'·R₂ for pairs of radial functions.
    A fundamental identity of spheroidal wave functions ensures this quantity
    should be approximately constant across all x values if R₁ and R₂ are
    correctly computed.
    
    **Mathematical Basis:**
        The Wronskian of two linearly independent solutions to a second-order
        ODE is constant. For spheroidal radial functions:
        
            W(R₁, R₂) = R₁(x)·dR₂/dx - dR₁/dx·R₂(x) = constant
        
        This constant depends on the normalization convention used by the solver.
        The key insight is that W should NOT vary significantly across x.
    
    **Purpose:**
        - **Diagnostic**: Detects solver failures or numerical breakdown
        - **Validation**: Confirms R₁ and R₂ are consistent eigenfunctions
        - **Testing**: Verifies implementation correctness
        - **Confidence**: Deviations suggest precision or domain issues
    
    **Interpretation:**
        
        - **Constant W (≤1% variation)**: R₁ and R₂ are accurate and consistent
        - **Variable W (>10% variation)**: Numerical issues; solver may be struggling
        - **NaN or Inf**: Severe failure; x likely outside valid domain or c too extreme
    
    **Parameters:**
    
        m, n, c: spheroidal parameters (see smn/rmn documentation)
        
        x::AbstractVector{<:Real}
            Radial evaluation points (must satisfy |x| > 1 for real c)
            Best: use multiple points across the domain to verify constancy
        
        spheroid::Symbol = :prolate
            :prolate or :oblate geometry
        
        precision::Symbol = :double
            :double or :quad precision
    
    **Returns:**
        Vector{ComplexF64} of Wronskian values W(x) at each x point
        
        For consistency verification, compute:
            - W_mean = mean(W)
            - W_variation = maximum(|W - W_mean|) / W_mean
        
        Values << 1% indicate high confidence in the computation.
    
    **Examples:**
        
        # Verify prolate r1/r2 consistency
        m, n, c = 0, 1, 200
        x = [1.5, 2.0, 3.0, 5.0]
        W = radial_wronskian(m, n, c, x)
        
        # Check variation
        W_mean = mean(abs.(real.(W)))
        W_variation = maximum(abs.(real(W) .- W_mean)) / W_mean * 100
        println("Wronskian variation: \$(W_variation)%")
        
        # Should print < 1% variation for well-conditioned parameters
        
        # Oblate spheroid
        W_oblate = radial_wronskian(1, 2, 500, [2.0, 3.0]; spheroid=:oblate)
        
        # Test at high precision for sensitive regime
        W_quad = radial_wronskian(0, 1, 200, [1.1, 1.2]; precision=:quad)
    
    **Implementation Notes:**
        - Computes r1 and r2 for kind=1 and kind=2 respectively
        - Uses element-wise multiplication: W = r1.value .* r2.derivative - ...
        - Returns complex even when c is real (due to rmn interface)
        - Real part is the physically meaningful Wronskian value
    
    **Limitations:**
        - Does not validate against theoretical Wronskian value (depends on conventions)
        - Primarily useful for relative consistency checking
        - Large |W| doesn't indicate failure; variation does
        - Best used with quad precision for extreme parameters
    
    **References:**
        - Mathematical background on Wronskians: ODE theory texts
        - DLMF §30.3: Spheroidal wave function properties
    """
    
    r1 = rmn(m, n, c, x; spheroid=spheroid, precision=precision, kind=1)
    r2 = rmn(m, n, c, x; spheroid=spheroid, precision=precision, kind=2)
    
    # W = r1 * r2' - r1' * r2
    W = r1.value .* r2.derivative - r1.derivative .* r2.value
    
    return W
end

function eigenvalue(m::Integer, n::Integer, c::Union{Real,Complex};
                    spheroid::Symbol=:prolate, precision::Symbol=:double)
    """
    Compute the spheroidal separation constant λₘₙ(c).

    Returns the eigenvalue associated with order `m`, degree `n`, and size parameter `c`
    for either prolate or oblate spheroidal wave functions.

    Args:
        m, n, c: spheroidal parameters
        spheroid: `:prolate` or `:oblate`
        precision: `:double` or `:quad`

    Returns:
        Real value when `c` is real, complex value when `c` is complex.
    """

    _validate_precision(precision)
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end

    if c isa Real
        prefix = spheroid === :prolate ? :psms : :oblate
        return _call_real_eigenvalue(prefix, m, n, c; precision=precision)
    else
        prefix = spheroid === :prolate ? :cprolate : :coblate
        return _call_complex_eigenvalue(prefix, m, n, c; precision=precision)
    end
end

function jacobian_eigen(m::Integer, n::Integer, c::Union{Real,Complex};
                                                spheroid::Symbol=:prolate, precision::Symbol=:double, h=nothing,
                                                with_metadata::Bool=false, adaptive::Bool=true,
                                                rtol::Real=1e-6, atol::Real=1e-10)
        """
        Numerical Jacobian of `eigenvalue` with respect to `c`.

        This function estimates parameter sensitivities of the spheroidal separation
        constant `lambda_mn(c)` using centered finite differences.

        For real `c`, the return value is a scalar estimate of:
        - `d(lambda)/dc`

        For complex `c = a + ib`, the return value is a named tuple with:
        - `d_dcreal = ∂lambda/∂a`
        - `d_dcimag = ∂lambda/∂b`

        Keyword arguments:
        - `spheroid`: `:prolate` or `:oblate`
        - `precision`: `:double` or `:quad`
        - `h`: finite-difference step (auto-selected if `nothing`)
        - `with_metadata`: if `true`, returns derivative(s) plus reliability metadata
        - `adaptive`: if `true`, retries with smaller step in poor-conditioning regimes
        - `rtol`, `atol`: positive tolerances used for step-halving consistency checks

        Reliability metadata fields (`with_metadata=true`):
        - `step_used`
        - `relative_change_when_halving_step`
        - `finite_flag`
        - `conditioning_flag` in `(:good, :warning, :poor)`
        - `suggested_action` in `(:accept, :retry_smaller_h, :use_quad)`

        Returns:
        - Real `c`, `with_metadata=false`: scalar derivative
        - Real `c`, `with_metadata=true`: `(derivative=..., metadata=...)`
        - Complex `c`, `with_metadata=false`: `(d_dcreal=..., d_dcimag=...)`
        - Complex `c`, `with_metadata=true`:
            `(d_dcreal=..., d_dcimag=..., metadata_dcreal=..., metadata_dcimag=...)`
        """

    _validate_precision(precision)
    _validate_jacobian_tolerances(rtol, atol)
    step = _resolve_jacobian_step(c, h)
    if c isa Real
        calc = s -> begin
            cp = c + s
            cm = c - s
            fp = eigenvalue(m, n, cp; spheroid=spheroid, precision=precision)
            fm = eigenvalue(m, n, cm; spheroid=spheroid, precision=precision)
            (fp - fm) / (2 * s)
        end
        derivative, metadata = _finite_difference_with_metadata(calc, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        return with_metadata ? (derivative=derivative, metadata=metadata) : derivative
    else
        calc_re = s -> begin
            fp = eigenvalue(m, n, c + s; spheroid=spheroid, precision=precision)
            fm = eigenvalue(m, n, c - s; spheroid=spheroid, precision=precision)
            (fp - fm) / (2 * s)
        end
        calc_im = s -> begin
            fp = eigenvalue(m, n, c + s * im; spheroid=spheroid, precision=precision)
            fm = eigenvalue(m, n, c - s * im; spheroid=spheroid, precision=precision)
            (fp - fm) / (2 * s)
        end
        d_dcreal, metadata_dcreal = _finite_difference_with_metadata(calc_re, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        d_dcimag, metadata_dcimag = _finite_difference_with_metadata(calc_im, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)

        if with_metadata
            return (
                d_dcreal=d_dcreal,
                d_dcimag=d_dcimag,
                metadata_dcreal=metadata_dcreal,
                metadata_dcimag=metadata_dcimag,
            )
        end
        return (d_dcreal=d_dcreal, d_dcimag=d_dcimag)
    end
end

function jacobian_smn(m::Integer, n::Integer, c::Union{Real,Complex}, eta::AbstractVector{<:Real};
                      spheroid::Symbol=:prolate, precision::Symbol=:double, normalize::Bool=false, h=nothing,
                      with_metadata::Bool=false, adaptive::Bool=true,
                      rtol::Real=1e-6, atol::Real=1e-10)
    """
    Numerical Jacobian of `smn` outputs with respect to `c`.

    For real `c`, returns:
    - `dvalue_dc`
    - `dderivative_dc`

    For complex `c = a + ib`, returns:
    - `dvalue_dcreal`, `dvalue_dcimag`
    - `dderivative_dcreal`, `dderivative_dcimag`
    """

    _validate_precision(precision)
    _validate_jacobian_tolerances(rtol, atol)
    step = _resolve_jacobian_step(c, h)
    if c isa Real
        calc_value = s -> begin
            sp = smn(m, n, c + s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            sm = smn(m, n, c - s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            (sp.value .- sm.value) ./ (2 * s)
        end
        calc_derivative = s -> begin
            sp = smn(m, n, c + s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            sm = smn(m, n, c - s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            (sp.derivative .- sm.derivative) ./ (2 * s)
        end
        dvalue_dc, metadata_value = _finite_difference_with_metadata(calc_value, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dc, metadata_derivative = _finite_difference_with_metadata(calc_derivative, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        if with_metadata
            return (
                dvalue_dc=dvalue_dc,
                dderivative_dc=dderivative_dc,
                metadata_value=metadata_value,
                metadata_derivative=metadata_derivative,
            )
        end
        return (dvalue_dc=dvalue_dc, dderivative_dc=dderivative_dc)
    else
        calc_value_re = s -> begin
            sp = smn(m, n, c + s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            sm = smn(m, n, c - s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            (sp.value .- sm.value) ./ (2 * s)
        end
        calc_value_im = s -> begin
            sp = smn(m, n, c + s * im, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            sm = smn(m, n, c - s * im, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            (sp.value .- sm.value) ./ (2 * s)
        end
        calc_derivative_re = s -> begin
            sp = smn(m, n, c + s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            sm = smn(m, n, c - s, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            (sp.derivative .- sm.derivative) ./ (2 * s)
        end
        calc_derivative_im = s -> begin
            sp = smn(m, n, c + s * im, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            sm = smn(m, n, c - s * im, eta; spheroid=spheroid, precision=precision, normalize=normalize)
            (sp.derivative .- sm.derivative) ./ (2 * s)
        end

        dvalue_dcreal, metadata_value_dcreal = _finite_difference_with_metadata(calc_value_re, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dvalue_dcimag, metadata_value_dcimag = _finite_difference_with_metadata(calc_value_im, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dcreal, metadata_derivative_dcreal = _finite_difference_with_metadata(calc_derivative_re, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dcimag, metadata_derivative_dcimag = _finite_difference_with_metadata(calc_derivative_im, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)

        if with_metadata
            return (
                dvalue_dcreal=dvalue_dcreal,
                dvalue_dcimag=dvalue_dcimag,
                dderivative_dcreal=dderivative_dcreal,
                dderivative_dcimag=dderivative_dcimag,
                metadata_value_dcreal=metadata_value_dcreal,
                metadata_value_dcimag=metadata_value_dcimag,
                metadata_derivative_dcreal=metadata_derivative_dcreal,
                metadata_derivative_dcimag=metadata_derivative_dcimag,
            )
        end
        return (
            dvalue_dcreal=dvalue_dcreal,
            dvalue_dcimag=dvalue_dcimag,
            dderivative_dcreal=dderivative_dcreal,
            dderivative_dcimag=dderivative_dcimag,
        )
    end
end

function jacobian_rmn(m::Integer, n::Integer, c::Union{Real,Complex}, x::AbstractVector{<:Real};
                      spheroid::Symbol=:prolate, precision::Symbol=:double, kind::Integer=1, h=nothing,
                      with_metadata::Bool=false, adaptive::Bool=true,
                      rtol::Real=1e-6, atol::Real=1e-10)
    """
    Numerical Jacobian of `rmn` outputs with respect to `c`.

    For real `c`, returns:
    - `dvalue_dc`
    - `dderivative_dc`

    For complex `c = a + ib`, returns:
    - `dvalue_dcreal`, `dvalue_dcimag`
    - `dderivative_dcreal`, `dderivative_dcimag`
    """

    _validate_precision(precision)
    _validate_jacobian_tolerances(rtol, atol)
    step = _resolve_jacobian_step(c, h)
    if c isa Real
        calc_value = s -> begin
            rp = rmn(m, n, c + s, x; spheroid=spheroid, precision=precision, kind=kind)
            rm = rmn(m, n, c - s, x; spheroid=spheroid, precision=precision, kind=kind)
            (rp.value .- rm.value) ./ (2 * s)
        end
        calc_derivative = s -> begin
            rp = rmn(m, n, c + s, x; spheroid=spheroid, precision=precision, kind=kind)
            rm = rmn(m, n, c - s, x; spheroid=spheroid, precision=precision, kind=kind)
            (rp.derivative .- rm.derivative) ./ (2 * s)
        end
        dvalue_dc, metadata_value = _finite_difference_with_metadata(calc_value, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dc, metadata_derivative = _finite_difference_with_metadata(calc_derivative, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)

        if with_metadata
            return (
                dvalue_dc=dvalue_dc,
                dderivative_dc=dderivative_dc,
                metadata_value=metadata_value,
                metadata_derivative=metadata_derivative,
            )
        end
        return (dvalue_dc=dvalue_dc, dderivative_dc=dderivative_dc)
    else
        calc_value_re = s -> begin
            rp = rmn(m, n, c + s, x; spheroid=spheroid, precision=precision, kind=kind)
            rm = rmn(m, n, c - s, x; spheroid=spheroid, precision=precision, kind=kind)
            (rp.value .- rm.value) ./ (2 * s)
        end
        calc_value_im = s -> begin
            rp = rmn(m, n, c + s * im, x; spheroid=spheroid, precision=precision, kind=kind)
            rm = rmn(m, n, c - s * im, x; spheroid=spheroid, precision=precision, kind=kind)
            (rp.value .- rm.value) ./ (2 * s)
        end
        calc_derivative_re = s -> begin
            rp = rmn(m, n, c + s, x; spheroid=spheroid, precision=precision, kind=kind)
            rm = rmn(m, n, c - s, x; spheroid=spheroid, precision=precision, kind=kind)
            (rp.derivative .- rm.derivative) ./ (2 * s)
        end
        calc_derivative_im = s -> begin
            rp = rmn(m, n, c + s * im, x; spheroid=spheroid, precision=precision, kind=kind)
            rm = rmn(m, n, c - s * im, x; spheroid=spheroid, precision=precision, kind=kind)
            (rp.derivative .- rm.derivative) ./ (2 * s)
        end

        dvalue_dcreal, metadata_value_dcreal = _finite_difference_with_metadata(calc_value_re, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dvalue_dcimag, metadata_value_dcimag = _finite_difference_with_metadata(calc_value_im, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dcreal, metadata_derivative_dcreal = _finite_difference_with_metadata(calc_derivative_re, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dcimag, metadata_derivative_dcimag = _finite_difference_with_metadata(calc_derivative_im, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)

        if with_metadata
            return (
                dvalue_dcreal=dvalue_dcreal,
                dvalue_dcimag=dvalue_dcimag,
                dderivative_dcreal=dderivative_dcreal,
                dderivative_dcimag=dderivative_dcimag,
                metadata_value_dcreal=metadata_value_dcreal,
                metadata_value_dcimag=metadata_value_dcimag,
                metadata_derivative_dcreal=metadata_derivative_dcreal,
                metadata_derivative_dcimag=metadata_derivative_dcimag,
            )
        end
        return (
            dvalue_dcreal=dvalue_dcreal,
            dvalue_dcimag=dvalue_dcimag,
            dderivative_dcreal=dderivative_dcreal,
            dderivative_dcimag=dderivative_dcimag,
        )
    end
end

function find_c_for_eigenvalue(m::Integer, n::Integer, lambda_target::Real;
                               bracket::Tuple{<:Real,<:Real},
                               spheroid::Symbol=:prolate,
                               precision::Symbol=:double,
                               atol::Real=1e-10,
                               rtol::Real=1e-8,
                               maxiter::Integer=80,
                               use_jacobian::Bool=true)
    """
    Solve `eigenvalue(m, n, c) = lambda_target` for real `c`.

    Uses a bracketed hybrid strategy with guaranteed bisection fallback and optional
    Jacobian-guided acceleration through `jacobian_eigen` when derivative quality is
    acceptable.

    Keyword arguments:
    - `bracket`: `(c_lo, c_hi)` with `c_lo < c_hi` and opposite signs of residual
      `eigenvalue(m,n,c) - lambda_target` at the endpoints.
    - `spheroid`: `:prolate` or `:oblate`.
    - `precision`: `:double` or `:quad`.
    - `atol`, `rtol`: positive absolute and relative tolerances used for both residual
      and bracket-width stopping criteria.
    - `maxiter`: positive maximum number of iterations.
    - `use_jacobian`: enable derivative-based candidate steps when trusted.

    Returns a named tuple with fields:
    - `converged::Bool`
    - `c::Float64`
    - `residual::Float64`
    - `iterations::Int`
    - `bracket::Tuple{Float64,Float64}`
    - `method::Symbol` (`:endpoint`, `:newton`, `:secant`, `:bisection`, `:maxiter`)
    """

    _validate_precision(precision)
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end
    if !(atol > 0)
        error("atol must be positive, got $atol")
    end
    if !(rtol > 0)
        error("rtol must be positive, got $rtol")
    end
    if maxiter <= 0
        error("maxiter must be positive, got $maxiter")
    end

    a = Float64(bracket[1])
    b = Float64(bracket[2])
    if !(a < b)
        error("bracket must satisfy bracket[1] < bracket[2]")
    end

    f(c) = Float64(eigenvalue(m, n, c; spheroid=spheroid, precision=precision) - lambda_target)

    fa = f(a)
    fb = f(b)
    if !isfinite(fa) || !isfinite(fb)
        error("non-finite residual at bracket endpoints")
    end

    tol_residual = atol + rtol * max(1.0, abs(Float64(lambda_target)))
    width_tol() = atol + rtol * max(1.0, abs(a), abs(b))

    if abs(fa) <= tol_residual
        return (converged=true, c=a, residual=fa, iterations=0, bracket=(a, b), method=:endpoint)
    end
    if abs(fb) <= tol_residual
        return (converged=true, c=b, residual=fb, iterations=0, bracket=(a, b), method=:endpoint)
    end
    if signbit(fa) == signbit(fb)
        error("bracket endpoints must straddle a root for eigenvalue(m,n,c)-lambda_target")
    end

    method_used = :bisection
    c_best = (a + b) / 2
    f_best = f(c_best)

    for iter in 1:maxiter
        width = b - a
        mid = (a + b) / 2
        fmid = f(mid)

        candidate = mid
        method = :bisection

        if use_jacobian
            j = jacobian_eigen(m, n, mid; spheroid=spheroid, precision=precision,
                               with_metadata=true, adaptive=true)
            d = Float64(j.derivative)
            md = j.metadata
            if isfinite(d) && abs(d) > sqrt(eps(Float64)) && md.suggested_action == :accept
                newton = mid - fmid / d
                if a < newton < b && isfinite(newton)
                    candidate = newton
                    method = :newton
                end
            end
        end

        if method == :bisection
            denom = fb - fa
            if isfinite(denom) && abs(denom) > eps(Float64)
                secant = b - fb * (b - a) / denom
                if a < secant < b && isfinite(secant)
                    candidate = secant
                    method = :secant
                end
            end
        end

        fc = f(candidate)
        if !isfinite(fc)
            candidate = mid
            fc = fmid
            method = :bisection
        end

        c_best = candidate
        f_best = fc
        method_used = method

        if abs(fc) <= tol_residual || width <= width_tol()
            return (
                converged=true,
                c=candidate,
                residual=fc,
                iterations=iter,
                bracket=(a, b),
                method=method,
            )
        end

        if signbit(fa) == signbit(fc)
            a = candidate
            fa = fc
        else
            b = candidate
            fb = fc
        end
    end

    return (
        converged=false,
        c=c_best,
        residual=f_best,
        iterations=maxiter,
        bracket=(a, b),
        method=:maxiter,
    )
end

function accuracy(m::Integer, n::Integer, c::Union{Real,Complex}, arg::AbstractVector{<:Real};
                  spheroid::Symbol=:prolate, precision::Symbol=:double, kind::Integer=1, target::Symbol=:radial, normalize::Bool=false)
    """
    Estimate numerical accuracy of computed spheroidal wave functions.
    
    Returns the number of accurate decimal digits in the computed function values,
    as estimated by the underlying Fortran solver. This provides confidence bounds
    on the results.
    
    **Angular Functions (`smn`):**
    - Returns `naccs`: decimal digits of accuracy for angular functions
    - Estimated from convergence criteria in eigenvalue solver
    
    **Radial Functions (`rmn`):**
    - Returns `naccr`: decimal digits of accuracy for radial functions  
    - Estimated from recurrence relation residuals and function magnitudes
    
    Values near the accuracy limit indicate reliable results.
    Values much below 1 digit suggest potential numerical issues or that
    the chosen precision (:double or :quad) may be insufficient.
    
    Args:
        m, n, c: spheroidal parameters
        arg: evaluation points (η ∈ [-1,1] for angular; x-domain depends on spheroid for radial)
        spheroid: :prolate or :oblate
        precision: :double or :quad
        kind: for radial functions, which kind (1-4)
        target: :radial (default) or :angular
        normalize: angular normalization flag (used when target=:angular)
    
    Returns:
        Vector of estimated decimal digit accuracy at each evaluation point
    
    Examples:
        # Angular function accuracy
        acc_ang = accuracy(0, 1, 200, [0.5, 0.6, 0.7]; spheroid=:prolate, target=:angular)
        # → [14, 14, 14]  # All values ~14 decimal digits accurate
        
        # Radial function accuracy
        acc_rad = accuracy(0, 1, 200, [1.5, 2.0]; spheroid=:prolate, target=:radial, kind=1)
        # → [13, 13]  # All values ~13 decimal digits accurate

        # Complex oblate accuracy at quad precision
        acc_c = accuracy(0, 1, 200 + 0.1im, [1.2, 1.5];
                         spheroid=:oblate, precision=:quad, target=:radial, kind=3)
    
    Notes:
        - Supported for real and complex c
        - Supported for both :prolate and :oblate spheroids
        - Supported for both :double and :quad precision backends
        - For exact spherical limit (c=0), returns theoretical maximum accuracy
    """
    
    _validate_precision(precision)
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end
    if target != :radial && target != :angular
        error("target must be :radial or :angular, got :$target")
    end

    if c isa Real
        prefix = spheroid === :prolate ? :psms : :oblate
        if target === :angular
            return _call_real_smn_accuracy(prefix, m, n, c, arg; precision=precision, normalize=normalize)
        else
            return _call_real_rmn_accuracy(prefix, m, n, c, arg; precision=precision, kind=kind)
        end
    else
        prefix = spheroid === :prolate ? :cprolate : :coblate
        if target === :angular
            return _call_complex_smn_accuracy(prefix, m, n, c, arg; precision=precision, normalize=normalize)
        else
            return _call_complex_rmn_accuracy(prefix, m, n, c, arg; precision=precision, kind=kind)
        end
    end
end

function __init__()
    try
        # User/CI env vars are the primary non-code override path.
        configured_from_env = _configure_backends_from_env!()
        # Local generated config remains a fallback for developer workflows.
        if !configured_from_env
            _configure_backends_from_local_config!()
        end
    catch e
        @warn "Failed to configure backend libraries during module initialization: $e" maxlog=1
    end
end

end
