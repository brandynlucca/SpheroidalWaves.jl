module SpheroidalWaves

using Artifacts
using Libdl
using LinearAlgebra: Tridiagonal

export smn, rmn, radial_wronskian, accuracy, eigenvalue, eigenvalue_sweep, jacobian_eigen, jacobian_smn, jacobian_rmn, find_c_for_eigenvalue

const _backend_libraries = Dict{Symbol,Union{Nothing,String}}(
    :double => nothing,
    :quad => nothing,
)

const _backend_handles = Dict{String,Ptr{Cvoid}}()
const _backend_registry_lock = ReentrantLock()

const _ENV_BACKEND_DOUBLE = "SPHEROIDALWAVES_LIBRARY_DOUBLE"
const _ENV_BACKEND_QUAD = "SPHEROIDALWAVES_LIBRARY_QUAD"
const _ARTIFACT_DOUBLE = "spheroidal_backend_double"
const _ARTIFACT_QUAD = "spheroidal_backend_quad"

function _with_backend_registry_lock(f::F) where {F}
    lock(_backend_registry_lock)
    try
        return f()
    finally
        unlock(_backend_registry_lock)
    end
end

function _validate_precision(precision::Symbol)
    if precision != :double && precision != :quad
        error("precision must be :double or :quad, got :$precision")
    end
end

function set_backend_library!(path::AbstractString; precision::Symbol=:double)
    _validate_precision(precision)
    return _with_backend_registry_lock() do
        _backend_libraries[precision] = String(path)
    end
end

function backend_library(; precision::Symbol=:double)
    _validate_precision(precision)
    return _with_backend_registry_lock() do
        _backend_libraries[precision]
    end
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
    if haskey(ENV, _ENV_BACKEND_DOUBLE)
        configured_any |= _set_backend_from_candidate(ENV[_ENV_BACKEND_DOUBLE], :double, "ENV[$_ENV_BACKEND_DOUBLE]")
    end
    if haskey(ENV, _ENV_BACKEND_QUAD)
        configured_any |= _set_backend_from_candidate(ENV[_ENV_BACKEND_QUAD], :quad, "ENV[$_ENV_BACKEND_QUAD]")
    end
    return configured_any
end

function _backend_filename(stem::AbstractString)
    if Sys.iswindows()
        return "$stem.dll"
    elseif Sys.isapple()
        return "lib$stem.dylib"
    else
        return "lib$stem.so"
    end
end

function _configure_one_backend_from_artifact!(artifact_name::String, precision::Symbol)
    artifacts_toml = joinpath(dirname(@__FILE__), "..", "Artifacts.toml")
    if !isfile(artifacts_toml)
        return false
    end
    hash = Artifacts.artifact_hash(artifact_name, artifacts_toml)
    if hash === nothing
        return false
    end

    if !Artifacts.artifact_exists(hash)
        @warn "Backend artifact is bound but not installed" artifact=artifact_name precision maxlog=1
        return false
    end

    root = Artifacts.artifact_path(hash)
    stems = precision === :double ? ["spheroidal_batch_double"] : ["spheroidal_batch_quad"]
    candidates = String[]
    for stem in stems
        libname = _backend_filename(stem)
        append!(candidates, [
            joinpath(root, libname),
            joinpath(root, "lib", libname),
            joinpath(root, "bin", libname),
        ])
    end

    for path in candidates
        isfile(path) || continue
        if precision === :quad && !_has_required_quad_abi(path)
            continue  # A local source build can replace an older artifact.
        end
        precision === :double && !_has_scaled_abi(path) && continue
        if _set_backend_from_candidate(path, precision, "artifact $artifact_name")
            return true
        end
    end
    return false
end

function _configure_backends_from_artifacts!()
    configured_any = false
    configured_any |= _configure_one_backend_from_artifact!(_ARTIFACT_DOUBLE, :double)
    configured_any |= _configure_one_backend_from_artifact!(_ARTIFACT_QUAD, :quad)
    return configured_any
end

function _configure_backends_from_local_build!()
    build_root = normpath(joinpath(dirname(@__FILE__), "..", "build"))
    configured_any = false
    for precision in (:double, :quad)
        backend_library(precision=precision) === nothing || continue
        stem = "spheroidal_batch_$(precision)"
        filename = _backend_filename(stem)
        for path in (
            joinpath(build_root, "bin", filename),
            joinpath(build_root, "lib", filename),
            joinpath(build_root, filename),
        )
            if isfile(path)
                configured_any |= _set_backend_from_candidate(path, precision, "local build directory")
                break
            end
        end
    end
    return configured_any
end



function _require_backend_library(precision::Symbol)
    _validate_precision(precision)
    lib = _with_backend_registry_lock() do
        _backend_libraries[precision]
    end
    if lib === nothing
        error("""
        No backend library configured for precision :$precision.
        
        To resolve this, try one of:
        1. Ensure artifacts are available (they should download automatically on first use).
        2. Set environment variable: SPHEROIDALWAVES_LIBRARY_$(uppercase(String(precision))) = /path/to/lib
        3. Run: julia> import Pkg; Pkg.build("SpheroidalWaves")
        
        """)
    end
    return lib
end

function _require_backend_handle(lib::String)
    return _with_backend_registry_lock() do
        get!(_backend_handles, lib) do
            Libdl.dlopen(lib)
        end
    end
end

function _symbol_pointer(lib::String, symbol::Symbol)
    return Libdl.dlsym(_require_backend_handle(lib), symbol)
end

function _real_suffix(precision::Symbol)
    return precision === :double ? "_double" : "_quad"
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

function _combine_split_parts(hi::Vector{Float64}, lo::Vector{Float64})
    out = Vector{BigFloat}(undef, length(hi))
    @inbounds for i in eachindex(hi)
        out[i] = BigFloat(hi[i]) + BigFloat(lo[i])
    end
    return out
end

function _split_real_to_double_pair(x::Real)
    bx = BigFloat(x)
    hi = Float64(bx)
    lo = Float64(bx - BigFloat(hi))
    return hi, lo
end

function _split_real_vector_to_double_pairs(xs::AbstractVector{<:Real})
    n = length(xs)
    hi = Vector{Float64}(undef, n)
    lo = Vector{Float64}(undef, n)
    @inbounds for i in eachindex(xs)
        h, l = _split_real_to_double_pair(xs[i])
        hi[i] = h
        lo[i] = l
    end
    return hi, lo
end

function _combine_split_complex_parts(re_hi::Vector{Float64}, re_lo::Vector{Float64}, im_hi::Vector{Float64}, im_lo::Vector{Float64})
    out = Vector{Complex{BigFloat}}(undef, length(re_hi))
    @inbounds for i in eachindex(re_hi)
        re = BigFloat(re_hi[i]) + BigFloat(re_lo[i])
        im = BigFloat(im_hi[i]) + BigFloat(im_lo[i])
        out[i] = Complex{BigFloat}(re, im)
    end
    return out
end

const _QUAD_TEXT_WIDTH = 96

# Bound the decimal payload independently of the caller's BigFloat precision.
# 192 bits retains ample guard bits above the native binary128 significand.
_quad_input_text(x::Real) = string(BigFloat(x; precision=192))

function _encode_real_text_scalar(x::Real; width::Int=_QUAD_TEXT_WIDTH)
    s = _quad_input_text(x)
    if ncodeunits(s) > width
        error("Quad text payload overflow for scalar input; increase _QUAD_TEXT_WIDTH.")
    end
    out = fill(UInt8(' '), width)
    copyto!(out, 1, codeunits(s), 1, ncodeunits(s))
    return out
end

function _encode_real_text_vector(xs::AbstractVector{<:Real}; width::Int=_QUAD_TEXT_WIDTH)
    n = length(xs)
    out = fill(UInt8(' '), width * n)
    @inbounds for i in eachindex(xs)
        s = _quad_input_text(xs[i])
        if ncodeunits(s) > width
            error("Quad text payload overflow for vector input; increase _QUAD_TEXT_WIDTH.")
        end
        off = (i - 1) * width + 1
        copyto!(out, off, codeunits(s), 1, ncodeunits(s))
    end
    return out
end

function _decode_real_text_vector(buf::Vector{UInt8}, n::Integer; width::Int=_QUAD_TEXT_WIDTH)
    out = Vector{BigFloat}(undef, n)
    @inbounds for i in 1:n
        off = (i - 1) * width + 1
        s = strip(String(buf[off:(off + width - 1)]))
        out[i] = BigFloat(s)
    end
    return out
end

function _decode_scaled_real_text_vector(buf::Vector{UInt8}, exponents::Vector{Cint}, n::Integer; width::Int=_QUAD_TEXT_WIDTH)
    out = Vector{BigFloat}(undef, n)
    @inbounds for i in 1:n
        off = (i - 1) * width + 1
        s = strip(String(buf[off:(off + width - 1)]))
        mant = BigFloat(s)
        out[i] = iszero(exponents[i]) ? mant : mant * (BigFloat(10) ^ Int(exponents[i]))
    end
    return out
end

function _is_exact_spherical_limit(c::Real)
    return iszero(c)
end

function _is_exact_spherical_limit(c::Complex)
    return iszero(real(c)) && iszero(imag(c))
end

function _validate_radial_parameter(c)
    iszero(c) && throw(DomainError(c,
        "radial spheroidal functions require c != 0; the c=0 radial normalization is undefined"))
    return nothing
end

function _resolve_jacobian_step(c::Union{Real,Complex}, h, precision::Symbol=:double)
    T = precision === :quad ? BigFloat : Float64
    epsilon = precision === :quad ? T(2)^(-112) : eps(T)
    step = h === nothing ? cbrt(epsilon)*max(one(T),T(abs(c))) : T(h)
    isfinite(step) && step > 0 || error("h must be finite and positive, got $h")
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

function _jacobian_metadata(dh, dh2, step_used::Real; precision::Symbol, rtol::Real, atol::Real)
    finite_flag = _all_finite(dh) && _all_finite(dh2)
    rel_change = _max_relative_change(dh, dh2; atol=atol)
    conditioning_flag = !finite_flag ? :poor : (rel_change <= 10 * rtol ? :good : (rel_change <= 1000 * rtol ? :warning : :poor))
    suggested_action = _jacobian_suggested_action(conditioning_flag, finite_flag, precision)
    return (
        method=:finite_difference,
        step_used=step_used,
        relative_change_when_halving_step=rel_change,
        finite_flag=finite_flag,
        conditioning_flag=conditioning_flag,
        suggested_action=suggested_action,
    )
end

function _finite_difference_with_metadata(calc::Function, step::Real; precision::Symbol, adaptive::Bool, rtol::Real, atol::Real)
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

include("spherical_angular.jl")
include("angular_endpoints.jl")
include("complex_quad.jl")
include("angular_continuation.jl")
include("accuracy.jl")
include("diagnostics.jl")
include("scaled_values.jl")
include("angular_coefficients.jl")
include("angular_second_kind.jl")
include("parameter_derivatives.jl")
include("radial_parameter_derivatives.jl")
include("integral_eigenvalues.jl")
include("coordinate_zeros.jl")

function _call_real_smn(prefix::Symbol, m::Integer, n::Integer, c::Real, eta::AbstractVector{<:Real}; precision::Symbol=:double, normalize::Bool=false, with_accuracy::Bool=false)
    if _is_exact_spherical_limit(c)
        return _spherical_smn_real(m, n, eta; normalize, precision)
    end
    if _use_small_parameter_expansion(c)
        T = precision === :quad ? BigFloat : Float64
        spheroid = prefix === :psms ? :prolate : :oblate
        result = map(v -> T.(v),_small_parameter_smn(m,n,c,eta,spheroid,precision,normalize))
        return with_accuracy ? (;result...,accuracy=fill(-1,length(eta))) : result
    end

    lib = _require_backend_library(precision)
    symbol = if precision === :quad
        Symbol(String(prefix) * "_smn_batch_quad_text" * (with_accuracy ? "_acc" : ""))
    else
        suffix = _real_suffix(precision)
        Symbol(String(prefix) * "_smn_batch" * suffix)
    end
    fnptr = precision === :quad ? _quad_symbol_pointer(lib, symbol) : _symbol_pointer(lib, symbol)

    n_eta = Cint(length(eta))
    status = Ref{Cint}(0)

    if precision === :quad
        c_text = _encode_real_text_scalar(c)
        endpoint = _angular_endpoint_plan(m,n,c,eta,prefix === :psms ? :prolate : :oblate)
        eta_text = _encode_real_text_vector(endpoint === nothing ? eta : endpoint.native_points)
        value_text = fill(UInt8(' '), _QUAD_TEXT_WIDTH * Int(n_eta))
        derivative_text = fill(UInt8(' '), _QUAD_TEXT_WIDTH * Int(n_eta))
        value_exp = zeros(Cint, Int(n_eta))
        derivative_exp = zeros(Cint, Int(n_eta))

        estimate = fill(Cint(-1), length(eta))
        if with_accuracy
            ccall(fnptr, Cvoid,
                  (Cint, Cint, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}, Ptr{Cint}),
                  Cint(m), Cint(n), n_eta, _bool_to_cint(normalize), c_text, Cint(_QUAD_TEXT_WIDTH), eta_text,
                  value_text, value_exp, derivative_text, derivative_exp, status, estimate)
        else
            ccall(fnptr, Cvoid,
                  (Cint, Cint, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
                  Cint(m), Cint(n), n_eta, _bool_to_cint(normalize), c_text, Cint(_QUAD_TEXT_WIDTH), eta_text,
                  value_text, value_exp, derivative_text, derivative_exp, status)
        end

        _check_scalar_status(status[])
        value = _decode_scaled_real_text_vector(value_text, value_exp, Int(n_eta))
        derivative = _decode_scaled_real_text_vector(derivative_text, derivative_exp, Int(n_eta))
        _angular_endpoint_reconstruct!(value,derivative,endpoint,m,n)
        endpoint !== nothing && (estimate[endpoint.indices] .= -1)
        return with_accuracy ? (; value, derivative, accuracy=_reported_accuracy(estimate, value, precision)) : (; value, derivative)
    end

    eta64 = Float64.(eta)
    value = zeros(Float64, n_eta)
    derivative = zeros(Float64, n_eta)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}),
          Cint(m), Cint(n), Cdouble(c), n_eta, eta64, _bool_to_cint(normalize), value, derivative, status)

    _check_scalar_status(status[])
    return (; value, derivative)
end

function _call_real_rmn(prefix::Symbol, m::Integer, n::Integer, c::Real, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1, with_accuracy::Bool=false)
    # The native prolate kernel documents x=1 only for kind 1; calling its
    # second-kind path there leaves outputs uninitialized. Represent the
    # undefined real-prolate second/Hankel kinds consistently as NaN.
    if prefix === :psms && kind >= 2 && any(isone,x)
        T = precision === :quad ? BigFloat : Float64
        value = fill(complex(T(NaN),zero(T)),length(x))
        derivative = copy(value)
        estimate = fill(-1,length(x))
        interior = findall(z -> !isone(z),x)
        if !isempty(interior)
            result = _call_real_rmn(prefix,m,n,c,x[interior];precision,kind,with_accuracy)
            value[interior] = result.value
            derivative[interior] = result.derivative
            with_accuracy && (estimate[interior] = result.accuracy)
        end
        return with_accuracy ? (;value,derivative,accuracy=estimate) : (;value,derivative)
    end
    lib = _require_backend_library(precision)
    offset_input = precision === :quad && prefix === :psms
    symbol = if offset_input
        :psms_rmn_batch_quad_offset_acc
    elseif precision === :quad
        Symbol(String(prefix) * "_rmn_batch_quad_fullsplit" * (with_accuracy ? "_acc" : ""))
    else
        suffix = _real_suffix(precision)
        Symbol(String(prefix) * "_rmn_batch" * suffix)
    end
    fnptr = precision === :quad ? _quad_symbol_pointer(lib, symbol) : _symbol_pointer(lib, symbol)

    n_x = Cint(length(x))
    status = zeros(Cint, n_x)

    if precision === :quad
        c_hi, c_lo = _split_real_to_double_pair(c)
        # Preserve the small boundary distance before converting to the native
        # representation. Sending x first can lose most digits of x-1.
        coordinates = offset_input ? BigFloat.(x) .- 1 : x
        x_hi, x_lo = _split_real_vector_to_double_pairs(coordinates)
        value_re_hi = zeros(Float64, n_x)
        value_re_lo = zeros(Float64, n_x)
        value_im_hi = zeros(Float64, n_x)
        value_im_lo = zeros(Float64, n_x)
        deriv_re_hi = zeros(Float64, n_x)
        deriv_re_lo = zeros(Float64, n_x)
        deriv_im_hi = zeros(Float64, n_x)
        deriv_im_lo = zeros(Float64, n_x)

        estimate = fill(Cint(-1), length(x))
        if with_accuracy || offset_input
            ccall(fnptr, Cvoid,
                  (Cint, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}),
                  Cint(m), Cint(n), Cdouble(c_hi), Cdouble(c_lo), n_x, x_hi, x_lo, Cint(kind),
                  value_re_hi, value_re_lo, value_im_hi, value_im_lo,
                  deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, estimate)
        else
            ccall(fnptr, Cvoid,
                  (Cint, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
                  Cint(m), Cint(n), Cdouble(c_hi), Cdouble(c_lo), n_x, x_hi, x_lo, Cint(kind),
                  value_re_hi, value_re_lo, value_im_hi, value_im_lo,
                  deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status)
        end

        _check_vector_status(status)
        value = _combine_split_complex_parts(value_re_hi, value_re_lo, value_im_hi, value_im_lo)
        derivative = _combine_split_complex_parts(deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo)
        return with_accuracy ? (; value, derivative, accuracy=_reported_accuracy(estimate, value, precision)) : (; value, derivative)
    end

    x64 = Float64.(x)
    value_re = zeros(Float64, n_x)
    value_im = zeros(Float64, n_x)
    deriv_re = zeros(Float64, n_x)
    deriv_im = zeros(Float64, n_x)

    ccall(fnptr, Cvoid,
          (Cint, Cint, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
          Cint(m), Cint(n), Cdouble(c), n_x, x64, Cint(kind), value_re, value_im, deriv_re, deriv_im, status)

    _check_vector_status(status)
    value = _complex_parts(value_re, value_im)
    derivative = _complex_parts(deriv_re, deriv_im)
    return (; value, derivative)
end

function _call_complex_smn_raw(prefix::Symbol, m::Integer, n::Integer, c::Complex, eta::AbstractVector{<:Real}; precision::Symbol=:double, normalize::Bool=false)
    if _is_exact_spherical_limit(c)
        return _spherical_smn_complex(m, n, eta; normalize, precision)
    end
    if _use_small_parameter_expansion(c)
        T = precision === :quad ? BigFloat : Float64
        spheroid = prefix === :cprolate ? :prolate : :oblate
        result = _small_parameter_smn(m,n,c,eta,spheroid,precision,normalize)
        return map(v -> Complex{T}.(v),result)
    end
    if iszero(real(c))
        opposite=prefix===:cprolate ? :oblate : :psms
        result=_call_real_smn(opposite,m,n,abs(imag(c)),eta;precision,normalize)
        return (;value=complex.(result.value),derivative=complex.(result.derivative))
    end

    points = precision === :quad ? BigFloat.(eta) : Float64.(eta)
    if precision === :quad
        result = _call_complex_quad(prefix, m, n, c, points, 1, _bool_to_cint(normalize))
        d = result.data
        value = complex.(d[1, :], d[2, :])
        derivative = complex.(d[3, :], d[4, :])
    else
        lib = _require_backend_library(precision)
        fnptr = _symbol_pointer(lib, Symbol(String(prefix) * "_smn_batch_c8"))
        n_eta = Cint(length(points))
        value_re = zeros(Float64, n_eta)
        value_im = zeros(Float64, n_eta)
        deriv_re = zeros(Float64, n_eta)
        deriv_im = zeros(Float64, n_eta)
        status = Ref{Cint}(0)
        ccall(fnptr, Cvoid,
              (Cint, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}),
              Cint(m), Cint(n), Cdouble(real(c)), Cdouble(imag(c)), n_eta, points, _bool_to_cint(normalize),
              value_re, value_im, deriv_re, deriv_im, status)
        _check_scalar_status(status[])
        value = _complex_parts(value_re, value_im)
        derivative = _complex_parts(deriv_re, deriv_im)
    end
    return (; value, derivative)
end

function _call_complex_rmn_raw(prefix::Symbol, m::Integer, n::Integer, c::Complex, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1)
    if _radial_needs_analytic(c,x,kind)
        spheroid=prefix===:cprolate ? :prolate : :oblate
        seed=_call_complex_eigenvalue(prefix,m,n,c;precision)
        return _radial_analytic_values(m,n,c,x,spheroid,precision,kind;eigenvalue_seed=seed)
    end
    precision === :quad && return _complex_quad_radial(prefix, m, n, c, x, kind)
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
    if precision === :quad
        return _call_real_smn(prefix,m,n,c,eta;precision,normalize,with_accuracy=true).accuracy
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
    return _reported_accuracy(naccs, value, precision)
end

function _call_real_rmn_accuracy(prefix::Symbol, m::Integer, n::Integer, c::Real, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1)
    if precision === :quad
        return _call_real_rmn(prefix,m,n,c,x;precision,kind,with_accuracy=true).accuracy
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
    return _reported_accuracy(naccr, complex.(value_re, value_im), precision)
end

function _call_complex_smn_accuracy(prefix::Symbol, m::Integer, n::Integer, c::Complex, eta::AbstractVector{<:Real}; precision::Symbol=:double, normalize::Bool=false)
    if precision === :quad
        result = _call_complex_quad(prefix,m,n,c,eta,1,_bool_to_cint(normalize))
        values = complex.(result.data[1,:], result.data[2,:])
        return _reported_accuracy(result.accuracy, values, precision)
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
    return _reported_accuracy(naccs, complex.(value_re, value_im), precision)
end

function _call_complex_rmn_accuracy(prefix::Symbol, m::Integer, n::Integer, c::Complex, x::AbstractVector{<:Real}; precision::Symbol=:double, kind::Integer=1)
    if precision === :quad
        result = _call_complex_quad(prefix,m,n,c,x,2,kind)
        values = complex.(result.data[5,:], result.data[6,:])
        return _reported_accuracy(result.accuracy, values, precision)
    end

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
    return _reported_accuracy(naccr, complex.(value_re, value_im), precision)
end

function _call_real_eigenvalue(prefix::Symbol, m::Integer, n::Integer, c::Real; precision::Symbol=:double)
    if _is_exact_spherical_limit(c)
        return precision === :quad ? BigFloat(n * (n + 1)) : Float64(n * (n + 1))
    end
    if _use_small_parameter_expansion(c)
        T = precision===:quad ? BigFloat : Float64
        spheroid = prefix===:oblate ? :oblate : :prolate
        plan = _small_parameter_plan(m,n,c,spheroid,precision)
        return T(plan.lambda)
    end

    lib = _require_backend_library(precision)
    symbol = if precision === :quad
        Symbol(String(prefix) * "_eigenvalue_quad_fullsplit")
    else
        suffix = _real_suffix(precision)
        Symbol(String(prefix) * "_eigenvalue" * suffix)
    end
    fnptr = _symbol_pointer(lib, symbol)

    if precision === :quad
        c_hi, c_lo = _split_real_to_double_pair(c)
        eig_hi = Ref{Cdouble}(0.0)
        eig_lo = Ref{Cdouble}(0.0)
        status = Ref{Cint}(0)

        ccall(fnptr, Cvoid,
              (Cint, Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}),
              Cint(m), Cint(n), Cdouble(c_hi), Cdouble(c_lo), eig_hi, eig_lo, status)

        _check_scalar_status(status[])
        return BigFloat(eig_hi[]) + BigFloat(eig_lo[])
    end

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
        T = precision === :quad ? BigFloat : Float64
        return complex(T(n) * (T(n) + 1), zero(T))
    end
    if _use_small_parameter_expansion(c)
        T = precision === :quad ? BigFloat : Float64
        spheroid = prefix === :cprolate ? :prolate : :oblate
        return Complex{T}(_small_parameter_plan(m,n,c,spheroid,precision).lambda)
    end
    if iszero(real(c))
        opposite=prefix===:cprolate ? :oblate : :psms
        return complex(_call_real_eigenvalue(opposite,m,n,abs(imag(c));precision))
    end

    if precision === :quad
        result = _call_complex_quad(prefix, m, n, c, [0], 0, 0)
        return complex(result.data[1, 1], result.data[2, 1])
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

"""
Spheroidal angular wave functions of the first or second kind.

Computes prolate or oblate angular functions and their coordinate derivatives.
`kind=1` (default) selects Ps; `kind=2` selects the Ferrers-branch Qs of
DLMF 30.5 and 30.8(ii), on the real interval -1 < η < 1.

The first-kind functions are eigenfunction solutions to the angular wave equation,
orthogonal on the interval η ∈ [-1, 1]. They closely resemble associated Legendre
polynomials when c is small.

**Mathematical Background:**
- The spheroidal wave equation separates into radial and angular parts
- Prolate angular equation: d/dη[(1-η²)dS/dη] + [λ - c²η² - m²/(1-η²)]S = 0
- For oblate angular functions, replace -c²η² by +c²η²
- Here λ(c,m,n) is the eigenvalue (separation constant)
- Angular values and derivatives include the Condon–Shortley phase (-1)^m
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
    η::Union{Real,AbstractVector{<:Real}}
        Evaluation point(s) in [-1, 1]
        Represents cos(θ) where θ is angle in spherical coordinates
        Scalar inputs are forwarded as a one-point batch
    
    Keyword Arguments:
        spheroid::Symbol = :prolate
            :prolate or :oblate geometry
            - :prolate: elongated spheroid (ζ > 0 in ξ,η,ϕ coordinates)
            - :oblate: flattened spheroid (ζ < 0)
        
        precision::Symbol = :double
            :double (default) or :quad precision
            - :double: Float64 outputs (about 16 significant decimal digits)
            - :quad: BigFloat outputs, targeting quad precision (slower)
        
        normalize::Bool = false
            Whether to scale by sqrt((2n+1)(n-m)! / (2(n+m)!))
            - false: Meixner-Schafke normalization (norm → ∞ as m → ∞)
            - true: Unity normalization (constant norm, convenient for applications)
            - Only available for kind=1; kind=2 rejects normalize=true

        kind::Integer = 1
            1: regular first-kind Ps, with the existing normalization
            2: second-kind Qs, normalized by DLMF 30.5.4 and 30.11.4
               Qs has parity (-1)^(n-m+1), and equals Ferrers Q_n^m at c=0.
               Endpoint values and derivatives are signed, one-sided infinities.
               Complex c uses the same eigenmode branch as kind=1.

        scaled::Bool = false
            Return value and derivative as (mantissa, exponent) arrays, with
            value = mantissa * 10^exponent, preserving native scaling.
        logderivative::Bool = false
            Append S′/S, calculated before output conversion; NaN at a computed zero.
        second_derivative::Bool = false
            Append S″ from the ODE, using one-sided limits at singular endpoints.

**Returns:**
    NamedTuple with fields:
    - .value::Vector{Real or Complex}
        Function values Smₙ(η, c) at each η point
    - .derivative::Vector{Real or Complex}
        First derivatives dSmₙ/dη at each η point

**Accuracy:**
    Accuracy depends on parameter, order and coordinate. The precision setting
    selects the working/output precision, not a guaranteed number of accurate digits.

**Special Cases:**
    - c = 0: Returns Ferrers P_n^m (kind=1) or Q_n^m (kind=2)
    - kind=1, c = 0, η = ±1: Uses one-sided endpoint derivatives; m = 1 gives signed infinities
    - c → ∞: Functions oscillate rapidly; use appropriate η resolution
    - Parity in η is (-1)^(n-m) for kind=1 and (-1)^(n-m+1) for kind=2
    - Qs is logarithmically singular at ±1 for m=0 and algebraically singular
      for m>0. Its logarithmic derivative has limits -Inf at -1 and Inf at +1.
      Singular or unresolved second-kind normalization raises an error.

**Performance:**
    - First-kind batches use the native Fortran backend.
    - Second-kind batches reuse coefficient data and Taylor propagation across
      sorted coordinates. They are slower, particularly near singular endpoints.

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

    # Second kind, with quad precision
    qs = smn(1, 2, 1.25, 0.3; kind=2, precision=:quad)

**Notes:**
    - Order and degree must satisfy 0 ≤ m ≤ n
    - All η values must be in [-1, 1]
    - With kind=1 and normalize=false, uses Meixner–Schäfke normalization
    - Angular values and derivatives include the Condon–Shortley phase (-1)^m
    - Complex c continues the eigenmode and normalization sign vertically from
      real(c), matching neighboring native labels of the same parity when needed.
      Unresolved continuation paths raise an error. See `eigenvalue_sweep` for
      eigenvalue continuation along an explicit complex path.

**References:**
    - DLMF §30.2: https://dlmf.nist.gov/30.2
    - DLMF §30.4(i), phase after Eq. 30.4.1 and exact-zero identity Eq. 30.4.2:
      https://dlmf.nist.gov/30.4#i
    - Second-kind definition and normalization: https://dlmf.nist.gov/30.5
      and https://dlmf.nist.gov/30.8#ii
    - Van Buren & Boisvert (2004): Accurate calculation of prolate spheroidal wave functions
"""
function smn(m::Integer, n::Integer, c::Union{Real,Complex}, eta::Real;
             kwargs...)
    return smn(m, n, c, [eta]; kwargs...)
end

function smn(m::Integer, n::Integer, c::Union{Real,Complex}, eta::AbstractVector{<:Real};
             spheroid::Symbol=:prolate, precision::Symbol=:double, normalize::Bool=false,
             kind::Integer=1, second_derivative::Bool=false, scaled::Bool=false, logderivative::Bool=false)

    _validate_wave_arguments(m,n,c,eta,spheroid,precision,:angular;kind)
    if c isa Complex && iszero(real(c)) && !iszero(c)
        opposite=spheroid===:prolate ? :oblate : :prolate
        r=smn(m,n,abs(imag(c)),eta;spheroid=opposite,precision,normalize,kind,scaled,logderivative,second_derivative)
        return map(v->v isa NamedTuple ? (;mantissa=complex.(v.mantissa),v.exponent) : complex.(v),r)
    end
    if kind == 2
        return _angular_second_kind(m,n,c,eta,spheroid,precision,normalize,scaled,logderivative,second_derivative)
    end
    if _use_angular_expansion(n,c,spheroid)
        return _coefficient_smn(m,n,c,eta,spheroid,precision,normalize,scaled,logderivative,second_derivative)
    end
    if scaled || logderivative
        return _extended_wave(m,n,c,eta,spheroid,precision,:angular,Int(normalize),scaled,logderivative,second_derivative)
    end

    result = if c isa Real
        prefix = spheroid === :prolate ? :psms : :oblate
        _call_real_smn(prefix, m, n, c, eta; precision=precision, normalize=normalize)
    else
        prefix = spheroid === :prolate ? :cprolate : :coblate
        _call_complex_smn(prefix, m, n, c, eta; precision=precision, normalize=normalize)
    end
    result = _angular_phase!(result, m)
    _fix_regular_endpoints!(result,m,n,c,eta,spheroid,precision,:angular,Int(normalize))
    return second_derivative ? _wave_second_derivative(result,m,n,c,eta,spheroid,precision,:angular;option=Int(normalize)) : result
end

"""
Spheroidal radial wave functions and coordinate derivatives.

Computes prolate or oblate spheroidal radial wave functions and their derivatives.
Ordinary calls reconstruct the values. `scaled=true` instead preserves native
mantissas and decimal exponents to avoid reconstruction overflow and underflow.

These are eigenfunction solutions to the radial part of the spheroidal wave equation,
defined for x >= 1 (prolate) or x >= 0 (oblate). Complex prolate calls require
x > 1. At the real prolate endpoint, only the first kind is defined.
Four kinds of radial functions are available.

**Mathematical Background:**
- Prolate radial functions satisfy: d/dξ[(ξ²-1)dRmₙ/dξ] + [c²ξ² - λ - m²/(ξ²-1)]Rmₙ = 0
- λ(c,m,n) is the eigenvalue from the angular function problem
- Functions are related to spherical Bessel and Hankel functions when c→0
- Characteristic-exponent representation: R = R_characteristic × 10^exponent

**Parameters:**

    m, n, c: See smn() documentation
    
    x::Union{Real,AbstractVector{<:Real}}
        Radial evaluation points
        - For prolate: x >= 1
        - For oblate: x >= 0
        - Vectorized computation: all points evaluated in single Fortran call
    
    kind::Integer = 1
        Which radial function kind:
        
        1: First kind R₁ (propagating, finite at x→∞)
           - Regular solution; used for most applications
           - Analogous to spherical Bessel j_ℓ(kr)
        
        2: Second kind R₂
           - Complementary solution
           - Analogous to spherical Bessel y_ℓ(kr)
        
        3: Third kind (linear combination) = R₁ + i·R₂
           - Hankel function analog (outgoing wave)
        
        4: Fourth kind (linear combination) = R₁ - i·R₂
           - Alternate Hankel function analog

**Output options:**
    - `scaled=true`: return `(mantissa, exponent)` arrays for values and derivatives,
      representing `mantissa*10^exponent`.
    - `logderivative=true`: append `R′/R`, calculated before output conversion;
      return NaN at a computed zero.
    - `second_derivative=true`: append R″ from the ODE, with one-sided limits for
      the regular prolate first kind at x=1; undefined kinds retain NaN.

**Returns (default options):**
    NamedTuple with fields:
    - .value::Vector{ComplexF64 or Complex{BigFloat}}
        Function values Rmₙ(x, c) (always complex, even for real input)
        Real parameters give real kinds 1/2 and complex combinations for kinds 3/4.
        Complex parameters can give complex values for every kind.
    
    - .derivative::Vector{ComplexF64 or Complex{BigFloat}}
        First derivatives dRmₙ/dx

**Accuracy & Overflow Protection:**
    - Complex kinds 3/4 use direct spherical Hankel expansions, preserving the
      outgoing/incoming exponential and its coordinate derivatives
    - Uses characteristic-exponent representation internally
    - Ordinary output can overflow or underflow during reconstruction
    - `scaled=true` transfers mantissas and exponents before reconstruction
    - Scaling does not improve the native solver's numerical accuracy

**Special Cases:**
    - c = 0: Throws DomainError; no zero-parameter radial normalization is defined
    - x = 1 (real prolate boundary): kinds 2–4 return NaN
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
    
    # Verification: scaled prolate Wronskian should be one (c=200 here)
    W = r1.value .* r2.derivative - r1.derivative .* r2.value
    println(200 .* (x.^2 .- 1) .* W)  # Approximately one at each point
    
    # Oblate spheroid
    r_oblate = rmn(1, 2, 500, x; spheroid=:oblate, kind=1)

**Notes:**
    - Return values are complex even when c and x are real
    - For real c, kind=1 has real values; kind=2 has real values; kinds 3,4 are complex
    - Scalar coordinates return length-one arrays
    - For real c > 0: W = 1/(c*(x²-1)) (prolate), 1/(c*(x²+1)) (oblate)

**Performance:**
    - Vectorized batch computation in Fortran
    - Four kinds computed efficiently in single Fortran call
    - Typical: 1000 points in ~0.001s

**References:**
    - DLMF §30.3: https://dlmf.nist.gov/30.3
    - Van Buren & Boisvert (2004): Accurate calculation of prolate spheroidal wave functions
"""
function rmn(m::Integer, n::Integer, c::Union{Real,Complex}, x::AbstractVector{<:Real};
             spheroid::Symbol=:prolate, precision::Symbol=:double, kind::Integer=1,
             second_derivative::Bool=false, scaled::Bool=false, logderivative::Bool=false)

    _validate_wave_arguments(m,n,c,x,spheroid,precision,:radial;kind)
    if _radial_needs_analytic(c,x,kind)
        return _radial_analytic_wave(m,n,c,x,spheroid,precision,kind,scaled,logderivative,second_derivative)
    end
    if scaled || logderivative
        return _extended_wave(m,n,c,x,spheroid,precision,:radial,kind,scaled,logderivative,second_derivative)
    end

    result = if c isa Real
        prefix = spheroid === :prolate ? :psms : :oblate
        _call_real_rmn(prefix, m, n, c, x; precision=precision, kind=kind)
    else
        prefix = spheroid === :prolate ? :cprolate : :coblate
        _call_complex_rmn(prefix, m, n, c, x; precision=precision, kind=kind)
    end
    _fix_regular_endpoints!(result,m,n,c,x,spheroid,precision,:radial,kind)
    return second_derivative ? _wave_second_derivative(result,m,n,c,x,spheroid,precision,:radial;option=kind) : result
end

rmn(m::Integer,n::Integer,c::Union{Real,Complex},x::Real;kwargs...) = rmn(m,n,c,[x];kwargs...)

"""
    radial_wronskian(m, n, c, x; spheroid=:prolate, precision=:double, form=:raw)

Return `W = R1 .* R2′ - R1′ .* R2` at the radial coordinates `x`.
The result uses the same complex element type as `rmn`.

For real `c > 0`, the standard normalization gives `W(x) = 1/(c*(x^2-1))`
for prolate functions and `W(x) = 1/(c*(x^2+1))` for oblate functions.
The raw Wronskian therefore varies with position. At `c=0`, this function
throws `DomainError`, as does `rmn`.

Use the scaled Wronskian to check consistency:

```julia
c = 200.0
x = [1.5, 2.0, 3.0, 5.0]
W = radial_wronskian(0, 1, c, x)
scaled_W = c .* (x.^2 .- 1) .* W
maximum(abs.(scaled_W .- 1))  # Small for an accurate prolate solution pair
```

For oblate functions, replace `x.^2 .- 1` by `x.^2 .+ 1`.
This is a consistency check, not an independent proof of accuracy: correlated
errors in the two functions can preserve the Wronskian. Precision loss can also
occur when its two product terms nearly cancel. `form=:raw` is the default;
`form=:normalized` returns `c*(x^2∓1)*W` and `form=:error` returns its absolute
difference from one. Normalized forms keep native scaling during intermediate
products. No form supplies a rigorous accuracy bound.

"""
function radial_wronskian(m::Integer, n::Integer, c::Union{Real,Complex}, x::AbstractVector{<:Real};
                          spheroid::Symbol=:prolate, precision::Symbol=:double, form::Symbol=:raw)

    form in (:raw,:normalized,:error) || throw(ArgumentError("form must be :raw, :normalized or :error"))

    _validate_wave_arguments(m,n,c,x,spheroid,precision,:radial)
    r1 = form === :raw ? rmn(m,n,c,x;spheroid,precision,kind=1) :
                        _scaled_native_values(m,n,c,x,spheroid,precision,:radial,1)
    r2 = form === :raw ? rmn(m,n,c,x;spheroid,precision,kind=2) :
                        _scaled_native_values(m,n,c,x,spheroid,precision,:radial,2)
    
    # W = r1 * r2' - r1' * r2
    W = r1.value .* r2.derivative - r1.derivative .* r2.value
    
    form === :raw && return W
    T = precision === :quad ? BigFloat : Float64
    xs = BigFloat.(x)
    factor = spheroid === :prolate ? (xs .- 1).*(xs .+ 1) : xs.^2 .+ 1
    normalized = c .* factor .* W
    return form === :normalized ? Complex{T}.(normalized) : T.(abs.(normalized .- 1))
end

radial_wronskian(m::Integer,n::Integer,c::Union{Real,Complex},x::Real;kwargs...) =
    radial_wronskian(m,n,c,[x];kwargs...)

"""
    eigenvalue(m, n, c; spheroid=:prolate, precision=:double,
               operator=:separation, form=:value)

Compute the spheroidal separation constant λₘₙ(c), or an integral-operator eigenvalue.

Returns the eigenvalue associated with order `m`, degree `n`, and size parameter `c`
for either prolate or oblate spheroidal wave functions.

`operator=:separation` is the default. The integral operators require `m=0`,
`n>=0`, `spheroid=:prolate`, and finite real `c>=0`:

- `operator=:concentration`: eigenvalue Λ of the kernel
  `sin(c*(x-t))/(pi*(x-t))` integrated over `[-1,1]` (diagonal value `c/pi`).
  `form=:value` returns Λ, `:log` returns `log(Λ)`, and `:complement` returns
  `1-Λ` computed before output rounding. For positive c, `0<Λ<1`.
- `operator=:fourier`: eigenvalue μ of the kernel `exp(-im*c*x*t)` integrated
  over `[-1,1]`. Its phase is `(-im)^n`, and `Λ=c*abs2(μ)/(2pi)` for c>0.
  Only `form=:value` is supported. Returns `ComplexF64` or `Complex{BigFloat}`.

At c=0, concentration eigenvalues are zero; Fourier eigenvalues are 2 for
n=0 and zero otherwise. Concentration results are `Float64` or `BigFloat`.
Ordinary output may round to zero or one; use the logarithm or complement
forms to retain small tails. Working precision and expansion size are refined
internally; failed convergence raises an error. Definitions: DLMF 30.15.3, 30.15.5.

For the separation constant with complex `c`, follow degree `n` from `real(c)` vertically to `c`, matching
the same mode used by `smn` and `rmn`. Use `eigenvalue_sweep` with a complex
vector for continuation along another path; closed paths can exchange modes.

Args:
    m, n, c: spheroidal parameters
    spheroid: `:prolate` or `:oblate`
    precision: `:double` or `:quad`

Returns:
    For the separation constant: real when `c` is real, complex when `c` is complex.
    Integral-operator return types are specified above.
"""
function eigenvalue(m::Integer, n::Integer, c::Union{Real,Complex};
                    spheroid::Symbol=:prolate, precision::Symbol=:double,
                    operator::Symbol=:separation, form::Symbol=:value)

    _validate_precision(precision)
    operator in (:separation,:concentration,:fourier) ||
        throw(ArgumentError("operator must be :separation, :concentration or :fourier"))
    operator!==:separation && return _integral_eigenvalue(m,n,c,spheroid,precision,operator,form)
    form===:value || throw(ArgumentError("the separation constant supports only form=:value"))
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end

    if c isa Real
        prefix = spheroid === :prolate ? :psms : :oblate
        return _call_real_eigenvalue(prefix, m, n, c; precision=precision)
    else
        prefix = spheroid === :prolate ? :cprolate : :coblate
        0 <= m <= n || error("require 0 <= m <= n")
        iszero(real(c)) && return _call_complex_eigenvalue(prefix,m,n,c;precision)
        return _complex_mode_state(prefix,m,n,c,precision).lambda
    end
end

include("eigenvalue_sweep.jl")

"""
    Numerical Jacobian of `eigenvalue` with respect to `c`.

    By default, differentiate the refined Legendre coefficient eigenproblem.
    Supplying `h` explicitly selects centered finite differences for comparison.

    `operator=:concentration` or `:fourier` differentiates the corresponding
    integral eigenvalue for real prolate `m=0`, `c>=0`. Concentration also
    accepts `form=:log` or `:complement`: these differentiate `log(Λ)` or `1-Λ`.
    The default analytic calculation retains small tails before output rounding.
    Integral results are `Float64`/`BigFloat` for concentration and
    `ComplexF64`/`Complex{BigFloat}` for Fourier, selected by `precision`.
    At zero, return right-hand limits: Λ₀′=2/pi, Λₙ′=0 for n>0;
    μ₁′=-2im/3 and all other μₙ′=0. The logarithmic derivative limit is `Inf`.
    Integral metadata uses `method=:integral_identity` or `:right_limit`;
    a divergent limit has `conditioning_flag=:singular`.
    An explicit `h` must be smaller than positive `c`; at zero it uses a
    forward stencil. The singular log derivative at zero requires `h=nothing`.

    For real `c`, the return value is a scalar estimate of:
    - `d(lambda)/dc`

    For complex `c = a + ib`, the return value is a named tuple with:
    - `d_dcreal = ∂lambda/∂a`
    - `d_dcimag = ∂lambda/∂b`

    Keyword arguments:
    - `spheroid`: `:prolate` or `:oblate`
    - `precision`: `:double` or `:quad`
    - `h`: explicit finite-difference step; `nothing` uses coefficient sensitivities
    - `with_metadata`: if `true`, returns derivative(s) plus reliability metadata
    - `adaptive`: if `true`, retries with smaller step in poor-conditioning regimes
    - `rtol`, `atol`: positive tolerances used for step-halving consistency checks

    Reliability metadata (`with_metadata=true`) identifies the `method`.
    Coefficient sensitivities report expansion refinement, tail size, and
    eigenproblem/sensitivity residuals. Their `step_used` and
    `relative_change_when_halving_step` are `nothing`. Finite differences report:
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
function jacobian_eigen(m::Integer, n::Integer, c::Union{Real,Complex};
                                                spheroid::Symbol=:prolate, precision::Symbol=:double, h=nothing,
                                                with_metadata::Bool=false, adaptive::Bool=true,
                                                rtol::Real=1e-6, atol::Real=1e-10,
                                                operator::Symbol=:separation, form::Symbol=:value)

    _validate_precision(precision)
    _validate_jacobian_tolerances(rtol, atol)
    operator in (:separation,:concentration,:fourier) ||
        throw(ArgumentError("operator must be :separation, :concentration or :fourier"))
    operator!==:separation && return _integral_eigen_jacobian(m,n,c,spheroid,precision,
        operator,form,h,with_metadata,adaptive,rtol,atol)
    form===:value || throw(ArgumentError("the separation constant supports only form=:value"))
    h === nothing && return _coefficient_eigen_jacobian(m,n,c,spheroid,precision,with_metadata)
    step = _resolve_jacobian_step(c, h, precision)
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
        evaluate = _complex_local_evaluator(m,n,c,spheroid,precision)
        calc_re = s -> begin
            fp = evaluate(c+s)
            fm = evaluate(c-s)
            (fp - fm) / (2 * s)
        end
        calc_im = s -> begin
            fp = evaluate(c+s*im)
            fm = evaluate(c-s*im)
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

"""
Numerical Jacobian of `smn` outputs with respect to `c`.

By default, differentiate the normalized Legendre coefficient eigenvector and
evaluate its expansion at the requested coordinates. Eigenvalue and angular
sensitivities reuse the same internal coefficient calculation. For complex
parameters this uses the analytic bilinear normalization and the local branch
selected by `smn`. Imaginary-direction partials follow analyticity on that branch.
Supplying `h` selects centered finite differences with local phase transport.
At first-kind angular endpoints, the mixed coordinate/parameter derivative uses its
one-sided limit, which can be infinite.

`kind=2` differentiates Qs, including its DLMF normalization and differential
equation, on the same local branch as `smn(...; kind=2)`. It rejects
`normalize=true`. At the singular coordinates `eta=±1`, both parameter
derivatives return `NaN`; nearby interior coordinates remain supported.
Real first-kind derivatives with `abs(c) > n + 1` use extra working precision to
protect exponentially small interior values from cancellation.

For real `c`, returns:
- `dvalue_dc`
- `dderivative_dc`

For complex `c = a + ib`, returns:
- `dvalue_dcreal`, `dvalue_dcimag`
- `dderivative_dcreal`, `dderivative_dcimag`
"""
function jacobian_smn(m::Integer, n::Integer, c::Union{Real,Complex}, eta::AbstractVector{<:Real};
                      spheroid::Symbol=:prolate, precision::Symbol=:double, normalize::Bool=false, h=nothing,
                      kind::Integer=1,
                      with_metadata::Bool=false, adaptive::Bool=true,
                      rtol::Real=1e-6, atol::Real=1e-10)

    _validate_precision(precision)
    _validate_jacobian_tolerances(rtol, atol)
    _validate_wave_arguments(m,n,c,eta,spheroid,precision,:angular;kind)
    kind == 2 && normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1; Qs uses the DLMF second-kind normalization"))
    if h === nothing
        return kind == 2 ? _qs_angular_jacobian(m,n,c,eta,spheroid,precision,normalize,with_metadata) :
                           _coefficient_angular_jacobian(m,n,c,eta,spheroid,precision,normalize,with_metadata)
    end
    step = _resolve_jacobian_step(c, h, precision)
    if c isa Real
        cache = Dict{Any,Any}()
        evaluate(parameter) = get!(cache,parameter) do
            smn(m,n,parameter,eta;spheroid,precision,normalize,kind)
        end
        calc_value = s -> begin
            sp = evaluate(c+s)
            sm = evaluate(c-s)
            (sp.value .- sm.value) ./ (2 * s)
        end
        calc_derivative = s -> begin
            sp = evaluate(c+s)
            sm = evaluate(c-s)
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
        evaluate = _angular_jacobian_evaluator(m, n, c, eta; spheroid, precision, normalize,kind)
        calc_value_re = s -> begin
            sp = evaluate(c + s)
            sm = evaluate(c - s)
            (sp.value .- sm.value) ./ (2 * s)
        end
        calc_value_im = s -> begin
            sp = evaluate(c + s * im)
            sm = evaluate(c - s * im)
            (sp.value .- sm.value) ./ (2 * s)
        end
        calc_derivative_re = s -> begin
            sp = evaluate(c + s)
            sm = evaluate(c - s)
            (sp.derivative .- sm.derivative) ./ (2 * s)
        end
        calc_derivative_im = s -> begin
            sp = evaluate(c + s * im)
            sm = evaluate(c - s * im)
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

"""
Numerical Jacobian of `rmn` outputs with respect to `c`.

By default, differentiate the normalized spherical-Bessel expansions and their
coefficient eigenproblem. Near boundaries where the second-kind expansion
converges slowly, propagate the differentiated radial equation from x=2.
Both geometries, all four kinds, and real/complex parameters are supported.
Supplying `h` selects centered finite differences with local mode transport.
The analytic path reports no finite-difference step. Its coefficient residuals
and series-tail checks are diagnostics, not rigorous error bounds.

At the prolate boundary x=1, real first-kind derivatives use regular one-sided
limits; other kinds return NaN. Complex prolate calls require x>1, as for `rmn`.
Exactly zero c remains undefined.

With explicit `h`, function values and coordinate derivatives reuse each native
evaluation. Metadata then reports step consistency, not a rigorous accuracy bound.

For real `c`, returns:
- `dvalue_dc`
- `dderivative_dc`

For complex `c = a + ib`, returns:
- `dvalue_dcreal`, `dvalue_dcimag`
- `dderivative_dcreal`, `dderivative_dcimag`
"""
function jacobian_rmn(m::Integer, n::Integer, c::Union{Real,Complex}, x::AbstractVector{<:Real};
                      spheroid::Symbol=:prolate, precision::Symbol=:double, kind::Integer=1, h=nothing,
                      with_metadata::Bool=false, adaptive::Bool=true,
                      rtol::Real=1e-6, atol::Real=1e-10)

    _validate_precision(precision)
    _validate_radial_parameter(c)
    _validate_jacobian_tolerances(rtol, atol)
    _validate_wave_arguments(m,n,c,x,spheroid,precision,:radial;kind)
    h === nothing && return _radial_analytic_jacobian(m,n,c,x,spheroid,precision,kind,with_metadata)
    step = _resolve_jacobian_step(c,h,precision)
    differentiate = _finite_difference_with_metadata
    if c isa Real
        T = precision === :quad ? BigFloat : Float64
        c = T(c)
        cache = Dict{Any,Any}()
        evaluate(parameter) = get!(cache,parameter) do
            rmn(m,n,parameter,x;spheroid,precision,kind)
        end
        calc_value = s -> begin
            rp = evaluate(c+s)
            rm = evaluate(c-s)
            (rp.value .- rm.value) ./ (2 * s)
        end
        calc_derivative = s -> begin
            rp = evaluate(c+s)
            rm = evaluate(c-s)
            (rp.derivative .- rm.derivative) ./ (2 * s)
        end
        dvalue_dc, metadata_value = differentiate(calc_value, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dc, metadata_derivative = differentiate(calc_derivative, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)

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
        _validate_wave_arguments(m,n,c,x,spheroid,precision,:radial;kind)
        evaluate = _complex_local_evaluator(m,n,c,spheroid,precision;points=x,kind)
        calc_value_re = s -> begin
            rp = evaluate(c+s)
            rm = evaluate(c-s)
            (rp.value .- rm.value) ./ (2 * s)
        end
        calc_value_im = s -> begin
            rp = evaluate(c+s*im)
            rm = evaluate(c-s*im)
            (rp.value .- rm.value) ./ (2 * s)
        end
        calc_derivative_re = s -> begin
            rp = evaluate(c+s)
            rm = evaluate(c-s)
            (rp.derivative .- rm.derivative) ./ (2 * s)
        end
        calc_derivative_im = s -> begin
            rp = evaluate(c+s*im)
            rm = evaluate(c-s*im)
            (rp.derivative .- rm.derivative) ./ (2 * s)
        end

        dvalue_dcreal, metadata_value_dcreal = differentiate(calc_value_re, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dvalue_dcimag, metadata_value_dcimag = differentiate(calc_value_im, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dcreal, metadata_derivative_dcreal = differentiate(calc_derivative_re, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)
        dderivative_dcimag, metadata_derivative_dcimag = differentiate(calc_derivative_im, step; precision=precision, adaptive=adaptive, rtol=rtol, atol=atol)

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

"""
Solve `eigenvalue(m, n, c) = lambda_target` for real `c`.

With `operator=:concentration`, solve for real prolate bandwidth with `m=0`.
`form=:value` interprets the target as Λ in [0,1); `:complement` as 1-Λ in
(0,1]; `:log` as log(Λ) in [-Inf,0). Exact unit concentration has no finite
bandwidth. A zero-concentration target requires a bracket containing zero.
Use complement or log targets to preserve extreme tails. The returned residual
is the computed value minus the target in the selected form (zero at an exact
zero-concentration match, including `form=:log`). Fourier inversion is unsupported.

Uses a bracketed hybrid strategy with guaranteed bisection fallback and optional
Jacobian-guided acceleration through `jacobian_eigen` when derivative quality is
acceptable.

Keyword arguments:
- `bracket`: `(c_lo, c_hi)` with `c_lo < c_hi` and opposite signs of residual
  `eigenvalue(m,n,c) - lambda_target` at the endpoints.
- `spheroid`: `:prolate` or `:oblate`.
- `precision`: `:double` or `:quad`.
- `atol`, `rtol`: positive finite tolerances. Separation defaults are 1e-10 and
  1e-8, used for residual and bracket-width checks. Concentration defaults are
  1e-12/1e-10 (double), 1e-30/1e-28 (quad); the coordinate tolerance is
  `atol+rtol*abs(c)`. Concentration additionally requires the residual in
  `log(Λ/(1-Λ))` to be at most `rtol`, avoiding false convergence in either tail.
- `maxiter`: positive maximum number of iterations.
- `use_jacobian`: enable derivative-based candidate steps when trusted.

Returns a named tuple with fields:
- `converged::Bool`
- `c::T`
- `residual::T`
- `iterations::Int`
- `bracket::Tuple{T,T}`
- `method::Symbol` (`:endpoint`, `:newton`, `:secant`, `:bisection`, `:maxiter`)

`T` is `Float64` for double precision and `BigFloat` for quad precision.
"""
function find_c_for_eigenvalue(m::Integer, n::Integer, lambda_target::Real;
                               bracket::Tuple{<:Real,<:Real},
                               spheroid::Symbol=:prolate,
                               precision::Symbol=:double,
                               atol::Union{Nothing,Real}=nothing,
                               rtol::Union{Nothing,Real}=nothing,
                               maxiter::Integer=160,
                               use_jacobian::Bool=true,
                               operator::Symbol=:separation, form::Symbol=:value)

    _validate_precision(precision)
    operator in (:separation,:concentration) ||
        throw(ArgumentError("inverse bandwidth supports operator=:separation or :concentration"))
    operator===:separation && form!==:value &&
        throw(ArgumentError("the separation constant supports only form=:value"))
    T = precision===:quad ? BigFloat : Float64
    atol = atol===nothing ? (operator===:separation ? T(1e-10) : precision===:quad ? T(10)^(-30) : T(1e-12)) : T(atol)
    rtol = rtol===nothing ? (operator===:separation ? T(1e-8) : precision===:quad ? T(10)^(-28) : T(1e-10)) : T(rtol)
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end
    if !(isfinite(atol) && atol > 0)
        error("atol must be positive, got $atol")
    end
    if !(isfinite(rtol) && rtol > 0)
        error("rtol must be positive, got $rtol")
    end
    if maxiter <= 0
        error("maxiter must be positive, got $maxiter")
    end
    operator===:concentration && return _find_integral_bandwidth(m,n,lambda_target,
        bracket,spheroid,precision,form,atol,rtol,maxiter,use_jacobian)

    a = T(bracket[1])
    b = T(bracket[2])
    if !(isfinite(a) && isfinite(b) && a < b)
        error("bracket must satisfy bracket[1] < bracket[2]")
    end

    f(c) = T(eigenvalue(m, n, c; spheroid=spheroid, precision=precision) - lambda_target)

    fa = f(a)
    fb = f(b)
    if !isfinite(fa) || !isfinite(fb)
        error("non-finite residual at bracket endpoints")
    end

    tol_residual = atol + rtol * max(one(T), abs(T(lambda_target)))
    width_tol() = atol + rtol * max(one(T), abs(a), abs(b))

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
            d = T(j.derivative)
            md = j.metadata
            if isfinite(d) && abs(d) > sqrt(eps(T)) && md.suggested_action == :accept
                newton = mid - fmid / d
                if a < newton < b && isfinite(newton)
                    candidate = newton
                    method = :newton
                end
            end
        end

        if method == :bisection
            denom = fb - fa
            if isfinite(denom) && abs(denom) > eps(T)
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

"""
    accuracy(m, n, c, arg; target=:radial, spheroid=:prolate,
             precision=:double, kind=1, normalize=false)

Return a vector of backend-estimated decimal digits for the requested function
values. These are solver diagnostics, not rigorous error bounds or statistical
confidence intervals, and do not certify the returned coordinate derivatives.

An entry of `-1` means **no estimate is available**. Zero means the backend
reports no reliable decimal digits. Positive entries are estimates, not guarantees.
Invalid inputs raise the same domain/argument errors as `smn` or `rmn`.

- Angular `c=0`: evaluate the associated Legendre recurrence, then return `-1`;
  there is no error estimator for that path.
- Angular `kind=2`: evaluate Qs, then return `-1`; its propagation and
  coefficient convergence checks are not calibrated decimal-digit estimates.
- Refined angular evaluations for real parameters with `abs(c) > n + 1` return `-1`;
  native estimates do not describe the coefficient-based result.
- Angular endpoint limits and quad endpoint reconstruction: return `-1`; the native estimate
  describes the normalization point rather than the requested coordinate.
- Radial `c=0`: throw `DomainError`, as for `rmn`.
- Radial `kind=2`: report the native second-kind estimate, which may use
  Wronskian consistency and cancellation diagnostics.
- Radial `kind=1`, `3`, or `4`: evaluate the requested function, then return `-1`.
  A second-kind/pair estimate does not establish the accuracy of these functions
  or their potentially cancelling Hankel combinations.
- Nonfinite values and unavailable or out-of-range native estimates return `-1`.

Quad estimates use the same precision-preserving inputs as the evaluation path.
For sensitive calculations, use independent references, identities, or precision
comparisons in addition to the backend diagnostic.

```julia
accuracy(0, 2, 0.0, [0.3]; target=:angular) # [-1]
accuracy(0, 1, 2.0, [1.5]; target=:radial, kind=2)
```

"""
function accuracy(m::Integer, n::Integer, c::Union{Real,Complex}, arg::AbstractVector{<:Real};
                  spheroid::Symbol=:prolate, precision::Symbol=:double, kind::Integer=1, target::Symbol=:radial, normalize::Bool=false)

    target in (:angular, :radial) || error("target must be :angular or :radial")
    _validate_wave_arguments(m,n,c,arg,spheroid,precision,target;kind)
    if target === :angular && (_use_small_parameter_expansion(c) || kind == 2 || _use_angular_expansion(n,c,spheroid) || (c isa Complex && iszero(real(c))))
        smn(m,n,c,arg;spheroid,precision,normalize,kind)
        return fill(-1, length(arg)) # Neither evaluation has a calibrated digit estimator.
    elseif target === :radial && (kind != 2 || _radial_needs_analytic(c,arg,kind))
        rmn(m,n,c,arg;spheroid,precision,kind)
        return fill(-1, length(arg)) # Native naccr does not certify these kinds.
    elseif target === :angular && any(x -> abs(x)==1,arg)
        smn(m,n,c,arg;spheroid,precision,normalize)
        estimates=fill(-1,length(arg))
        interior=findall(x -> abs(x)!=1,arg)
        if !isempty(interior)
            estimates[interior]=accuracy(m,n,c,arg[interior];spheroid,precision,target,normalize)
        end
        return estimates # Native endpoint diagnostics do not cover series limits.
    elseif target === :angular && c isa Complex
        # A diagnostic must not accept a continuation path that smn rejects.
        smn(m,n,c,arg;spheroid,precision,normalize)
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
        selected_n = _complex_mode_state(prefix,m,n,c,precision).n
        if target === :angular
            return _call_complex_smn_accuracy(prefix, m, selected_n, c, arg; precision=precision, normalize=normalize)
        else
            return _call_complex_rmn_accuracy(prefix, m, selected_n, c, arg; precision=precision, kind=kind)
        end
    end
end

include("degree_ranges.jl")

function __init__()
    try
        # Default path for end users: shipped artifacts.
        _configure_backends_from_artifacts!()
        # Developer fallback when no artifact is available.
        _configure_backends_from_local_build!()
        # Overrides for CI/power users.
        _configure_backends_from_env!()
    catch e
        @warn "Failed to configure backend libraries during module initialization: $e" maxlog=1
    end
end

end

