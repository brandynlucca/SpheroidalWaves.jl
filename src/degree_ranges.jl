function smn(m::Integer, n::AbstractUnitRange{<:Integer},
        c::Union{Real, Complex}, z::Complex; kwargs...)
    smn(m, n, c, [z]; kwargs...)
end
function rmn(m::Integer, n::AbstractUnitRange{<:Integer},
        c::Union{Real, Complex}, z::Complex; kwargs...)
    rmn(m, n, c, [z]; kwargs...)
end

function smn(m::Integer, n::AbstractUnitRange{<:Integer},
        c::Union{Real, Complex}, z::AbstractVector{<:Number}; kwargs...)
    _validate_degree_range(
        m, n, z, get(kwargs, :spheroid, :prolate), get(kwargs, :precision, :double))
    _stack_wave_results([smn(m, degree, c, z; kwargs...) for degree in n])
end

function rmn(m::Integer, n::AbstractUnitRange{<:Integer},
        c::Union{Real, Complex}, z::AbstractVector{<:Number}; kwargs...)
    _validate_degree_range(
        m, n, z, get(kwargs, :spheroid, :prolate), get(kwargs, :precision, :double))
    _stack_wave_results([rmn(m, degree, c, z; kwargs...) for degree in n])
end

function _validate_degree_range(m, n, points, spheroid, precision)
    _validate_precision(precision)
    spheroid in (:prolate, :oblate) || throw(ArgumentError("invalid spheroid: $spheroid"))
    !isempty(n) && 0 <= m <= first(n) ||
        throw(ArgumentError("require a nonempty degree range with 0 ≤ m ≤ first(n)"))
    isempty(points) && throw(ArgumentError("evaluation points must not be empty"))
    all(isfinite, points) || throw(ArgumentError("evaluation points must be finite"))
    return nothing
end

function rmn(m::Integer, n::AbstractUnitRange{<:Integer},
        c::Union{Real, Complex}, x::Real; kwargs...)
    rmn(m, n, c, [x]; kwargs...)
end

function _degree_range_pointer(symbol, c, spheroid, precision)
    if spheroid === :prolate && precision === :quad && c isa Real && c > 0
        lib = _require_backend_library(precision)
        return Libdl.dlsym_e(_require_backend_handle(lib), symbol)
    end
    return C_NULL
end

# One native expansion returns every requested degree. Keep mantissas scaled
# until combinations and optional ratios have been formed in a wide type.
function _shared_real_degree_range(m, n, c, points, spheroid, precision, target, option,
        scaled, logderivative, second_derivative)
    c isa Real && c>0 || return nothing
    target === :angular && _use_small_parameter_expansion(c) && return nothing
    # The earlier quad prolate ABI already shares degrees and transfers fewer
    # channels for ordinary output. Retain it when no extended fields are needed.
    if spheroid===:prolate && precision===:quad && !scaled && !logderivative &&
       !second_derivative &&
       (target===:radial || !_use_angular_expansion(first(n), c, spheroid))
        return nothing
    end
    if target===:angular
        any(x->_near_angular_endpoint(x) || abs(x)==1, points) && return nothing
    else
        _radial_needs_analytic(c, points, option) && return nothing
        spheroid===:prolate && any(x->x<=1, points) && return nothing
    end
    lib = _require_backend_library(precision)
    symbol = precision===:quad ? :spheroidal_degrees_scaled_text :
             :spheroidal_degrees_scaled_double
    pointer = Libdl.dlsym_e(_require_backend_handle(lib), symbol)
    pointer==C_NULL && return nothing
    native_first = first(n)
    if target===:angular
        while native_first<=last(n) && _use_angular_expansion(native_first, c, spheroid)
            native_first+=1
        end
        native_first>last(n) && return nothing
    end
    native_degrees = native_first:last(n)
    native_points = _input_bigfloat.(points)
    target===:radial && spheroid===:prolate && (native_points .-= 1)
    ctext = _format_fortran_input([c, zero(c)])
    xtext = _format_fortran_input(native_points)
    count = 8length(points)*length(native_degrees)
    exponents = zeros(Cint, count)
    status = Ref{Cint}(0)
    if precision===:quad
        output = fill(UInt8(' '), _QUAD_TEXT_WIDTH*count)
        ccall(pointer,
            Cvoid,
            (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
                Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
            spheroid===:oblate, target===:angular ? 1 : 2, option, m, native_first, last(n), length(points),
            _QUAD_TEXT_WIDTH, ctext, xtext, output, exponents, status)
        _check_scalar_status(status[])
        data = _parse_fortran_output(output, exponents, count)
    else
        output = zeros(Float64, count)
        ccall(pointer,
            Cvoid,
            (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
                Ptr{UInt8}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cint}, Ref{Cint}),
            spheroid===:oblate, target===:angular ? 1 : 2, option, m, native_first, last(n), length(points),
            _QUAD_TEXT_WIDTH, ctext, xtext, output, exponents, status)
        _check_scalar_status(status[])
        data = [BigFloat(v)*BigFloat(10)^e for (v, e) in zip(output, exponents)]
    end
    channels = reshape(data, 8, length(points), :)
    R = precision===:quad ? BigFloat : Float64
    T = target===:angular ? R : Complex{R}
    results = map(enumerate(native_degrees)) do (column, degree)
        value, derivative = channels[1, :, column], channels[3, :, column]
        if target===:radial
            if option==2
                value, derivative = channels[5, :, column], channels[7, :, column]
            elseif option>=3
                sign = option==3 ? 1 : -1
                value = complex.(value, sign .* channels[5, :, column])
                derivative = complex.(derivative, sign .* channels[7, :, column])
            end
            value, derivative = complex.(value), complex.(derivative)
        end
        result = (; value, derivative)
        target===:angular && _angular_phase!(result, m)
        _fix_regular_endpoints!(
            result, m, degree, c, points, spheroid, precision, target, option)
        second_derivative && (result = _wave_second_derivative(
            result, m, degree, c, points, spheroid, precision, target; option))
        converted = map(v->scaled ? _decimal_scaled(v, T) : T.(v), result)
        if logderivative
            ratio = T[_wave_logderivative(
                          v, d, x, m, spheroid, target, target===:angular ? 1 : option)
                      for (v, d, x) in zip(result.value, result.derivative, points)]
            converted = (; converted..., logderivative = ratio)
        end
        return converted
    end
    if native_first>first(n)
        refined = [smn(m, degree, c, points; spheroid, precision, normalize = option!=0,
                       scaled, logderivative, second_derivative)
                   for degree in first(n):(native_first - 1)]
        return _stack_wave_results(vcat(refined, results))
    end
    return _stack_wave_results(results)
end

"""
    smn(m, n::AbstractUnitRange{<:Integer}, c, eta; kwargs...)

Evaluate angular wave functions with `n` given as a unit range of degrees. The `value` and
`derivative` matrices have rows corresponding to `eta` and columns corresponding
to `n`. Normalization, precision, and spheroid keywords match the
single-degree method; a scalar `eta` gives one row.

Use `precision=:quad` to request quad precision with the same numeric arguments.
"""
function smn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real, Complex},
        eta::AbstractVector{<:Real}; spheroid::Symbol = :prolate,
        precision::Symbol = :double, normalize::Bool = false, kind::Integer = 1, scaled::Bool = false,
        logderivative::Bool = false, second_derivative::Bool = false, derivatives::Integer = 1)
    _validate_degree_range(m, n, eta, spheroid, precision)
    isfinite(c) || throw(ArgumentError("c must be finite"))
    all(x -> abs(x) <= 1, eta) || throw(ArgumentError("eta must lie in [-1, 1]"))
    kind in (1, 2) || throw(ArgumentError("angular kind must be 1 or 2"))
    derivatives in 1:4 || throw(ArgumentError("derivatives must be in 1:4"))
    derivatives > 2 && return _stack_wave_results([smn(m, degree, c, eta;
                                    spheroid, precision, normalize, kind, scaled,
                                    logderivative, derivatives) for degree in n])
    second_derivative |= derivatives == 2
    if kind==1
        shared = _shared_real_degree_range(
            m, n, c, eta, spheroid, precision, :angular, Int(normalize),
            scaled, logderivative, second_derivative)
        shared!==nothing && return shared
    end
    if kind == 2 || scaled || logderivative || second_derivative ||
       _use_angular_expansion(first(n), c, spheroid) || _use_small_parameter_expansion(c)
        return _stack_wave_results([smn(m, degree, c, eta; spheroid, precision, normalize,
                                        kind, scaled, logderivative, second_derivative)
                                    for degree in n])
    end
    pointer = _degree_range_pointer(:psms_smn_degrees_quad_text, c, spheroid, precision)
    any(x -> _near_angular_endpoint(x) || abs(x)==1, eta) && (pointer = C_NULL)
    if pointer == C_NULL
        results = [smn(m, degree, c, eta; spheroid, precision, normalize) for degree in n]
        return (; value = hcat((r.value for r in results)...),
            derivative = hcat((r.derivative for r in results)...))
    end
    count = length(eta)*length(n)
    c_text = _format_fortran_input(c)
    eta_text = _format_fortran_input(eta)
    value_text = fill(UInt8(' '), _QUAD_TEXT_WIDTH*count)
    derivative_text = similar(value_text)
    value_exp = zeros(Cint, count)
    derivative_exp = similar(value_exp)
    status = Ref{Cint}(0)
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8},
            Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        Cint(m), Cint(first(n)), Cint(last(n)), Cint(length(eta)),
        _bool_to_cint(normalize), c_text, Cint(_QUAD_TEXT_WIDTH), eta_text,
        value_text, value_exp, derivative_text, derivative_exp, status)
    _check_scalar_status(status[])
    result = (;
        value = reshape(_parse_fortran_output(value_text, value_exp, count), length(eta), :),
        derivative = reshape(_parse_fortran_output(derivative_text, derivative_exp, count), length(eta), :))
    return _angular_phase!(result, m)
end

function smn(m::Integer, n::AbstractUnitRange{<:Integer},
        c::Union{Real, Complex}, eta::Real; kwargs...)
    smn(m, n, c, [eta]; kwargs...)
end

"""
    rmn(m, n::AbstractUnitRange{<:Integer}, c, x; kwargs...)

Evaluate radial wave functions with `n` given as a unit range of degrees. The `value` and
`derivative` matrices have rows corresponding to `x` and columns corresponding
to `n`. All four radial kinds and the spheroid and precision keywords are
supported. A scalar `x` gives one row. Use `precision=:quad` to request quad
precision with the same numeric arguments.
"""
function rmn(m::Integer, n::AbstractUnitRange{<:Integer}, c::Union{Real, Complex},
        x::AbstractVector{<:Real}; spheroid::Symbol = :prolate,
        precision::Symbol = :double, kind::Integer = 1, scaled::Bool = false,
        logderivative::Bool = false, second_derivative::Bool = false, derivatives::Integer = 1,
        normalization::Symbol = :standard)
    _validate_degree_range(m, n, x, spheroid, precision)
    isfinite(c) || throw(ArgumentError("c must be finite"))
    kind in 1:4 || throw(ArgumentError("kind must be in 1:4"))
    _validate_wave_arguments(
        m, first(n), c, x, spheroid, precision, :radial; kind, normalization)
    derivatives in 1:4 || throw(ArgumentError("derivatives must be in 1:4"))
    (derivatives > 2 || normalization === :static) &&
        return _stack_wave_results([rmn(m, degree, c, x;
                                        spheroid, precision, kind, scaled, logderivative,
                                        second_derivative, derivatives, normalization)
                                    for degree in n])
    second_derivative |= derivatives == 2
    shared = _shared_real_degree_range(m, n, c, x, spheroid, precision, :radial, kind,
        scaled, logderivative, second_derivative)
    shared!==nothing && return shared
    if scaled || logderivative || second_derivative
        return _stack_wave_results([rmn(m, degree, c, x; spheroid, precision, kind,
                                        scaled, logderivative, second_derivative)
                                    for degree in n])
    end
    pointer = all(t -> t > 1, x) ?
              _degree_range_pointer(:psms_rmn_degrees_quad_text, c, spheroid, precision) :
              C_NULL
    if pointer == C_NULL
        results = [rmn(m, degree, c, x; spheroid, precision, kind) for degree in n]
        return (; value = hcat((r.value for r in results)...),
            derivative = hcat((r.derivative for r in results)...))
    end
    count = 4*length(x)*length(n)
    c_text = _format_fortran_input(c)
    x_text = _format_fortran_input(x)
    output = fill(UInt8(' '), _QUAD_TEXT_WIDTH*count)
    exponents = zeros(Cint, count)
    status = Ref{Cint}(0)
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8}, Cint,
            Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        Cint(m), Cint(first(n)), Cint(last(n)), Cint(kind), c_text,
        Cint(length(x)), x_text, Cint(_QUAD_TEXT_WIDTH), output, exponents, status)
    _check_scalar_status(status[])
    data = reshape(_parse_fortran_output(output, exponents, count), 4, length(x), :)
    if kind == 1 || kind == 2
        channel = kind == 1 ? 1 : 3
        return (; value = complex.(data[channel, :, :]),
            derivative = complex.(data[channel + 1, :, :]))
    end
    sign = kind == 3 ? 1 : -1
    return (; value = complex.(data[1, :, :], sign .* data[3, :, :]),
        derivative = complex.(data[2, :, :], sign .* data[4, :, :]))
end

# Share the real-prolate expansion when the native ABI and domains permit it.
function _native_prolate_block(m, n, c, eta, x, precision, normalize, kind,
        derivatives, terms)
    # Endpoints, complex inputs, and small parameters retain their existing evaluators.
    c isa Real && c > 0 && !_use_small_parameter_expansion(c) || return nothing
    isempty(eta) && isempty(x) && terms == 0 && return nothing
    all(t -> t isa Real && abs(t) < 1 && !_near_angular_endpoint(t), eta) || return nothing
    all(t -> t isa Real && t > 1, x) || return nothing
    !isempty(x) && _radial_needs_analytic(c, x, kind) && return nothing
    library = _require_backend_library(precision)
    pointer = Libdl.dlsym_e(_require_backend_handle(library), :spheroidal_prolate_block_text)
    # Older artifacts remain usable without the additive native entry point.
    pointer == C_NULL && return nothing
    ctext = _format_fortran_input(c)
    etatext = _format_fortran_input(Real[t for t in eta])
    # Form x-1 before native conversion to preserve distance from the prolate boundary.
    xtext = _format_fortran_input(BigFloat[_input_bigfloat(t)-1 for t in x])
    # Native wave arrays place channels first, then coordinates, then degrees.
    angular_count = (1+derivatives)*length(eta)*length(n)
    radial_count = 2(1+derivatives)*length(x)*length(n)
    angular = fill(UInt8(' '), _QUAD_TEXT_WIDTH*angular_count)
    radial = fill(UInt8(' '), _QUAD_TEXT_WIDTH*radial_count)
    # Coefficients use term rows and degree columns, with one count per degree.
    coefficients = fill(UInt8(' '), _QUAD_TEXT_WIDTH*terms*length(n))
    eigenvalues = fill(UInt8(' '), _QUAD_TEXT_WIDTH*length(n))
    aexp, rexp = zeros(Cint, angular_count), zeros(Cint, radial_count)
    counts, status = zeros(Cint, length(n)), Ref{Cint}(0)
    # Decimal mantissas and separate exponents preserve quad precision and range.
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
            Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint},
            Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ref{Cint}),
        m, first(n), last(n), length(eta), length(x), normalize, kind, derivatives, terms,
        _QUAD_TEXT_WIDTH, ctext, etatext, xtext, angular, aexp, radial, rexp,
        coefficients, counts, eigenvalues, status)
    # Reject invalid native output before parsing uninitialized text buffers.
    _check_scalar_status(status[])
    # Parse into wide values before the public API applies its requested result type.
    return (;
        angular = reshape(_parse_fortran_output(angular, aexp, angular_count),
            1+derivatives, length(eta), length(n)),
        radial = reshape(_parse_fortran_output(radial, rexp, radial_count),
            2(1+derivatives), length(x), length(n)),
        coefficients = reshape(_parse_fortran_output(coefficients, terms*length(n)), terms, length(n)),
        counts, eigenvalues = _parse_fortran_output(eigenvalues, length(n)))
end

"""
    spheroidal(m, n, c; eta=nothing, x=nothing, coefficients=false,
               spheroid=:prolate, precision=:double, normalize=false, kind=1,
               derivatives=true, parameter_derivatives=false, basis=:ferrers,
               rtol=nothing, max_terms=512)

Evaluate requested angular functions, radial functions, and expansion
coefficients together. Supply `eta` for first-kind angular functions, `x` for
radial functions of the selected `kind`, and `coefficients=true` for `dmn`
results. Coordinates can be scalars or vectors. The result contains only the
requested `angular`, `radial`, and `coefficients` fields. Both
`spheroid=:prolate` and `spheroid=:oblate` are supported.

For integer `n`, wave results have coordinate vectors and the coefficient
result is a named tuple. For a unit range of degrees, wave results have
coordinate rows and degree columns, and coefficients are a vector of named
tuples. All conventions and domains match `smn`, `rmn`, and `dmn`.

`derivatives=false` omits coordinate derivative outputs. Coefficient parameter
sensitivities are separately requested with `parameter_derivatives=true`.
`normalize` controls angular normalization, while `basis` selects the
coefficient basis. `rtol` and `max_terms` control the coefficient expansion.

Compatible native libraries share one real-prolate expansion between angular,
radial, and coefficient results at fixed `(m,c)`. Additional radial coordinates
require separate radial work. Guarded angular evaluation and coefficient
accuracy checks remain active. Oblate evaluations use the existing oblate
evaluators. Other supported parameters, coordinates, and older libraries also
use the existing individual evaluators.
"""
function spheroidal(m::Integer, n::Union{Integer, AbstractUnitRange{<:Integer}},
        c::Union{Real, Complex}; eta = nothing, x = nothing, coefficients::Bool = false,
        spheroid::Symbol = :prolate, precision::Symbol = :double, normalize::Bool = false,
        kind::Integer = 1, derivatives::Bool = true, parameter_derivatives::Bool = false,
        basis::Symbol = :ferrers, rtol = nothing, max_terms::Integer = 512)
    eta === nothing && x === nothing && !coefficients &&
        throw(ArgumentError("request eta, x, or coefficients"))
    # Validate geometry before applying its angular and radial coordinate domains.
    degrees = n isa Integer ? (n:n) : n
    _validate_degree_range(m, degrees, [0], spheroid, precision)
    _validate_wave_arguments(m, first(degrees), c, [0], spheroid, precision, :angular)
    kind in 1:4 || throw(ArgumentError("kind must be in 1:4"))
    # Empty internal arrays suppress omitted outputs without changing scalar result shapes.
    angular_points = eta === nothing ? Real[] : eta isa Number ? [eta] : eta
    radial_points = x === nothing ? Real[] : x isa Number ? [x] : x
    if eta !== nothing
        _validate_degree_range(m, degrees, angular_points, spheroid, precision)
        if eltype(angular_points) <: Real
            _validate_wave_arguments(m, first(degrees), c, angular_points,
                spheroid, precision, :angular)
        end
    end
    if x !== nothing
        _validate_degree_range(m, degrees, radial_points, spheroid, precision)
        if eltype(radial_points) <: Real
            _validate_wave_arguments(m, first(degrees), c, radial_points,
                spheroid, precision, :radial; kind)
        end
    end
    coefficients && _validate_coefficient_options(m, last(degrees), c, spheroid,
        precision, rtol, max_terms, basis)
    # Native coefficients contain values only, so sensitivities use the guarded solver.
    terms = coefficients && !parameter_derivatives ?
            min(max_terms, max(32, div(last(degrees)-m, 2)+ceil(Int, abs(c))+32)) : 0
    # The public API covers both geometries, while this shared native ABI is prolate-only.
    block = spheroid === :prolate ?
            _native_prolate_block(m, degrees, c, angular_points, radial_points,
        precision, normalize, kind, derivatives, terms) : nothing
    result = NamedTuple()
    T = precision === :quad ? BigFloat : Float64
    if eta !== nothing
        if block === nothing
            # Forward geometry unchanged to preserve oblate and analytic continuation behavior.
            angular = smn(m, n, c, angular_points; spheroid, precision, normalize)
            !derivatives && (angular = (; value = angular.value))
        else
            angular = (; value = T.(block.angular[1, :, :]))
            derivatives &&
                (angular = (; angular..., derivative = T.(block.angular[2, :, :])))
            # Native angular values omit the Condon-Shortley phase used by the public API.
            if isodd(m)
                for field in angular
                    field .*= -1
                end
            end
            # Cancellation-prone modes retain guarded pointwise angular evaluation.
            for (i, degree) in enumerate(degrees)
                if _use_angular_expansion(degree, c, spheroid) ||
                   any(v -> !isfinite(v), angular.value[:, i])
                    guarded = smn(
                        m, degree, c, angular_points; spheroid, precision, normalize)
                    angular.value[:, i] = guarded.value
                    derivatives && (angular.derivative[:, i] = guarded.derivative)
                end
            end
            # Scalar degrees return coordinate vectors, while degree ranges retain columns.
            n isa Integer && (angular = map(v -> vec(v), angular))
        end
        result = (; result..., angular)
    end
    if x !== nothing
        if block === nothing
            # The oblate fallback includes x=0 and uses the established radial conventions.
            radial = rmn(m, n, c, radial_points; spheroid, precision, kind)
            !derivatives && (radial = (; value = radial.value))
        else
            # Native channels hold R1 and R2, with adjacent derivatives only when requested.
            second = 2+derivatives
            # Construct outgoing and incoming waves as R1 plus or minus i*R2.
            combine(a, b) = kind == 1 ? complex.(a) :
                            kind == 2 ? complex.(b) :
                            complex.(a, (kind == 3 ? 1 : -1) .* b)
            value = combine(block.radial[1, :, :], block.radial[second, :, :])
            radial = (; value = Complex{T}.(value))
            if derivatives
                slope = combine(block.radial[2, :, :], block.radial[second + 1, :, :])
                radial = (; radial..., derivative = Complex{T}.(slope))
            end
            n isa Integer && (radial = map(v -> vec(v), radial))
        end
        result = (; result..., radial)
    end
    if coefficients
        # Reuse the shared expansion only after coefficient residual and tail checks pass.
        native = parameter_derivatives || block === nothing ? nothing :
                 _native_prolate_coefficient_block(m, degrees, c, spheroid, precision,
            normalize, basis, rtol, max_terms; block)
        # A rejected degree falls back independently without discarding accepted neighbors.
        values = map(enumerate(degrees)) do (i, degree)
            native !== nothing && native[i] !== nothing && return native[i]
            _julia_dmn(m, degree, c; spheroid, precision, normalize, rtol, max_terms,
                derivatives = parameter_derivatives, basis)
        end
        result = (; result..., coefficients = n isa Integer ? only(values) : values)
    end
    return result
end
