# Guard precision must not change Base's process-wide BigFloat defaults on
# Julia 1.10/1.11. This private number type allocates every MPFR destination
# with an explicit precision. Only owned operands participate in its methods.
# Public results are converted back to ordinary BigFloat or Float64 values.
struct _SWFloat <: AbstractFloat
    value::BigFloat
    # Adopt an MPFR result at its existing precision.
    _SWFloat(value::BigFloat, ::Val{:raw}) = new(value)
end

# The private number type identifies its task-local arithmetic settings.
function _swsettings()
    settings = get(task_local_storage(), _SWFloat, nothing)
    settings === nothing &&
        return (Base.precision(BigFloat), Base.Rounding.rounding_raw(BigFloat))
    return settings::Tuple{Int, Base.MPFR.MPFRRoundingMode}
end
_swprecision() = first(_swsettings())
_swrounding() = convert(RoundingMode, last(_swsettings()))
_swmpfrrounding() = last(_swsettings())

function _with_swprecision(f::F, bits::Integer) where {F}
    bits > 0 || throw(ArgumentError("precision must be positive"))
    task_local_storage(f, _SWFloat, (Int(bits), _swmpfrrounding()))
end

function _SWFloat(x::Real)
    _SWFloat(BigFloat(x, _swmpfrrounding(); precision = _swprecision()), Val(:raw))
end
_SWFloat(x::_SWFloat) = _SWFloat(x.value)
function _SWFloat(x::AbstractString)
    _SWFloat(parse(BigFloat, strip(x); precision = _swprecision(), rounding = _swmpfrrounding()), Val(:raw))
end
function _SWFloat(x::Rational{T}) where {T}
    # Base's rational conversion also changes global precision on Julia 1.10.
    # Convert both integers exactly, then round only the quotient.
    bits = _swprecision()
    a = BigFloat(numerator(x); precision = max(bits, ndigits(numerator(x); base = 2)))
    b = BigFloat(denominator(x); precision = max(bits, ndigits(denominator(x); base = 2)))
    _SWFloat(a, Val(:raw)) / _SWFloat(b, Val(:raw))
end

Base.convert(::Type{_SWFloat}, x::_SWFloat) = x
Base.convert(::Type{_SWFloat}, x::Real) = _SWFloat(x)
Base.promote_rule(::Type{_SWFloat}, ::Type{T}) where {T <: Real} = _SWFloat
Base.promote_rule(::Type{_SWFloat}, ::Type{BigFloat}) = _SWFloat
Base.promote_rule(::Type{BigFloat}, ::Type{_SWFloat}) = _SWFloat
function Base.BigFloat(x::_SWFloat, r::RoundingMode = rounding(BigFloat);
        precision::Integer = Base.precision(BigFloat))
    BigFloat(x.value, r; precision)
end
Base.Float64(x::_SWFloat) = Float64(x.value)
Base.Float32(x::_SWFloat) = Float32(x.value)
Base.BigInt(x::_SWFloat) = BigInt(x.value)
Base.float(x::_SWFloat) = x
Base.float(::Type{_SWFloat}) = _SWFloat
Base.zero(::Type{_SWFloat}) = _SWFloat(0)
Base.one(::Type{_SWFloat}) = _SWFloat(1)
Base.zero(::_SWFloat) = zero(_SWFloat)
Base.one(::_SWFloat) = one(_SWFloat)
Base.precision(::Type{_SWFloat}) = _swprecision()
Base.precision(x::_SWFloat) = precision(x.value)
Base.rounding(::Type{_SWFloat}) = _swrounding()
Base.eps(::Type{_SWFloat}) = ldexp(one(_SWFloat), 1-_swprecision())
Base.eps(x::_SWFloat) = _SWFloat(eps(x.value), Val(:raw))
Base.decompose(x::_SWFloat) = Base.decompose(x.value)
Base.nextfloat(x::_SWFloat, n::Integer = 1) = _SWFloat(nextfloat(x.value, n), Val(:raw))
Base.prevfloat(x::_SWFloat, n::Integer = 1) = _SWFloat(prevfloat(x.value, n), Val(:raw))
Base.hash(x::_SWFloat, h::UInt) = hash(x.value, h)
Base.show(io::IO, x::_SWFloat) = show(io, x.value)
Base.string(x::_SWFloat) = string(x.value)
for predicate in (:isfinite, :isinf, :isnan, :iszero, :isone, :signbit, :isinteger, :exponent)
    @eval Base.$predicate(x::_SWFloat) = Base.$predicate(x.value)
end
for comparison in (:(==), :(<), :(<=), :isless)
    @eval Base.$comparison(x::_SWFloat, y::_SWFloat) = Base.$comparison(x.value, y.value)
end

for (operation, symbol) in ((:+, :mpfr_add), (:-, :mpfr_sub), (:*, :mpfr_mul),
    (:/, :mpfr_div), (:^, :mpfr_pow), (:hypot, :mpfr_hypot),
    (:atan, :mpfr_atan2), (:copysign, :mpfr_copysign))
    @eval function Base.$operation(x::_SWFloat, y::_SWFloat)
        result = BigFloat(; precision = _swprecision())
        ccall(($(QuoteNode(symbol)), Base.MPFR.libmpfr), Cint,
            (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
            result, x.value, y.value, _swmpfrrounding())
        _SWFloat(result, Val(:raw))
    end
end
for (operation, symbol) in ((:-, :mpfr_neg), (:abs, :mpfr_abs), (:sqrt, :mpfr_sqrt),
    (:exp, :mpfr_exp), (:expm1, :mpfr_expm1), (:log, :mpfr_log),
    (:log1p, :mpfr_log1p), (:log2, :mpfr_log2), (:log10, :mpfr_log10),
    (:sin, :mpfr_sin), (:cos, :mpfr_cos), (:sinh, :mpfr_sinh),
    (:cosh, :mpfr_cosh), (:sinpi, :mpfr_sinpi), (:cospi, :mpfr_cospi))
    @eval function Base.$operation(x::_SWFloat)
        result = BigFloat(; precision = _swprecision())
        ccall(($(QuoteNode(symbol)), Base.MPFR.libmpfr), Cint,
            (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
            result, x.value, _swmpfrrounding())
        _SWFloat(result, Val(:raw))
    end
end
Base.:+(x::_SWFloat) = x
Base.abs2(x::_SWFloat) = x*x
Base.inv(x::_SWFloat) = one(x)/x
Base.sincos(x::_SWFloat) = (sin(x), cos(x))
Base.sincospi(x::_SWFloat) = (sinpi(x), cospi(x))
function Base.:^(x::_SWFloat, y::Integer)
    result = BigFloat(; precision = _swprecision())
    exponent = BigInt(y)
    ccall((:mpfr_pow_z, Base.MPFR.libmpfr), Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Ref{BigInt}, Base.MPFR.MPFRRoundingMode),
        result, x.value, exponent, _swmpfrrounding())
    _SWFloat(result, Val(:raw))
end
Base.:^(x::_SWFloat, y::Rational) = x^_SWFloat(y)
Base.:^(x::BigFloat, y::_SWFloat) = _SWFloat(x)^y
function Base.ldexp(x::_SWFloat, exponent::Integer)
    result = BigFloat(; precision = _swprecision())
    ccall((:mpfr_mul_2si, Base.MPFR.libmpfr), Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Clong, Base.MPFR.MPFRRoundingMode),
        result, x.value, exponent, _swmpfrrounding())
    _SWFloat(result, Val(:raw))
end
for operation in (:ceil, :floor, :trunc, :round)
    @eval Base.$operation(::Type{T}, x::_SWFloat) where {T <: Integer} = Base.$operation(T, x.value)
    @eval Base.$operation(::Type{Bool}, x::_SWFloat) = Bool(Base.$operation(BigInt, x.value))
end
function _swround(x::_SWFloat, mode::RoundingMode)
    result = BigFloat(; precision = _swprecision())
    ccall((:mpfr_rint, Base.MPFR.libmpfr), Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        result, x.value, convert(Base.MPFR.MPFRRoundingMode, mode))
    _SWFloat(result, Val(:raw))
end
Base.round(x::_SWFloat, mode::RoundingMode = RoundNearest) = _swround(x, mode)
Base.round(x::_SWFloat, mode::RoundingMode{:FromZero}) = _swround(x, mode)
function Base.round(x::_SWFloat, ::RoundingMode{:NearestTiesAway})
    result = BigFloat(; precision = _swprecision())
    ccall((:mpfr_round, Base.MPFR.libmpfr), Cint,
        (Ref{BigFloat}, Ref{BigFloat}), result, x.value)
    _SWFloat(result, Val(:raw))
end
function Base.round(x::_SWFloat, ::RoundingMode{:NearestTiesUp})
    result = round(x, RoundNearest)
    # Removing the integer part is exact at the operand's own precision.
    fraction = BigFloat(; precision = precision(x))
    ccall((:mpfr_frac, Base.MPFR.libmpfr), Cint,
        (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode),
        fraction, x.value, convert(Base.MPFR.MPFRRoundingMode, RoundNearest))
    return (fraction == 0.5 || fraction == -0.5) && result < x ? result+one(x) : result
end
Base.ceil(x::_SWFloat) = round(x, RoundUp)
Base.floor(x::_SWFloat) = round(x, RoundDown)
Base.trunc(x::_SWFloat) = round(x, RoundToZero)

# Preserve guard precision in private diagnostics without exposing the wrapper.
_unwrap_swfloat(x::_SWFloat) = x.value
_unwrap_swfloat(x) = x

# Julia 1.10/1.11's BigFloat(::Rational) changes global arithmetic settings.
# Use explicit destinations when accepting exact user inputs as well.
_input_float(::Type{T}, x) where {T} = T(x)
function _input_float(::Type{BigFloat}, x::Rational)
    _with_swprecision(Base.precision(BigFloat)) do
        _SWFloat(x).value
    end
end
function _input_float(::Type{Complex{T}}, x::Number) where {T}
    complex(_input_float(T, real(x)), _input_float(T, imag(x)))
end
_input_bigfloat(x) = _input_float(BigFloat, x)

function _stack_wave_results(results)
    fields = keys(first(results))
    values = map(fields) do field
        entries = [getproperty(r, field) for r in results]
        if first(entries) isa NamedTuple
            (; mantissa = hcat((v.mantissa for v in entries)...),
                exponent = hcat((v.exponent for v in entries)...))
        else
            hcat(entries...)
        end
    end
    return NamedTuple{fields}(values)
end

function _scaled_native_values(m, n, c, points, spheroid, precision, target, option)
    if target===:radial && _radial_needs_analytic(c, points, option)
        return _radial_analytic_values(m, n, c, points, spheroid, precision, option)
    end
    if target === :radial && spheroid === :prolate && c isa Complex && any(isone, points)
        return _radial_analytic_values(m, n, c, points, spheroid, precision, option)
    end
    if target === :angular && iszero(c)
        return _angular_phase!(
            _spherical_smn_real(
                m, n, _input_bigfloat.(points); normalize = option!=0, precision = :quad),
            m)
    end
    if target === :angular && _use_small_parameter_expansion(c)
        return _angular_phase!(
            _small_parameter_smn(
                m, n, c, points, spheroid, precision, option!=0), m)
    end
    requested_n = n
    state = c isa Complex ?
            _complex_mode_state(
        spheroid === :prolate ? :cprolate : :coblate, m, n, c, precision) : nothing
    state !== nothing && target === :angular && _require_angular_anchor(state)
    state !== nothing && (n = state.n)
    boundary = target === :radial && spheroid === :prolate && option >= 2 ?
               findall(isone, points) : Int[]
    if !isempty(boundary)
        value = fill(complex(BigFloat(NaN), big"0"), length(points))
        derivative = copy(value)
        interior = findall(x -> !isone(x), points)
        if !isempty(interior)
            r = _scaled_native_values(
                m, requested_n, c, points[interior], spheroid, precision, target, option)
            value[interior], derivative[interior] = r.value, r.derivative
        end
        return (; value, derivative)
    end
    lib = _require_backend_library(precision)
    fn = Libdl.dlsym_e(_require_backend_handle(lib), :spheroidal_scaled_text)
    fn == C_NULL &&
        error("Backend lacks scaled output support; rebuild the native libraries with Pkg.build(\"SpheroidalWaves\").")
    endpoint = target === :angular ?
               _angular_endpoint_plan(m, n, c, points, spheroid; precision) : nothing
    native_points = endpoint === nothing ? _input_bigfloat.(points) : endpoint.native_points
    target === :radial && spheroid === :prolate && (native_points = native_points .- 1)
    ctext = _format_fortran_input([real(c), imag(c)])
    xtext = _format_fortran_input(native_points)
    count = 8length(points)
    output = fill(UInt8(' '), _QUAD_TEXT_WIDTH*count)
    exponents = zeros(Cint, count)
    status = Ref{Cint}(0)
    ccall(fn,
        Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
            Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        spheroid === :oblate, c isa Complex, target === :angular ? 1 : 2, option, m, n, length(points),
        _QUAD_TEXT_WIDTH, ctext, xtext, output, exponents, status)
    _check_scalar_status(status[])
    data = reshape(_parse_fortran_output(output, exponents, count), 8, :)
    value = complex.(data[1, :], data[2, :])
    derivative = complex.(data[3, :], data[4, :])
    if target === :angular
        _angular_endpoint_reconstruct!(value, derivative, endpoint, m, n)
        if c isa Complex
            factor = _angular_mode_factor(m, requested_n, state, option!=0, BigFloat)
            value .*= factor
            derivative .*= factor
        end
        return _angular_phase!((; value, derivative), m)
    end
    if option != 1
        second, dsecond = complex.(data[5, :], data[6, :]), complex.(data[7, :], data[8, :])
        if option == 2
            value, derivative = second, dsecond
        else
            phase = option == 3 ? im : -im
            value .+= phase .* second
            derivative .+= phase .* dsecond
        end
    end
    result = (; value, derivative)
    return state === nothing ? result :
           _scale_mode_result!(result, _radial_mode_factor(requested_n, state))
end

function _decimal_scaled(values, T)
    exponents = [iszero(v) || !isfinite(v) ? 0 : floor(Int, log10(abs(v))) for v in values]
    mantissas = T[v/BigFloat(10)^e for (v, e) in zip(values, exponents)]
    for i in eachindex(mantissas)
        if isfinite(mantissas[i]) && abs(mantissas[i]) >= 10
            mantissas[i] /= 10
            exponents[i] += 1
        end
    end
    return (; mantissa = mantissas, exponent = exponents)
end

function _extended_wave(m, n, c, points, spheroid, precision, target,
        option, scaled, logderivative, second_derivative)
    result = _scaled_native_values(m, n, c, points, spheroid, precision, target, option)
    _fix_regular_endpoints!(result, m, n, c, points, spheroid, precision, target, option)
    second_derivative && (result = _wave_second_derivative(
        result, m, n, c, points, spheroid, precision, target; option))
    R = precision === :quad ? BigFloat : Float64
    T = target === :angular && c isa Real ? R : Complex{R}
    # Real angular calls discard only the exactly-zero imaginary storage channel.
    convert_values(v) = T <: Real ? real.(v) : v
    output = map(
        v -> scaled ? _decimal_scaled(convert_values(v), T) :
             T.(convert_values(v)), result)
    if logderivative
        ratio = [_wave_logderivative(
                     v, d, x, m, spheroid, target, target===:angular ? 1 : option)
                 for (v, d, x) in zip(result.value, result.derivative, points)]
        output = (; output..., logderivative = T.(convert_values(ratio)))
    end
    return output
end
