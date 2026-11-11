function _precision_wave_change(previous, current; nonfinite_indices = ())
    keys(previous)==keys(current) || return _SWFloat(Inf)
    change = zero(_SWFloat)
    for name in keys(current)
        left, right = getproperty(previous, name), getproperty(current, name)
        length(left)==length(right) || return _SWFloat(Inf)
        for i in eachindex(right)
            a, b = left[i], right[i]
            if !isfinite(a) || !isfinite(b)
                matching = isequal(real(a), real(b)) && isequal(imag(a), imag(b))
                i in nonfinite_indices && matching || return _SWFloat(Inf)
            elseif a!=b
                scale = max(abs(a), abs(b))
                change = max(change, abs(a-b)/scale)
            end
        end
    end
    return change
end

function _precision_wave_output(result, bits; scaled = false)
    scaled || return map(channel -> _round_to_bigfloat.(channel, bits), result)
    stored_bits = maximum(channel -> maximum(v -> Base.precision(real(v)), channel; init = bits),
        values(result); init = bits)
    return _with_swprecision(max(bits+32, stored_bits)) do
        map(result) do channel
            exponents = [iszero(v) || !isfinite(v) ? 0 : floor(Int, log10(abs(v))) for v in channel]
            mantissas = [_round_to_bigfloat(v/_SWFloat(10)^e, bits) for (v, e) in zip(channel, exponents)]
            for i in eachindex(mantissas)
                value = mantissas[i] isa Real ? _SWFloat(mantissas[i]) : _SWFloat_complex(mantissas[i])
                if isfinite(value) && abs(value)>=10
                    mantissas[i] = _round_to_bigfloat(channel[i]/_SWFloat(10)^(exponents[i]+1), bits)
                    exponents[i] += 1
                end
            end
            (; mantissa = mantissas, exponent = exponents)
        end
    end
end

_precision_input_bits(x::Real) = 0
_precision_input_bits(x::BigFloat) = Int(ccall((:mpfr_min_prec, Base.MPFR.libmpfr), Clong, (Ref{BigFloat},), x))
_precision_input_bits(x::_SWFloat) = _precision_input_bits(x.value)
_precision_input_bits(x::Complex) = max(_precision_input_bits(real(x)), _precision_input_bits(imag(x)))
_precision_input_bits(x::Rational) = max(ndigits(numerator(x); base = 2), ndigits(denominator(x); base = 2))

function _angular_values_at_precision(eigenpair, points, kind, normalize, order;
        sensitivity = false, rtol, max_terms, step_fraction = 1)
    m, n, c, spheroid = eigenpair.m, eigenpair.n, eigenpair.c, eigenpair.spheroid
    coordinates = [x isa Real ? _SWFloat(x) : _SWFloat_complex(x) for x in points]
    if kind==1
        phase = _coefficient_phase_at_precision(eigenpair)
        scale = phase*(normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(m, n, _SWFloat)))
        value = _evaluate_coefficient_vector(eigenpair, scale .* eigenpair.v, coordinates;
            derivatives = order)
        sensitivity || return value
        tangent = _evaluate_coefficient_vector(eigenpair, scale .* eigenpair.dv, coordinates)
        return (; value..., dvalue_dc = tangent.value, dderivative_dc = tangent.derivative)
    end
    normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1"))
    data = _qs_initial_data(m, n, c, spheroid, :mpfr; eigenpair, rtol,
        tail_rtol = rtol, max_terms)
    V = c isa Real && all(x -> x isa Real, points) ? _SWFloat : Complex{_SWFloat}
    arrays = [zeros(V, length(points)) for _ in 0:order]
    tangents, dtangents = zeros(V, length(points)), zeros(V, length(points))
    equation = (; eigenpair..., spheroid = :prolate, q = data.g, dq = data.dg)
    initial = sensitivity ? (data.y, data.dy, data.z, data.dz) : (data.y, data.dy)
    for (i, x) in enumerate(coordinates)
        if isreal(x) && abs(real(x))==1
            parity = real(x)<0 && iseven(n-m) ? -1 : 1
            divergent = _directed_infinity(data.amplitude, V<:Real, _SWFloat)
            for j in 0:order
                arrays[j+1][i] = parity*real(x)^j*divergent
            end
            tangents[i] = dtangents[i] = V(NaN)
            continue
        end
        state = if iszero(x)
            initial
        else
            _coordinate_continuation(equation, V.(initial), zero(V), [V(x)];
                step_fraction, rtol, max_terms)
        end
        y, dy = state[1:2]
        arrays[1][i], arrays[2][i] = y, dy
        if sensitivity
            tangents[i], dtangents[i] = state[3:4]
        end
        order==1 && continue
        a, b, d = _wave_equation(m, c, x, eigenpair.lambda, spheroid, :angular)
        ddy = -(b*dy+d*y)/a
        arrays[3][i] = ddy
        if order>=3
            dprime = -2data.g*x+m^2*b/a^2
            third = -(2b*ddy+(d-2)*dy+dprime*y)/a
            arrays[4][i] = third
            if order==4
                dsecond = -2data.g+m^2*(-2/a^2-2b^2/a^3)
                arrays[5][i] = -(3b*third+(d-6)*ddy+2dprime*dy+dsecond*y)/a
            end
        end
    end
    names = (:value, :derivative, :second_derivative, :third_derivative, :fourth_derivative)
    result = NamedTuple{names[1:(order+1)]}(Tuple(arrays))
    return sensitivity ? (; result..., dvalue_dc = tangents, dderivative_dc = dtangents) : result
end

function _refine_angular_values_at_precision(m, n, c, points, spheroid, precision,
        bits, max_terms, kind, normalize, order; sensitivity = false, max_refinements = 4)
    previous = nothing
    min_terms = 0
    input_bits = max(_precision_input_bits(c), maximum(_precision_input_bits, points; init = 0))
    requested = max(bits+32, input_bits)
    endpoints = findall(x -> isreal(x) && abs(real(x))==1, points)
    for attempt in 1:max_refinements
        eigenpair = _solve_coefficient_eigenpair_at_precision(m, n, c, spheroid, precision,
            requested, nothing, max_terms, sensitivity || kind==2; min_terms)
        result, change = try
            _with_swprecision(Base.precision(real(eigenpair.lambda))) do
                tolerance = ldexp(one(_SWFloat), -requested-8)
                values = _angular_values_at_precision(eigenpair, points, kind, normalize, order;
                    sensitivity, rtol = tolerance, max_terms, step_fraction = 1//2^(attempt-1))
                for name in keys(values), i in eachindex(getproperty(values, name))
                    i in endpoints || isfinite(getproperty(values, name)[i]) ||
                        error("Angular evaluation produced a nonfinite interior value")
                end
                delta = previous===nothing ? _SWFloat(Inf) :
                    _precision_wave_change(previous, values; nonfinite_indices = endpoints)
                values, delta
            end
        catch exception
            refinable = exception isa ErrorException &&
                (occursin("Qs joining sums did not converge", exception.msg) ||
                 occursin("Differentiated radial equation did not converge", exception.msg))
            refinable || rethrow()
            nothing, _SWFloat(Inf)
        end
        if change<=_with_swprecision(() -> ldexp(one(_SWFloat), -bits-8), bits+32)
            return (; result, eigenpair, relative_change = change)
        end
        result===nothing || (previous = result)
        min_terms = min(max_terms, eigenpair.terms+max(16, eigenpair.terms÷2))
        requested += max(32, bits÷2)
    end
    error("Angular values did not converge at the requested precision within max_terms=$max_terms")
end

function _smn_at_precision(m, n, c, points; spheroid = :prolate, precision = :double,
        precision_bits, max_terms = 512, normalize = false, kind = 1, derivatives = 1,
        second_derivative = false, scaled = false, logderivative = false)
    bits = _validate_precision_bits(precision_bits)
    _validate_coefficient_options(m, n, c, spheroid, precision, nothing, max_terms, :orthonormal)
    if eltype(points)<:Real
        _validate_wave_arguments(m, n, c, points, spheroid, precision, :angular; kind)
    else
        _validate_complex_coordinates(m, n, c, points, spheroid, precision, :angular; kind)
    end
    kind==2 && normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1"))
    derivatives in 1:4 || throw(ArgumentError("derivatives must be in 1:4"))
    kind in (1, 2) || throw(ArgumentError("angular kind must be 1 or 2"))
    order = max(derivatives, second_derivative ? 2 : 1)
    data = _refine_angular_values_at_precision(m, n, c, points, spheroid, precision,
        bits, max_terms, kind, normalize, order)
    return _with_swprecision(Base.precision(real(data.eigenpair.lambda))) do
        output = _precision_wave_output(data.result, bits; scaled)
        if logderivative
            ratios = map(eachindex(points)) do i
                x, y, dy = points[i], data.result.value[i], data.result.derivative[i]
                if isreal(x) && abs(real(x))==1 && (kind==2 || m>0)
                    oftype(y, _SWFloat(copysign(Inf, real(x)))*(kind==1 ? -1 : 1))
                else
                    iszero(y) ? oftype(y, NaN) : dy/y
                end
            end
            return (; output..., logderivative = _round_to_bigfloat.(ratios, bits))
        end
        return output
    end
end

function _angular_parameter_derivatives_at_precision(m, n, c, points;
        spheroid = :prolate, precision = :double, precision_bits, max_terms = 512,
        normalize = false, kind = 1, diagnostics = false)
    bits = _validate_precision_bits(precision_bits)
    _validate_coefficient_options(m, n, c, spheroid, precision, nothing, max_terms, :orthonormal)
    if eltype(points)<:Real
        _validate_wave_arguments(m, n, c, points, spheroid, precision, :angular; kind)
    else
        _validate_complex_coordinates(m, n, c, points, spheroid, precision, :angular; kind)
    end
    kind==2 && normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1"))
    data = _refine_angular_values_at_precision(m, n, c, points, spheroid, precision,
        bits, max_terms, kind, normalize, 1; sensitivity = true)
    return _with_swprecision(Base.precision(real(data.eigenpair.lambda))) do
        value = _round_to_bigfloat.(data.result.dvalue_dc, bits)
        derivative = _round_to_bigfloat.(data.result.dderivative_dc, bits)
        finite = all(isfinite, value) && all(isfinite, derivative)
        metadata = (; method = kind==1 ? :coefficients : :differentiated_equation,
            step_used = nothing, relative_change_when_halving_step = nothing,
            coefficient_relative_change = _round_to_bigfloat(data.relative_change, bits),
            coefficient_tail = _round_to_bigfloat(data.eigenpair.tail, bits),
            eigenproblem_residual = _round_to_bigfloat(data.eigenpair.residual, bits),
            sensitivity_residual = _round_to_bigfloat(data.eigenpair.derivative_residual, bits),
            eigenvalue_condition = _round_to_bigfloat(data.eigenpair.eigenvalue_condition, bits),
            finite_flag = finite,
            conditioning_flag = !finite ? :poor : data.eigenpair.eigenvalue_condition>100 ? :warning : :good,
            suggested_action = finite ? :accept : :unresolved)
        if c isa Real
            result = (; dvalue_dc = value, dderivative_dc = derivative)
            return diagnostics ? (; result..., metadata_value = metadata, metadata_derivative = metadata) : result
        end
        imaginary_value = _round_to_bigfloat.(im .* data.result.dvalue_dc, bits)
        imaginary_derivative = _round_to_bigfloat.(im .* data.result.dderivative_dc, bits)
        result = (; dvalue_dcreal = value, dvalue_dcimag = imaginary_value,
            dderivative_dcreal = derivative, dderivative_dcimag = imaginary_derivative)
        return diagnostics ? (; result..., metadata_value_dcreal = metadata,
            metadata_value_dcimag = metadata, metadata_derivative_dcreal = metadata,
            metadata_derivative_dcimag = metadata) : result
    end
end
