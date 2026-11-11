function _radial_states_at_precision(eigenpair, points, kind, normalization, tolerance, step_fraction)
    m, n, c, spheroid = eigenpair.m, eigenpair.n, eigenpair.c, eigenpair.spheroid
    if iszero(c)
        normalization === :static || throw(DomainError(c, "standard radial normalization requires nonzero c"))
        sigma = spheroid === :prolate ? 1 : -1
        anchor = _SWFloat(2)
        initial_y, initial_dy = kind==1 ? _static_regular_radial(m, n, anchor, sigma) :
                               _static_irregular_radial(m, n, anchor, sigma)
        initial = complex.((initial_y, initial_dy, zero(initial_y), zero(initial_dy)))
        states = map(points) do point
            x = point isa Real ? _SWFloat(point) : _SWFloat_complex(point)
            spheroid===:prolate && x==1 && kind!=1 && return ntuple(_ -> complex(_SWFloat(NaN)), 4)
            if point isa Real && (kind==1 || x>=2)
                y, dy = kind==1 ? _static_regular_radial(m, n, x, sigma) :
                        _static_irregular_radial(m, n, x, sigma)
                return complex.((y, dy, zero(y), zero(dy)))
            end
            _coordinate_continuation(eigenpair, initial, complex(anchor),
                _coordinate_path(complex(x), :radial, spheroid); step_fraction, rtol = tolerance)
        end
        return (; states, tail = zero(_SWFloat))
    end
    if all(x -> x isa Real, points)
        data = _radial_expansion_batch(eigenpair, _SWFloat.(points), kind;
            rtol = tolerance, step_fraction)
        states, tail = data.states, data.tail
    else
        anchor = _SWFloat(2)
        data = _radial_expansion_data(eigenpair, anchor, kind)
        tail = data.relative_tail
        states = map(points) do point
            if isreal(point) && real(point)>= (spheroid===:prolate ? 1 : 0)
                real_data = _radial_expansion_batch(eigenpair, [_SWFloat(real(point))], kind;
                    rtol = tolerance, step_fraction)
                tail = max(tail, real_data.tail)
                return complex.(only(real_data.states))
            end
            z = _SWFloat_complex(point)
            _coordinate_continuation(eigenpair, complex.(data.state), complex(anchor),
                _coordinate_path(z, :radial, spheroid); step_fraction, rtol = tolerance)
        end
    end
    if normalization===:static
        power = kind==1 ? -n : n+1
        factor = c^power
        states = map(states) do state
            y, dy, tangent, dtangent = state
            (factor*y, factor*dy, factor*(tangent+power*y/c), factor*(dtangent+power*dy/c))
        end
    end
    (; states, tail)
end

function _radial_coordinate_derivatives_at_precision(eigenpair, points, states, kind, normalization, order)
    m, n, c, spheroid = eigenpair.m, eigenpair.n, eigenpair.c, eigenpair.spheroid
    arrays = [Complex{_SWFloat}[state[k] for state in states] for k in 1:2]
    append!(arrays, [zeros(Complex{_SWFloat}, length(points)) for _ in 2:order])
    sigma = spheroid===:prolate ? 1 : -1
    for (i, point) in enumerate(points)
        x = point isa Real ? _SWFloat(point) : _SWFloat_complex(point)
        order==1 && continue
        if spheroid===:prolate && x==1
            if kind!=1
                for k in 3:(order+1)
                    arrays[k][i] = complex(_SWFloat(NaN))
                end
                continue
            end
            amplitude = if iszero(c)
                anchor = 1+ldexp(one(_SWFloat), -8)
                y, _ = _static_regular_radial(m, n, anchor, sigma)
                polynomial = _angular_endpoint_coefficients(m, eigenpair.lambda, c^2, anchor-1)
                y/((anchor^2-1)^(m//2)*_angular_endpoint_polynomial(polynomial, 1-anchor).value)
            else
                value = _radial_expansion_data(eigenpair, one(_SWFloat), 1; regular_factor = true).state[1]
                normalization===:static ? value*c^(-n) : value
            end
            polynomial = _angular_endpoint_coefficients(m, eigenpair.lambda, c^2, zero(_SWFloat))
            for k in 2:order
                arrays[k+1][i] = _endpoint_derivative(m, k, amplitude, polynomial, one(_SWFloat),
                    :radial, c isa Real, _SWFloat)
            end
            continue
        end
        y, dy = arrays[1][i], arrays[2][i]
        a, b, d = _wave_equation(m, c, x, eigenpair.lambda, spheroid, :radial)
        ddy = -(b*dy+d*y)/a
        arrays[3][i] = ddy
        order==2 && continue
        dprime = 2c^2*x+sigma*m^2*b/a^2
        dsecond = 2c^2+sigma*m^2*(2/a^2-2b^2/a^3)
        third = -(2b*ddy+(2+d)*dy+dprime*y)/a
        arrays[4][i] = third
        order==4 && (arrays[5][i] = -(3b*third+(6+d)*ddy+2dprime*dy+dsecond*y)/a)
    end
    names = (:value, :derivative, :second_derivative, :third_derivative, :fourth_derivative)
    NamedTuple{names[1:(order+1)]}(Tuple(arrays))
end

function _refine_radial_values_at_precision(m, n, c, points, spheroid, precision, kind,
        normalization, bits, max_terms, order, sensitivity; max_refinements = 5, local_mode = nothing)
    input_bits = max(_precision_input_bits(c), maximum(_precision_input_bits, points; init = 0))
    workbits = max(bits+96, input_bits+64, bits+64+ceil(Int, 8abs(c)+4abs(imag(c))*maximum(abs, points)))
    previous = nothing
    minimum = min(max_terms, max(32, cld(bits, 2)+32))
    endpoints = findall(x -> spheroid===:prolate && isone(x), points)
    for attempt in 1:max_refinements
        result, stable = _with_swprecision(workbits) do
            tolerance = ldexp(one(_SWFloat), -bits-24)
            eigenpair = _solve_coefficient_eigenpair_at_precision(m, n, c, spheroid, precision,
                workbits-48, nothing, max_terms, true; min_terms = minimum, local_mode)
            _with_swprecision(Base.precision(real(eigenpair.lambda))) do
                data = _radial_states_at_precision(eigenpair, points, kind, normalization, tolerance, 1//attempt)
                result = sensitivity ?
                         (; dvalue_dc = Complex{_SWFloat}[s[3] for s in data.states],
                    dderivative_dc = Complex{_SWFloat}[s[4] for s in data.states]) :
                         _radial_coordinate_derivatives_at_precision(eigenpair, points, data.states, kind, normalization, order)
                stable = data.tail<=tolerance && previous!==nothing &&
                         _precision_wave_change(previous, result; nonfinite_indices = endpoints)<=ldexp(one(_SWFloat), -bits-8)
                result, stable
            end
        end
        stable && return result
        previous = result
        workbits += max(64, bits÷2)
        minimum = min(max_terms, minimum+max(32, minimum÷2))
    end
    error("Requested radial precision did not converge within max_terms=$max_terms")
end

function _rmn_at_precision(m, n, c, points; spheroid, precision, precision_bits, max_terms = 512,
        kind = 1, normalization = :standard, derivatives = 1, second_derivative = false,
        scaled = false, logderivative = false, local_mode = nothing)
    bits = _validate_precision_bits(precision_bits)
    _validate_coefficient_options(m, n, c, spheroid, precision, nothing, max_terms, :orthonormal)
    derivatives in 1:4 || throw(ArgumentError("derivatives must be in 1:4"))
    if eltype(points)<:Real
        _validate_wave_arguments(m, n, c, points, spheroid, precision, :radial; kind, normalization)
    else
        _validate_complex_coordinates(m, n, c, points, spheroid, precision, :radial; kind, normalization)
    end
    order = max(derivatives, second_derivative ? 2 : 1)
    result = _refine_radial_values_at_precision(m, n, c, points, spheroid, precision, kind,
        normalization, bits, max_terms, order, false; local_mode)
    output = _with_swprecision(max(bits+64, Base.precision(real(first(result.value))))) do
        _precision_wave_output(result, bits; scaled)
    end
    if logderivative
        ratio = _with_swprecision(bits+64) do
            [_round_to_bigfloat(isreal(x) ?
                _wave_logderivative(v, d, real(x), m, spheroid, :radial, kind) :
                iszero(v) ? oftype(v, NaN) : d/v, bits)
             for (x, v, d) in zip(points, result.value, result.derivative)]
        end
        output = (; output..., logderivative = ratio)
    end
    return output
end

function _radial_parameter_derivatives_at_precision(m, n, c, points; spheroid, precision,
        precision_bits, max_terms = 512, kind = 1, normalization = :standard, diagnostics = false,
        local_mode = nothing)
    bits = _validate_precision_bits(precision_bits)
    _validate_coefficient_options(m, n, c, spheroid, precision, nothing, max_terms, :orthonormal)
    if eltype(points)<:Real
        _validate_wave_arguments(m, n, c, points, spheroid, precision, :radial; kind, normalization)
    else
        _validate_complex_coordinates(m, n, c, points, spheroid, precision, :radial; kind, normalization)
    end
    result = _refine_radial_values_at_precision(m, n, c, points, spheroid, precision, kind,
        normalization, bits, max_terms, 1, true; local_mode)
    output = _precision_wave_output(result, bits)
    if c isa Complex
        stored_bits = maximum(channel -> maximum(z -> Base.precision(real(z)), channel), values(result))
        output = _with_swprecision(stored_bits) do
            rotate(z) = _round_to_bigfloat(complex(-imag(z), real(z)), bits)
            (; dvalue_dcreal = output.dvalue_dc, dvalue_dcimag = rotate.(result.dvalue_dc),
                dderivative_dcreal = output.dderivative_dc, dderivative_dcimag = rotate.(result.dderivative_dc))
        end
    end
    diagnostics || return output
    finite = all(channel -> all(isfinite, channel), values(output))
    metadata = _precision_derivative_metadata(finite ? 0 : NaN, :coefficients)
    if c isa Real
        return (; output..., metadata_value = metadata, metadata_derivative = metadata)
    end
    (; output..., metadata_value_dcreal = metadata, metadata_value_dcimag = metadata,
        metadata_derivative_dcreal = metadata, metadata_derivative_dcimag = metadata)
end
