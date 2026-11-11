function _precision_stencil_metadata(first, second, step, bits, rtol, atol, order)
    finite = all(isfinite, first) && all(isfinite, second)
    relative = finite ? maximum(abs.(second .- first) ./ max.(abs.(second), _SWFloat(atol))) : _SWFloat(Inf)
    condition = !finite ? :poor : relative<=10rtol ? :good : relative<=1000rtol ? :warning : :poor
    (; method = order==1 ? :finite_difference : :differenced_sensitivity,
        step_used = _round_to_bigfloat(step, bits),
        relative_change_when_halving_step = _round_to_bigfloat(relative, bits),
        finite_flag = finite, conditioning_flag = condition,
        suggested_action = condition===:poor ? :retry_smaller_h : :accept)
end

function _precision_parameter_stencil(evaluate, c, bits, h, order, adaptive, rtol, atol;
        nonnegative = false, nonzero = false)
    step, workbits = _with_swprecision(max(bits+64, _precision_input_bits(c)+32)) do
        parameter = c isa Real ? _SWFloat(c) : _SWFloat_complex(c)
        scale = nonzero && !iszero(parameter) ? abs(parameter) : max(one(_SWFloat), abs(parameter))
        step = h===nothing ? ldexp(scale, -cld(bits+24, 4)) : _SWFloat(h)
        isfinite(step) && step>0 || throw(ArgumentError("h must be finite and positive"))
        if !iszero(parameter) && (nonnegative || nonzero)
            h===nothing && (step=min(step, abs(parameter)/4))
            order*step<abs(parameter) || throw(ArgumentError("$(order)h must be smaller than abs(c)"))
        end
        lost = max(0, ceil(Int, log2(max(one(_SWFloat), abs(parameter))/step)))
        step, max(bits+lost+64, _precision_input_bits(c)+32, _precision_input_bits(h===nothing ? 0 : h)+32)
    end
    return _with_swprecision(workbits) do
        parameter = c isa Real ? _SWFloat(c) : _SWFloat_complex(c)
        cache = Dict{Any, Any}()
        sample(z) = get!(cache, z) do
            map(x -> x isa Real ? _SWFloat(x) : _SWFloat_complex(x), evaluate(z, workbits))
        end
        function estimate(s, direction)
            boundary = nonnegative && iszero(parameter)
            offsets, weights, denominator = if order==1
                boundary ? ((0, 1, 2), (-3, 4, -1), 2s) : ((-1, 1), (-1, 1), 2s)
            else
                boundary ? ((0, 1, 2, 3, 4), (-25, 48, -36, 16, -3), 12s) :
                    ((-2, -1, 1, 2), (1, -8, 8, -1), 12s)
            end
            grid = [parameter+direction*k*s for k in offsets]
            length(unique(grid))==length(grid) || throw(ArgumentError("h is below parameter resolution"))
            center = sample(parameter)
            sum(w .* (sample(z) .- center) for (w, z) in zip(weights, grid)) ./ denominator
        end
        function differentiate(direction)
            first = estimate(step, direction)
            for attempt in 1:(h===nothing ? 5 : 2)
                current_step = step/2^attempt
                second = estimate(current_step, direction)
                metadata = _precision_stencil_metadata(first, second,
                    adaptive || h===nothing ? current_step : step, bits, rtol, atol, order)
                if h!==nothing
                    !adaptive && return first, metadata
                    (metadata.conditioning_flag!==:poor || attempt==2) && return second, metadata
                else
                    change = _precision_wave_change((; value = vec(first)), (; value = vec(second));
                        nonfinite_indices = eachindex(second))
                    change<=ldexp(one(_SWFloat), -bits-4) && return second, metadata
                end
                first = second
            end
            error("Parameter derivative did not converge under step refinement")
        end
        real_value, real_metadata = differentiate(1)
        if c isa Complex && order==1
            imag_value, imag_metadata = differentiate(im)
            return (; real_value, real_metadata, imag_value, imag_metadata)
        end
        (; real_value, real_metadata)
    end
end

function _precision_endpoint_parameter(m, n, c, points, spheroid, precision, target,
        normalize, normalization, bits, max_terms, sensitivity)
    pair = _solve_coefficient_eigenpair_at_precision(m, n, c, spheroid, precision,
        bits, nothing, max_terms, true)
    _with_swprecision(Base.precision(real(pair.lambda))) do
        if target===:angular
            scale = _coefficient_scale_at_precision(pair, normalize)
            coefficients = sensitivity ? pair.dv : pair.v
            regular = _evaluate_coefficient_vector(pair, scale .* coefficients, real.(points); regular_factor = true)
            return -real.(points) .* regular.value
        end
        value = if iszero(c)
            sensitivity ? zero(_SWFloat) : one(_SWFloat)/prod(_SWFloat(2j+1) for j in 0:n)
        else
            state = _radial_expansion_data(pair, one(_SWFloat), 1; regular_factor = true).state
            if normalization===:static
                sensitivity ? pair.c^(-n)*(state[3]-n*state[1]/pair.c) : pair.c^(-n)*state[1]
            else
                state[sensitivity ? 3 : 1]
            end
        end
        fill(value, length(points))
    end
end

function _wave_jacobian_at_precision(m, n, c, points, spheroid, precision, target, kind,
        normalize, normalization, bits, max_terms, order, h, diagnostics, adaptive, rtol, atol)
    bits = _validate_precision_bits(bits)
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    analytic(z, p; diagnostics = false) = target===:angular ?
        _angular_parameter_derivatives_at_precision(m, n, z, points; spheroid, precision,
            precision_bits = p, max_terms, kind, normalize, diagnostics) :
        _radial_parameter_derivatives_at_precision(m, n, z, points; spheroid, precision,
            precision_bits = p, max_terms, kind, normalization, diagnostics)
    order==1 && h===nothing && return analytic(c, bits; diagnostics)
    endpoints = m==1 && kind==1 ? findall(x -> isreal(x) &&
        (target===:angular ? abs(x)==1 : spheroid===:prolate && x==1), points) : Int[]
    evaluate(z, p) = begin
        result = if order==2
            sensitivity = analytic(z, p)
            z isa Real ? (sensitivity.dvalue_dc, sensitivity.dderivative_dc) :
                (sensitivity.dvalue_dcreal, sensitivity.dderivative_dcreal)
        else
            wave = target===:angular ?
                smn(m, n, z, points; spheroid, precision, kind, normalize, precision_bits = p, max_terms) :
                rmn(m, n, z, points; spheroid, precision, kind, normalization, precision_bits = p, max_terms)
            (wave.value, wave.derivative)
        end
        values = hcat(result...)
        if !isempty(endpoints)
            values[endpoints, 2] = _precision_endpoint_parameter(m, n, z, points[endpoints],
                spheroid, precision, target, normalize, normalization, p, max_terms, order==2)
        end
        values
    end
    data = _precision_parameter_stencil(evaluate, c, bits, h, order, adaptive, rtol, atol;
        nonzero = target===:radial && normalization===:standard)
    return _with_swprecision(bits+64) do
        restore(values) = begin
            result = copy(values)
            for i in endpoints
                result[i, 2] = _directed_infinity(result[i, 2], c isa Real && eltype(points)<:Real, _SWFloat)
            end
            _round_to_bigfloat.(result, bits)
        end
        result = restore(data.real_value)
        value, derivative = result[:, 1], result[:, 2]
        metadata = data.real_metadata
        if any(isinf, derivative)
            metadata = (; metadata..., finite_flag = false, conditioning_flag = :singular, suggested_action = :singular_limit)
        end
        if order==2
            output = c isa Real ? (; d2value_dc2 = value, d2derivative_dc2 = derivative) :
                (; d2value_dcreal2 = value, d2value_dcreal_dcimag = restore(im .* data.real_value)[:, 1],
                    d2value_dcimag2 = restore(-data.real_value)[:, 1],
                    d2derivative_dcreal2 = derivative, d2derivative_dcreal_dcimag = restore(im .* data.real_value)[:, 2],
                    d2derivative_dcimag2 = restore(-data.real_value)[:, 2])
            return diagnostics ? (; output..., metadata) : output
        elseif c isa Real
            output = (; dvalue_dc = value, dderivative_dc = derivative)
            return diagnostics ? (; output..., metadata_value = metadata, metadata_derivative = metadata) : output
        end
        imaginary = restore(data.imag_value)
        output = (; dvalue_dcreal = value, dderivative_dcreal = derivative,
            dvalue_dcimag = imaginary[:, 1], dderivative_dcimag = imaginary[:, 2])
        diagnostics ? (; output..., metadata_value_dcreal = metadata, metadata_derivative_dcreal = metadata,
            metadata_value_dcimag = data.imag_metadata, metadata_derivative_dcimag = data.imag_metadata) : output
    end
end

function _eigen_difference_at_precision(m, n, c, spheroid, precision, operator, form,
        bits, max_terms, order, h, diagnostics, adaptive, rtol, atol)
    _validate_coefficient_options(m, n, c, spheroid, precision, nothing, max_terms, :orthonormal)
    operator===:separation || _validate_integral_precision_arguments(m, n, c, spheroid, operator, form)
    operator!==:separation && iszero(c) && form===:log &&
        throw(DomainError(c, "log concentration is singular at c=0"))
    evaluate(z, p) = begin
        value = order==1 ? eigenvalue(m, n, z; spheroid, precision, operator, form, precision_bits = p, max_terms) :
            jacobian_eigen(m, n, z; spheroid, precision, operator, form, precision_bits = p, max_terms)
        reshape([order==2 && z isa Complex ? value.d_dcreal : value], 1, 1)
    end
    data = _precision_parameter_stencil(evaluate, c, bits, h, order, adaptive, rtol, atol;
        nonnegative = operator!==:separation)
    _with_swprecision(bits+64) do
        value = only(data.real_value)
        rounded = _round_to_bigfloat(value, bits)
        metadata = data.real_metadata
        c isa Real && return diagnostics ? (; derivative = rounded, metadata) : rounded
        if order==2
            result = (; d2_dcreal2 = rounded, d2_dcreal_dcimag = _round_to_bigfloat(im*value, bits),
                d2_dcimag2 = _round_to_bigfloat(-value, bits))
            return diagnostics ? (; result..., metadata) : result
        end
        result = (; d_dcreal = rounded, d_dcimag = _round_to_bigfloat(only(data.imag_value), bits))
        diagnostics ? (; result..., metadata_dcreal = metadata, metadata_dcimag = data.imag_metadata) : result
    end
end
