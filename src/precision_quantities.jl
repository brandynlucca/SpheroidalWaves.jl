function _refine_quantity_at_precision(calculate, bits, rtol = nothing; valid = _ -> true)
    workbits = _with_swprecision(bits+64) do
        rtol === nothing ? bits+96 : max(bits+96, ceil(Int, -log2(_SWFloat(rtol)))+96)
    end
    previous = nothing
    for attempt in 1:5
        result, stable = _with_swprecision(workbits) do
            result = calculate(workbits-48)
            tolerance = ldexp(one(_SWFloat), -bits-8)
            rtol === nothing || (tolerance = min(tolerance, _SWFloat(rtol)/256))
            stable = valid(result) && all(isfinite, values(result)) && previous !== nothing &&
                     all(abs(a-b)<=tolerance*abs(a)
                         for (a, b) in zip(values(result), values(previous)))
            result, stable
        end
        stable && return result
        previous = result
        workbits += max(96, bits÷2)
    end
    error("Requested quantity did not converge under working-precision refinement")
end

function _coefficient_scale_at_precision(eigenpair, normalize)
    phase = _coefficient_phase_at_precision(eigenpair)
    phase*(normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(eigenpair.m, eigenpair.n, _SWFloat)))
end

function _amn_at_precision(m, n, c, spheroid, precision, bits, rtol, max_terms)
    bits = _validate_precision_bits(bits)
    p = abs(m)
    _validate_coefficient_options(p, n, c, spheroid, precision, rtol, max_terms, :orthonormal)
    result = _refine_quantity_at_precision(bits, rtol) do workbits
        eigenpair = _solve_coefficient_eigenpair_at_precision(
            p, n, c, spheroid, precision, workbits, rtol, max_terms, true)
        _with_swprecision(Base.precision(real(eigenpair.lambda))) do
            scale = _coefficient_scale_at_precision(eigenpair, false)
            sums = _coefficient_joining_sums(eigenpair, scale; extend = m>=0)
            value = m<0 ? sums.weighted/prod(_SWFloat(k)
                for k in (n-p+1):(n+p); init = one(_SWFloat)) : sums.joining
            (; value)
        end
    end
    _round_to_bigfloat(result.value, bits)
end

function _kmn_at_precision(m, n, c, spheroid, precision, normalize, bits, rtol, max_terms)
    bits = _validate_precision_bits(bits)
    _validate_coefficient_options(m, n, c, spheroid, precision, rtol, max_terms, :orthonormal)
    result = _refine_quantity_at_precision(bits, rtol) do workbits
        eigenpair = _solve_coefficient_eigenpair_at_precision(
            m, n, c, spheroid, precision, workbits, rtol, max_terms, false)
        _with_swprecision(Base.precision(real(eigenpair.lambda))) do
            (; value = _joining_factor_value(eigenpair, _coefficient_scale_at_precision(eigenpair, normalize)))
        end
    end
    _round_to_bigfloat(result.value, bits)
end

function _precision_derivative_metadata(value, method)
    finite = isfinite(value)
    (; method, step_used = nothing, relative_change_when_halving_step = nothing,
        finite_flag = finite, conditioning_flag = finite ? :good : :singular,
        suggested_action = finite ? :accept : :singular_limit)
end

function _radial_wronskian_at_precision(m, n, c, points, spheroid, precision, form,
        normalization, bits, max_terms)
    bits = _validate_precision_bits(bits)
    _with_swprecision(bits+64) do
        first = rmn(m, n, c, points; spheroid, precision, precision_bits = bits+64,
            max_terms, kind = 1, normalization)
        second = rmn(m, n, c, points; spheroid, precision, precision_bits = bits+64,
            max_terms, kind = 2, normalization)
        values = _SWFloat_complex.(first.value) .* _SWFloat_complex.(second.derivative) .-
                 _SWFloat_complex.(first.derivative) .* _SWFloat_complex.(second.value)
        if form!==:raw
            coordinates = [x isa Real ? _SWFloat(x) : _SWFloat_complex(x) for x in points]
            factor = coordinates .^ 2 .- (spheroid===:prolate ? 1 : -1)
            if normalization===:standard
                parameter = c isa Real ? _SWFloat(c) : _SWFloat_complex(c)
                factor = factor .* parameter
            end
            values .*= factor
            form===:error && return _round_to_bigfloat.(abs.(values .- 1), bits)
        end
        _round_to_bigfloat.(values, bits)
    end
end

function _eigen_jacobian_at_precision(m, n, c, spheroid, precision, operator, form,
        bits, max_terms, order, h, diagnostics, adaptive = true, rtol = 1e-6, atol = 1e-10;
        local_mode = nothing)
    bits = _validate_precision_bits(bits)
    h === nothing || return _eigen_difference_at_precision(m, n, c, spheroid, precision,
        operator, form, bits, max_terms, order, h, diagnostics, adaptive, rtol, atol)
    if operator !== :separation
        value = _integral_eigen_jacobian_at_precision(
            m, n, c, spheroid, precision, operator, form, bits, max_terms, order)
        metadata = _precision_derivative_metadata(value, iszero(c) ? :right_limit : :integral_identity)
    else
        _validate_coefficient_options(m, n, c, spheroid, precision, nothing, max_terms, :orthonormal)
        result = _refine_quantity_at_precision(bits) do workbits
            eigenpair = _solve_coefficient_eigenpair_at_precision(
                m, n, c, spheroid, precision, workbits, nothing, max_terms, true; local_mode)
            _with_swprecision(Base.precision(real(eigenpair.lambda))) do
                value = if order==1
                    eigenpair.dlambda
                else
                    _, _, dH = _coefficient_matrix(m, n, eigenpair.c, eigenpair.terms, spheroid)
                    _, _, ddH = _coefficient_matrix(m, n, one(eigenpair.c), eigenpair.terms, spheroid)
                    sum(eigenpair.v .* (ddH*eigenpair.v)) + 2sum(eigenpair.v .* (dH*eigenpair.dv))
                end
                (; value)
            end
        end
        value = _round_to_bigfloat(result.value, bits)
        metadata = _precision_derivative_metadata(value, :coefficients)
    end
    c isa Real && return diagnostics ? (; derivative = value, metadata) : value
    imaginary, negative = _with_swprecision(Base.precision(real(result.value))) do
        unrounded = result.value
        (_round_to_bigfloat(complex(-imag(unrounded), real(unrounded)), bits),
            _round_to_bigfloat(-unrounded, bits))
    end
    if order==2
        result = (; d2_dcreal2 = value, d2_dcreal_dcimag = imaginary, d2_dcimag2 = negative)
        return diagnostics ? (; result..., metadata) : result
    end
    result = (; d_dcreal = value, d_dcimag = imaginary)
    diagnostics ? (; result..., metadata_dcreal = metadata, metadata_dcimag = metadata) : result
end

function _find_c_at_precision(m, n, target, bracket, spheroid, precision, operator, form,
        bits, atol, rtol, maxiter, use_jacobian, max_terms)
    bits = _validate_precision_bits(bits)
    maxiter>0 || throw(ArgumentError("maxiter must be positive"))
    for tolerance in (atol, rtol)
        tolerance === nothing || (isfinite(tolerance) && tolerance>0) ||
            throw(ArgumentError("atol and rtol must be finite and positive"))
    end
    workbits = _with_swprecision(bits+96) do
        maximum(t === nothing ? bits+96 : max(bits+96, ceil(Int, -log2(_SWFloat(t)))+96)
            for t in (atol, rtol))
    end
    _with_swprecision(workbits) do
        a, b, y = _SWFloat(bracket[1]), _SWFloat(bracket[2]), _SWFloat(target)
        isfinite(a) && isfinite(b) && a<b || throw(ArgumentError("bracket must be finite and increasing"))
        tolerance = ldexp(one(_SWFloat), -bits-8)
        absolute = atol === nothing ? tolerance : min(tolerance, _SWFloat(atol))
        relative = rtol === nothing ? tolerance : min(tolerance, _SWFloat(rtol))
        output(converged, c, residual, iterations, method) =
            (; converged, c = _round_to_bigfloat(c, bits), residual = _round_to_bigfloat(residual, bits),
                iterations, bracket = (_round_to_bigfloat(a, bits), _round_to_bigfloat(b, bits)), method)
        if operator===:concentration
            _validate_integral_precision_arguments(m, n, a, spheroid, :concentration, form)
            _validate_integral_precision_arguments(m, n, b, spheroid, :concentration, form)
            _validate_coefficient_options(m, n, a, spheroid, precision, nothing, max_terms, :orthonormal)
            valid = form===:value ? 0<=y<1 : form===:complement ? 0<y<=1 : y<0
            valid && !isnan(y) || throw(DomainError(target, "invalid concentration target"))
            zero_target = form===:value ? iszero(y) : form===:complement ? isone(y) : y==-Inf
            if zero_target
                iszero(a) || throw(ArgumentError("zero concentration requires a bracket containing c=0"))
                return output(true, a, zero(_SWFloat), 0, :endpoint)
            end
        else
            _validate_coefficient_options(m, n, a, spheroid, precision, nothing, max_terms, :orthonormal)
            _validate_coefficient_options(m, n, b, spheroid, precision, nothing, max_terms, :orthonormal)
            isfinite(y) || throw(DomainError(target, "target must be finite"))
        end
        odds = operator!==:concentration ? zero(_SWFloat) : form===:value ? log(y)-log1p(-y) :
               form===:complement ? log1p(-y)-log(y) : y-log(-expm1(y))
        function evaluate(c)
            if operator===:separation
                eigenpair = _solve_coefficient_eigenpair_at_precision(
                    m, n, c, spheroid, precision, workbits-48, nothing, max_terms, use_jacobian)
                residual = eigenpair.lambda-y
                return (; f = residual, d = use_jacobian ? eigenpair.dlambda : _SWFloat(NaN), residual)
            end
            if iszero(c)
                value = form===:value ? zero(_SWFloat) : form===:complement ? one(_SWFloat) : _SWFloat(-Inf)
                return (; f = _SWFloat(-Inf), d = _SWFloat(Inf), residual = value-y)
            end
            data = _integral_data_at_precision(n, c, precision, workbits-48, max_terms; sensitivity = use_jacobian)
            _with_swprecision(Base.precision(data.concentration)) do
                logarithm = data.concentration>1//2 ? log1p(-data.complement) : log(data.concentration)
                value = form===:value ? data.concentration : form===:complement ? data.complement : logarithm
                (; f = logarithm-log(data.complement)-odds,
                    d = use_jacobian ? data.dlog/data.complement : _SWFloat(NaN), residual = value-y)
            end
        end
        left, right = evaluate(a), evaluate(b)
        residual_tolerance = operator===:concentration ? relative : absolute+relative*max(one(_SWFloat), abs(y))
        accepted(c, result) = iszero(result.f) || (abs(result.f)<=residual_tolerance &&
            (b-a<=absolute+relative*abs(c) ||
             (use_jacobian && isfinite(result.d) && !iszero(result.d) && abs(result.f/result.d)<=absolute+relative*abs(c))))
        accepted(a, left) && return output(true, a, left.residual, 0, :endpoint)
        accepted(b, right) && return output(true, b, right.residual, 0, :endpoint)
        signbit(left.f)!=signbit(right.f) || throw(ArgumentError("bracket endpoints must straddle the requested eigenvalue"))
        c, method = (a+b)/2, :bisection
        for iteration in 1:maxiter
            result = evaluate(c)
            accepted(c, result) && return output(true, c, result.residual, iteration, method)
            if signbit(result.f)==signbit(left.f)
                a, left = c, result
            else
                b, right = c, result
            end
            iteration==maxiter && return output(false, c, result.residual, iteration, :maxiter)
            candidate = use_jacobian ? c-result.f/result.d : b-right.f*(b-a)/(right.f-left.f)
            if iteration%3!=0 && isfinite(candidate) && a<candidate<b && candidate!=c
                c, method = candidate, use_jacobian ? :newton : :secant
            else
                c, method = (a+b)/2, :bisection
            end
            (c==a || c==b) && return output(false, c, (c==a ? left : right).residual, iteration, :unresolved)
        end
    end
end
