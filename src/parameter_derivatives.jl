function _coordinate_parameter_evaluator(
        m, n, c, points, spheroid, precision, target, kind, normalize, normalization)
    R = precision===:quad ? BigFloat : Float64
    center = _input_float(Complex{R}, c)
    prefix = spheroid===:prolate ? :cprolate : :coblate
    state = _complex_mode_state(prefix, m, n, center, precision)
    evaluate = _angular_phase_evaluator(prefix, m, n, precision;
        endpoint_anchor = state.endpoint_anchor, endpoint_coordinate = state.endpoint_point)
    endpoints = findall(z -> _coordinate_is_endpoint(z, target, spheroid), points)
    return parameter -> begin
        local_state = _transport_angular_phase(evaluate, center, state, _input_float(Complex{R}, parameter))
        data = _complex_coordinate_data(
            m, n, parameter, points, spheroid, precision, target,
            kind, normalize, normalization; mode = local_state)
        states = [s[1:2] for s in data.states]
        if !isempty(endpoints)
            # Endpoint stencil values must retain the same transported mode.
            result = _with_swprecision(_coordinate_precision(m, n, parameter, points)) do
                coordinates = real.(points[endpoints])
                if target===:angular && kind==1
                    scale = _coefficient_phase(data.plan, precision; mode = local_state)*(normalize ?
                                                                                          one(_SWFloat) :
                                                                                          sqrt(_ferrers_norm2(
                        m, n, _SWFloat)))
                    _evaluate_coefficient_vector(data.plan, scale .* data.plan.v, coordinates)
                elseif target===:angular
                    _angular_second_kind(m, n, parameter, coordinates, spheroid, precision,
                        false, false, false, false; mode = local_state)
                elseif normalization===:static
                    _static_radial_values(
                        m, n, parameter, coordinates, spheroid, precision,
                        kind; eigenvalue_seed = local_state.lambda)
                else
                    _radial_analytic_values(
                        m, n, parameter, coordinates, spheroid, precision,
                        kind; eigenvalue_seed = local_state.lambda)
                end
            end
            states[endpoints] = collect(zip(result.value, result.derivative))
        end
        return (; states)
    end
end

function _complex_coordinate_jacobian(
        m, n, c, points, spheroid, precision, target, kind, normalize, normalization,
        h, diagnostics, adaptive, rtol, atol, order)
    _validate_jacobian_tolerances(rtol, atol)
    _validate_complex_coordinates(
        m, n, c, points, spheroid, precision, target; kind, normalization)
    kind==2 && normalize &&
        throw(ArgumentError("normalize=true is only defined for angular kind=1"))
    order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
    order==2 && return _wave_curvature(
        m, n, c, points, spheroid, precision, kind, normalize, target,
        h, diagnostics, adaptive, rtol, atol; normalization)
    R = precision===:quad ? BigFloat : Float64
    T = Complex{R}
    if h===nothing
        data = _complex_coordinate_data(
            m, n, c, points, spheroid, precision, target, kind, normalize, normalization)
        value, derivative = ([T(s[j]) for s in data.states] for j in (3, 4))
        mv = (; _coefficient_derivative_metadata(data.plan, value)...,
            method = :differentiated_equation)
        md = (; _coefficient_derivative_metadata(data.plan, derivative)...,
            method = :differentiated_equation)
        iv, id = complex.(-imag.(value), real.(value)),
        complex.(-imag.(derivative), real.(derivative))
        miv, mid = mv, md
    else
        step = _resolve_jacobian_step(c, h, precision)
        evaluate = c isa Complex ?
                   _coordinate_parameter_evaluator(
            m, n, c, points, spheroid, precision, target, kind, normalize, normalization) :
                   parameter -> _complex_coordinate_data(
            m, n, parameter, points, spheroid, precision,
            target, kind, normalize, normalization)
        cache = Dict{Any, Any}()
        sample(parameter) = get!(() -> evaluate(parameter), cache, parameter)
        function difference(direction, index)
            estimate(s) = T[(a[index]-b[index])/(2s)
                            for (a, b) in zip(sample(c+direction*s).states, sample(c-direction*s).states)]
            _finite_difference_diagnostics(
                estimate, step; precision, adaptive, rtol, atol)
        end
        value, mv = difference(1, 1)
        derivative, md = difference(1, 2)
        if c isa Complex
            iv, miv = difference(im, 1)
            id, mid = difference(im, 2)
        end
    end
    if c isa Real
        result = (; dvalue_dc = value, dderivative_dc = derivative)
        return diagnostics ?
               (; result..., metadata_value = mv, metadata_derivative = md) : result
    end
    result = (; dvalue_dcreal = value, dvalue_dcimag = iv,
        dderivative_dcreal = derivative, dderivative_dcimag = id)
    return diagnostics ?
           (; result..., metadata_value_dcreal = mv, metadata_value_dcimag = miv,
        metadata_derivative_dcreal = md, metadata_derivative_dcimag = mid) : result
end

# Analytic sensitivities of the truncated, refined coefficient eigenproblem.
function _parameter_curvature(evaluate, c, precision, h, adaptive, rtol, atol;
        nonnegative = false, nonzero = false)
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? _input_float(T, c) : _input_float(Complex{T}, c)
    epsilon = precision === :quad ? max(T(2)^(-112), eps(T)) : eps(T)
    step = h === nothing ? epsilon^(T(1)/5)*max(one(T), abs(parameter)) :
           _resolve_jacobian_step(c, h, precision)
    # Stay away from the radial pole and the integral operator's left boundary.
    if !iszero(parameter) && (nonnegative || nonzero)
        h === nothing ? (step=epsilon^(T(1)/5)*abs(parameter)) :
        2step >= abs(parameter) && throw(ArgumentError("2h must be smaller than abs(c)"))
    end
    cache = Dict{Any, Any}()
    sample(z) = get!(cache, z) do
        evaluate(z)
    end
    function estimate(s)
        offsets = nonnegative && iszero(parameter) ? (0, 1, 2, 3, 4) : (-2, -1, 1, 2)
        weights = length(offsets)==5 ? (-25, 48, -36, 16, -3) : (1, -8, 8, -1)
        grid = [parameter+k*s for k in offsets]
        length(unique(grid)) == length(grid) ||
            throw(ArgumentError("h is below parameter resolution"))
        center = sample(parameter)
        sum(w .* (sample(z) .- center) for (w, z) in zip(weights, grid)) ./ (12s)
    end
    value, metadata = _finite_difference_diagnostics(
        estimate, step; precision, adaptive, rtol, atol)
    # Identical nonzero samples can hide a small second derivative after rounding.
    # Report this loss of resolution instead of treating a flat stencil as evidence.
    center = sample(parameter)
    unresolved = iszero.(value) .& .!iszero.(center)
    for result in values(cache)
        unresolved = unresolved .& (result .== center)
    end
    if unresolved isa Bool ? unresolved : any(unresolved)
        metadata = (; metadata..., conditioning_flag = :poor,
            suggested_action = precision===:double ? :use_quad : :unresolved)
    end
    return value, (; metadata..., method = :differenced_sensitivity)
end

function _eigen_curvature(m, n, c, spheroid, precision, operator, form,
        h, diagnostics, adaptive, rtol, atol)
    operator === :separation ?
    _validate_wave_arguments(m, n, c, [0], spheroid, precision, :angular) :
    _integral_parameter(m, n, c, spheroid, precision, operator, form)
    T = precision === :quad ? BigFloat : Float64
    # Exact right limits avoid differentiating a singular logarithm at zero.
    if operator !== :separation && iszero(c) && h === nothing
        value = operator === :fourier ?
                complex(n==0 ? -T(2)/9 : n==2 ? -T(8)/45 : zero(T)) :
                form === :log ? T(-Inf) : zero(T)
        metadata = (method = :right_limit, step_used = nothing,
            relative_change_when_halving_step = nothing,
            finite_flag = isfinite(value), conditioning_flag = isfinite(value) ? :good :
                                                               :singular,
            suggested_action = isfinite(value) ? :accept : :singular_limit)
    else
        operator !== :separation && iszero(c) && form === :log &&
            throw(DomainError(c, "log concentration is singular at c=0"))
        evaluate(z) = begin
            result = jacobian_eigen(m, n, z; spheroid, precision, operator, form)
            z isa Real ? result : result.d_dcreal
        end
        value, metadata = _parameter_curvature(
            evaluate, c, precision, h, adaptive, rtol, atol;
            nonnegative = operator!==:separation)
    end
    if c isa Real
        return diagnostics ? (; derivative = value, metadata) : value
    end
    result = (; d2_dcreal2 = value, d2_dcreal_dcimag = im*value, d2_dcimag2 = -value)
    return diagnostics ? (; result..., metadata) : result
end

function _wave_curvature(m, n, c, points, spheroid, precision, kind, normalize, target,
        h, diagnostics, adaptive, rtol, atol; normalization = :standard)
    endpoints = m==1 && kind==1 ?
                findall(
        x -> isreal(x) && (target === :angular ? abs(x)==1 : spheroid === :prolate && x==1),
        points) : Int[]
    evaluate(z) = begin
        result = target === :angular ?
                 jacobian_smn(m, n, z, points; spheroid, precision, kind, normalize) :
                 jacobian_rmn(m, n, z, points; spheroid, precision, kind, normalization)
        values = z isa Real ? hcat(result.dvalue_dc, result.dderivative_dc) :
                 hcat(result.dvalue_dcreal, result.dderivative_dcreal)
        if !isempty(endpoints)
            values[endpoints, 2] = _endpoint_parameter_sensitivity(
                m, n, z, real.(points[endpoints]), spheroid,
                precision, target, normalize, normalization)
        end
        values
    end
    values, metadata = _parameter_curvature(
        evaluate, c, precision, h, adaptive, rtol, atol;
        nonzero = target===:radial && normalization===:standard)
    value, derivative = values[:, 1], values[:, 2]
    for i in endpoints
        derivative[i] = _directed_infinity(derivative[i], c isa Real, precision === :quad ?
                                                                      BigFloat : Float64)
    end
    if any(isinf, derivative)
        metadata = (; metadata..., finite_flag = false, conditioning_flag = :singular,
            suggested_action = :singular_limit)
    end
    result = c isa Real ? (; d2value_dc2 = value, d2derivative_dc2 = derivative) :
             (; d2value_dcreal2 = value,
        d2value_dcreal_dcimag = complex.(-imag.(value), real.(value)),
        d2value_dcimag2 = -value,
        d2derivative_dcreal2 = derivative, d2derivative_dcreal_dcimag = complex.(
            -imag.(derivative), real.(derivative)),
        d2derivative_dcimag2 = -derivative)
    return diagnostics ? (; result..., metadata) : result
end

function _endpoint_parameter_sensitivity(
        m, n, c, points, spheroid, precision, target, normalize, normalization)
    if target === :angular
        result = _coefficient_angular_jacobian(
            m, n, c, points, spheroid, precision, normalize, false; regular_factor = true)
        values = c isa Real ? result.dvalue_dc : result.dvalue_dcreal
        return -points .* values
    end
    value = _with_swprecision(_static_precision(m, n, c)) do
        iszero(c) && return zero(_SWFloat)
        data, plan = _radial_analytic_data(m, n, c, [1], spheroid, precision, 1)
        state = _radial_expansion_data(plan, one(_SWFloat), 1; regular_factor = true).state
        normalization === :static ? plan.c^(-n)*(state[3]-n*state[1]/plan.c) : state[3]
    end
    return fill(value, length(points))
end

function _sensitivity_plan(m, n, c, spheroid, precision)
    plan = _coefficient_plan(m, n, c; spheroid, precision, max_terms = max(512, (n-m)÷2+33))
    plan.converged || error("Parameter sensitivity coefficient expansion did not converge")
    return plan
end

function _coefficient_derivative_metadata(plan, value)
    finite = _all_finite(value)
    good = finite && plan.converged
    # Transpose-normalized complex eigenvectors become self-orthogonal at an
    # exceptional point. A small residual alone must not hide that sensitivity.
    conditioning = !good ? :poor : plan.eigenvalue_condition > 100 ? :warning : :good
    return (method = :coefficients, step_used = nothing,
        relative_change_when_halving_step = nothing,
        coefficient_relative_change = _unwrap_swfloat(plan.relative_change), coefficient_tail = _unwrap_swfloat(plan.tail),
        eigenproblem_residual = _unwrap_swfloat(plan.residual), sensitivity_residual = _unwrap_swfloat(plan.derivative_residual),
        eigenvalue_condition = _unwrap_swfloat(plan.eigenvalue_condition),
        finite_flag = finite, conditioning_flag = conditioning,
        suggested_action = good ? :accept : :use_quad)
end

function _coefficient_eigen_jacobian(m, n, c, spheroid, precision, diagnostics)
    plan = _sensitivity_plan(m, n, c, spheroid, precision)
    T = precision === :quad ? BigFloat : Float64
    value = c isa Real ? T(plan.dlambda) : Complex{T}(plan.dlambda)
    metadata = _coefficient_derivative_metadata(plan, value)
    if c isa Real
        return diagnostics ? (; derivative = value, metadata) : value
    end
    result = (d_dcreal = value, d_dcimag = im*value)
    return diagnostics ?
           (; result..., metadata_dcreal = metadata, metadata_dcimag = metadata) : result
end

function _coefficient_angular_jacobian(
        m, n, c, eta, spheroid, precision, normalize, diagnostics; regular_factor = false)
    _validate_wave_arguments(m, n, c, eta, spheroid, precision, :angular)
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? _input_float(T, c) : _input_float(Complex{T}, c)
    coordinates = _SWFloat.(eta)
    guarded = _use_angular_expansion(n, parameter, spheroid) ||
              (parameter isa Complex && iszero(real(parameter)) && abs(parameter)>n+1)
    bits = guarded ? _angular_precision(parameter) : Base.precision(_SWFloat)
    plan, result = _with_swprecision(bits) do
        plan = guarded ? _guarded_angular_plan(m, n, parameter, spheroid, precision) :
               _sensitivity_plan(m, n, parameter, spheroid, precision)
        phase = _coefficient_phase(plan, precision)
        scale = phase*(normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(m, n, _SWFloat)))
        plan,
        _evaluate_coefficient_vector(plan, plan.dv .* scale, coordinates; regular_factor)
    end
    convert_values(v) = c isa Real ? T.(v) : complex.(T.(real.(v)), T.(imag.(v)))
    value, derivative = convert_values(result.value), convert_values(result.derivative)
    mv = _coefficient_derivative_metadata(plan, value)
    md = _coefficient_derivative_metadata(plan, derivative)
    if c isa Real
        result = (dvalue_dc = value, dderivative_dc = derivative)
        return diagnostics ?
               (; result..., metadata_value = mv, metadata_derivative = md) : result
    end
    result = (dvalue_dcreal = value, dvalue_dcimag = complex.(-imag.(value), real.(value)),
        dderivative_dcreal = derivative, dderivative_dcimag = complex.(-imag.(derivative), real.(derivative)))
    return diagnostics ?
           (; result..., metadata_value_dcreal = mv, metadata_value_dcimag = mv,
        metadata_derivative_dcreal = md, metadata_derivative_dcimag = md) : result
end

function _qs_angular_jacobian(m, n, c, eta, spheroid, precision, normalize, diagnostics)
    normalize &&
        throw(ArgumentError("normalize=true is only defined for angular kind=1; Qs uses the DLMF second-kind normalization"))
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? _input_float(T, c) : _input_float(Complex{T}, c)
    coordinates = _SWFloat.(eta)
    result = _with_swprecision(_qs_precision(parameter)) do
        _qs_values(m, n, parameter, coordinates, spheroid, precision; sensitivity = true)
    end
    convert_values(v) = c isa Real ? T.(v) : complex.(T.(real.(v)), T.(imag.(v)))
    value, derivative = convert_values(result.dvalue_dc),
    convert_values(result.dderivative_dc)
    metadata(v) = (; _coefficient_derivative_metadata(result.plan, v)...,
        method = :differentiated_equation)
    mv, md = metadata(value), metadata(derivative)
    if c isa Real
        output = (dvalue_dc = value, dderivative_dc = derivative)
        return diagnostics ?
               (; output..., metadata_value = mv, metadata_derivative = md) : output
    end
    output = (dvalue_dcreal = value, dvalue_dcimag = im .* value,
        dderivative_dcreal = derivative, dderivative_dcimag = im .* derivative)
    return diagnostics ?
           (; output..., metadata_value_dcreal = mv, metadata_value_dcimag = mv,
        metadata_derivative_dcreal = md, metadata_derivative_dcimag = md) : output
end
