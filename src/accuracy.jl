# Shared validation keeps diagnostics and evaluations in the same domain.
function _validate_wave_arguments(
        m, n, c, points, spheroid, precision, target; kind = 1, normalization = :standard)
    _validate_precision(precision)
    spheroid in (:prolate, :oblate) || error("spheroid must be :prolate or :oblate")
    0 <= m <= n || error("require 0 <= m <= n")
    isfinite(c) || error("c must be finite")
    isempty(points) && error("evaluation points must not be empty")
    all(isfinite, points) || error("evaluation points must be finite")
    if target === :angular
        normalization === :standard ||
            throw(ArgumentError("radial normalization does not apply to angular functions"))
        kind in (1, 2) || throw(ArgumentError("angular kind must be 1 or 2"))
        all(x -> abs(x) <= 1, points) || error("eta must lie in [-1, 1]")
    else
        normalization in (:standard, :static) ||
            throw(ArgumentError("normalization must be :standard or :static"))
        normalization === :standard && _validate_radial_parameter(c)
        normalization === :static && !(kind in (1, 2)) &&
            throw(ArgumentError("static normalization requires kind=1 or 2"))
        kind in 1:4 || error("kind must be in 1:4")
        boundary = spheroid === :prolate ? 1 : 0
        all(x -> x >= boundary, points) ||
            error("radial coordinates must be >= $boundary for $spheroid")
    end
    return nothing
end

# -1 means unavailable, never a claim of high precision. Reject nonsensical
# native estimates instead of converting them into apparently reliable digits.
function _reported_accuracy(estimates, values, precision)
    maxdigits = precision === :quad ? 33 : 16
    return [isfinite(value) && 0 <= estimate <= maxdigits ? Int(estimate) : -1
            for (estimate, value) in zip(estimates, values)]
end

function _accuracy_epsilon(precision)
    precision === :quad ?
    max(_SWFloat(2)^(-112), _SWFloat(eps(BigFloat))) : _SWFloat(eps(Float64))
end

function _accuracy_unavailable(status)
    return (; digits = -1, absolute_error = BigFloat(NaN),
        relative_error = BigFloat(NaN), status)
end

function _accuracy_singularities(m, point, target, spheroid, kind)
    boundary = isreal(point) && (target === :angular ? abs(real(point)) == 1 :
                spheroid === :prolate && real(point) == 1)
    return (boundary && kind >= 2, boundary && (kind >= 2 || m == 1))
end

# Estimates are deliberately heuristic. Keep absolute errors in BigFloat so an
# underflowed double result cannot turn a small, nonzero error estimate into zero.
function _accuracy_component(actual, coarse, refined, scale, precision;
        coefficient_error = zero(_SWFloat), supported = true, singular = false)
    singular && return _accuracy_unavailable(:singularity)
    isfinite(actual) || return _accuracy_unavailable(:nonfinite)
    supported || return _accuracy_unavailable(:unavailable)
    all(isfinite, (coarse, refined, scale, coefficient_error)) ||
        return _accuracy_unavailable(:refinement_failed)
    epsilon = _accuracy_epsilon(precision)
    change = abs(coarse-refined)
    change <= 32epsilon*scale && coefficient_error <= 32epsilon ||
        return _accuracy_unavailable(:refinement_failed)
    uncertainty = max(4abs(actual-refined), 4change,
        coefficient_error*scale, epsilon*scale)
    relative = iszero(refined) ? _SWFloat(NaN) : uncertainty/abs(refined)
    digits = isnan(relative) ? -1 :
             iszero(relative) ? (precision === :quad ? 33 : 16) :
             !isfinite(relative) ? 0 :
             clamp(floor(Int, -log10(relative)), 0, precision === :quad ? 33 : 16)
    return (; digits, absolute_error = BigFloat(uncertainty),
        relative_error = BigFloat(relative), status = :estimated)
end

function _accuracy_scales(state, plan, point, target, spheroid)
    value, derivative = state[1:2]
    z = Complex{_SWFloat}(point)
    a, b, d = _wave_equation(plan.m, plan.c, z, plan.lambda, spheroid, target)
    second = iszero(a) ? zero(value) : -(b*derivative+d*value)/a
    # Absolute sensitivities remain meaningful at a zero or stationary point.
    parameter_value = length(state) >= 4 ? abs(plan.c*state[3]) : zero(_SWFloat)
    parameter_derivative = length(state) >= 4 ? abs(plan.c*state[4]) : zero(_SWFloat)
    return (max(abs(value), abs(z*derivative), parameter_value),
        max(abs(derivative), abs(z*second), parameter_derivative))
end

function _accuracy_refinement(calculate)
    try
        return calculate()
    catch exception
        # Original evaluations and invalid arguments still throw. A failed
        # additional numerical comparison must not masquerade as an estimate.
        exception isa Union{ErrorException, DomainError, SingularException} || rethrow()
        return nothing
    end
end

function _real_accuracy_reference(m, n, c, points, spheroid, precision, target,
        kind, normalize; min_terms = 0, step_fraction = 1)
    tolerance = sqrt(eps(_SWFloat))
    if target === :angular
        plan = _coefficient_plan(m, n, c; spheroid, precision, rtol = tolerance,
            min_terms, max_terms = max(512, min_terms+16, (n-m)÷2+2ceil(Int, abs(c))+32))
        plan.converged || error("Accuracy reference coefficients did not converge")
        scale = _use_small_parameter_expansion(c) ?
                _SWFloat(real(plan.v[(n - m) ÷ 2 + 1]) < 0 ? -1 : 1) :
                _coefficient_phase(plan, precision)
        normalize || (scale *= sqrt(_ferrers_norm2(m, n, _SWFloat)))
        wave = _evaluate_coefficient_vector(plan, scale .* plan.v, points)
        tangent = _evaluate_coefficient_vector(plan, scale .* plan.dv, points)
        states = collect(zip(wave.value, wave.derivative, tangent.value, tangent.derivative))
        return (; states, plan)
    end
    required_terms = min_terms
    coefficients = (
        m, n, c; min_terms = 0, max_terms = 512, kwargs...) -> _coefficient_plan(
        m, n, c; kwargs..., min_terms = max(min_terms, required_terms),
        max_terms = max(max_terms, required_terms+16))
    expansion = (plan, points, kind) -> _radial_expansion_batch(
        plan, points, kind; step_fraction, rtol = tolerance)
    data, plan = _radial_analytic_data(m, n, c, points, spheroid, precision, kind;
        rtol = tolerance, coefficient_plan = coefficients, expansion_batch = expansion)
    return (; data.states, plan)
end

function _real_accuracy_diagnostics(m, n, c, points, spheroid, precision, target,
        kind, normalize, normalization; reference = _real_accuracy_reference)
    wave = target === :angular ?
           smn(m, n, c, points; spheroid, precision, kind, normalize) :
           rmn(m, n, c, points; spheroid, precision, kind, normalization)
    supported = (target === :angular && kind == 1) ||
                (target === :radial && normalization === :standard)
    interior = findall(eachindex(points)) do i
        supported && !_coordinate_is_endpoint(points[i], target, spheroid) &&
            isfinite(wave.value[i]) && isfinite(wave.derivative[i])
    end
    bits = max(384, _angular_precision(c), _static_precision(m, n, c))
    comparison = isempty(interior) ? nothing :
                 _accuracy_refinement() do
        coarse = _with_swprecision(bits) do
            reference(
                m, n, c, points[interior], spheroid, precision, target, kind, normalize)
        end
        refined = _with_swprecision(bits+64) do
            reference(
                m, n, c, points[interior], spheroid, precision, target, kind, normalize;
                min_terms = coarse.plan.terms+16, step_fraction = 0.5)
        end
        (; coarse, refined)
    end
    return _with_swprecision(bits+128) do
        map(eachindex(points)) do i
            singular = _accuracy_singularities(m, points[i], target, spheroid, kind)
            index = findfirst(==(i), interior)
            method = !supported || index === nothing ? :unavailable :
                     target === :angular ? :coefficient_refinement : :radial_refinement
            estimates = if comparison === nothing || index === nothing
                map((wave.value[i], wave.derivative[i]), singular) do value, is_singular
                    status = is_singular ? :singularity :
                             !isfinite(value) ? :nonfinite :
                             index === nothing ? :unavailable : :refinement_failed
                    _accuracy_unavailable(status)
                end
            else
                coarse, refined = comparison.coarse.states[index],
                comparison.refined.states[index]
                plan = comparison.refined.plan
                scales = _accuracy_scales(refined, plan, points[i], target, spheroid)
                error = max(plan.relative_change, plan.tail, plan.residual)*plan.eigenvalue_condition
                map((wave.value[i], wave.derivative[i]), coarse[1:2], refined[1:2], scales) do a,
                b, r, scale
                    _accuracy_component(Complex{_SWFloat}(a), b, r, scale, precision;
                        coefficient_error = error)
                end
            end
            terms = comparison === nothing || index === nothing ? nothing :
                    (comparison.coarse.plan.terms, comparison.refined.plan.terms)
            working_precision = terms === nothing ? nothing :
                                (Base.precision(real(comparison.coarse.plan.c)),
                Base.precision(real(comparison.refined.plan.c)))
            (; method, value = estimates[1],
                derivative = estimates[2], working_precision, terms)
        end
    end
end

"""
    _coordinate_accuracy_record(original, refined, plan, point, precision, endpoint)

Estimate value digits and conditioning at one complex coordinate `point`.
Both states contain `(value, dvalue_dz, dvalue_dc, d2value_dzdc)`. `original`
matches the requested evaluation. `refined` uses 64 extra arithmetic bits and
halved continuation steps. `plan` supplies the refined coefficient diagnostics
and eigenvalue condition number.

Return `(; digits, diagnostics)`. Endpoint limits, zero values, nonfinite results,
and unresolved refinement receive `digits = -1`. Digit counts are heuristic
estimates of value accuracy, not certified bounds or derivative accuracy estimates.
"""
function _coordinate_accuracy_record(original, refined, plan, point, precision, endpoint;
        target = :angular, spheroid = plan.spheroid, kind = 1, refinement_failed = false)
    R = precision===:quad ? BigFloat : Float64
    relative(a, b) = iszero(b) ? (iszero(a) ? zero(_SWFloat) : _SWFloat(Inf)) :
                     abs(a)/abs(b)
    changes = abs.((original[1]-refined[1], original[2]-refined[2]))
    relative_changes = (relative(changes[1], refined[1]), relative(changes[2], refined[2]))
    # |z*f_z/f| and |c*f_c/f| measure sensitivity to relative input perturbations.
    coordinate_condition = iszero(refined[1]) ? _SWFloat(Inf) :
                           relative(Complex{_SWFloat}(point)*refined[2], refined[1])
    parameter_condition = iszero(refined[1]) ? _SWFloat(Inf) :
                          relative(plan.c*refined[3], refined[1])
    eigenvalue_condition = plan.eigenvalue_condition
    condition = max(one(_SWFloat), coordinate_condition, parameter_condition, eigenvalue_condition)
    # Quad uses binary128 spacing, or the caller's BigFloat spacing if coarser.
    epsilon = _accuracy_epsilon(precision)
    # Include public output rounding, particularly overflow and underflow.
    rounded = Complex{_SWFloat}.(Complex{R}.(original[1:2]))
    finite = all(isfinite, rounded) && all(isfinite, refined)
    # Require relative value and slope changes <= 32*epsilon, a heuristic threshold.
    converged = !endpoint && !refinement_failed && finite &&
                all(r -> r<=32epsilon, relative_changes)
    coefficient_error = max(plan.relative_change, plan.tail, plan.residual)*eigenvalue_condition
    # Take the largest error indicator: refinement change with a heuristic 4x
    # margin, output rounding, coefficient errors, or rounding amplified by conditioning.
    uncertainty = max(4relative_changes[1], relative(rounded[1]-refined[1], refined[1]),
        coefficient_error, epsilon*condition)
    available = converged && !iszero(refined[1])
    digits = !available ? -1 :
             !isfinite(uncertainty) ? 0 :
             clamp(floor(Int, -log10(uncertainty)), 0, precision===:quad ? 33 : 16)
    # Warn when conditioning can cost more than two decimal digits.
    flag = endpoint ? :unavailable :
           !finite ? :nonfinite : !converged ? :poor : condition>100 ? :warning : :good
    singular = _accuracy_singularities(plan.m, point, target, spheroid, kind)
    scales = endpoint ? (one(_SWFloat), one(_SWFloat)) :
             _accuracy_scales(refined, plan, point, target, spheroid)
    estimates = map(1:2) do j
        if refinement_failed && !endpoint && !singular[j] && isfinite(rounded[j])
            return _accuracy_unavailable(:refinement_failed)
        end
        _accuracy_component(rounded[j], original[j], refined[j], scales[j], precision;
            coefficient_error, supported = !endpoint, singular = singular[j])
    end
    diagnostics = (method = :continuation_refinement, converged,
        finite_flag = finite, conditioning_flag = flag,
        absolute_change_value = R(changes[1]), absolute_change_derivative = R(changes[2]),
        relative_change_value = R(relative_changes[1]), relative_change_derivative = R(relative_changes[2]),
        coordinate_condition = R(coordinate_condition), parameter_condition = R(parameter_condition),
        eigenvalue_condition = R(eigenvalue_condition), coefficient_tail = R(plan.tail),
        value = estimates[1], derivative = estimates[2])
    return (; digits, diagnostics)
end

function _complex_coordinate_accuracy(m, n, c, points, spheroid, precision, target,
        kind, normalize, normalization, diagnostics)
    _validate_complex_coordinates(
        m, n, c, points, spheroid, precision, target; kind, normalization)
    target===:angular && kind==2 && normalize &&
        throw(ArgumentError("normalize=true is only defined for angular kind=1"))
    angular_normalization = target===:angular && normalize
    # Match wave evaluation, then independently tighten arithmetic and steps.
    bits = _coordinate_precision(m, n, c, points)
    original = _with_swprecision(bits) do
        _complex_coordinate_data(m, n, c, points, spheroid, precision, target,
            kind, angular_normalization, normalization)
    end
    refine() = _with_swprecision(bits+64) do
        _complex_coordinate_data(m, n, c, points, spheroid, precision, target, kind,
            angular_normalization, normalization; step_fraction = 0.5)
    end
    refined = diagnostics ? _accuracy_refinement(refine) : refine()
    refinement_failed = refined === nothing
    refinement_failed && (refined = original)
    records = _with_swprecision(bits+64) do
        [_coordinate_accuracy_record(a, b, refined.plan, z, precision,
             _coordinate_is_endpoint(z, target, spheroid); target, spheroid, kind, refinement_failed)
         for (a, b, z) in zip(original.states, refined.states, points)]
    end
    digits = [record.digits for record in records]
    working_precision = (
        Base.precision(real(original.plan.c)), Base.precision(real(refined.plan.c)))
    terms = (original.plan.terms, refined.plan.terms)
    return diagnostics ?
           (; digits,
        diagnostics = [(; record.diagnostics..., working_precision, terms)
                       for record in records]) : digits
end
