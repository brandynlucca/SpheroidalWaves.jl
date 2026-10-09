# Shared validation keeps diagnostics and evaluations in the same domain.
function _validate_wave_arguments(m, n, c, points, spheroid, precision, target; kind=1,normalization=:standard)
    _validate_precision(precision)
    spheroid in (:prolate, :oblate) || error("spheroid must be :prolate or :oblate")
    0 <= m <= n || error("require 0 <= m <= n")
    isfinite(c) || error("c must be finite")
    isempty(points) && error("evaluation points must not be empty")
    all(isfinite, points) || error("evaluation points must be finite")
    if target === :angular
        normalization === :standard || throw(ArgumentError("radial normalization does not apply to angular functions"))
        kind in (1,2) || throw(ArgumentError("angular kind must be 1 or 2"))
        all(x -> abs(x) <= 1, points) || error("eta must lie in [-1, 1]")
    else
        normalization in (:standard,:static) || throw(ArgumentError("normalization must be :standard or :static"))
        normalization === :standard && _validate_radial_parameter(c)
        normalization === :static && !(kind in (1,2)) && throw(ArgumentError("static normalization requires kind=1 or 2"))
        kind in 1:4 || error("kind must be in 1:4")
        boundary = spheroid === :prolate ? 1 : 0
        all(x -> x >= boundary, points) || error("radial coordinates must be >= $boundary for $spheroid")
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
function _coordinate_accuracy_record(original,refined,plan,point,precision,endpoint)
    R = precision===:quad ? BigFloat : Float64
    relative(a,b) = iszero(b) ? (iszero(a) ? zero(_SWFloat) : _SWFloat(Inf)) : abs(a)/abs(b)
    changes = abs.((original[1]-refined[1],original[2]-refined[2]))
    relative_changes = (relative(changes[1],refined[1]),relative(changes[2],refined[2]))
    # |z*f_z/f| and |c*f_c/f| measure sensitivity to relative input perturbations.
    coordinate_condition = iszero(refined[1]) ? _SWFloat(Inf) : relative(Complex{_SWFloat}(point)*refined[2],refined[1])
    parameter_condition = iszero(refined[1]) ? _SWFloat(Inf) : relative(plan.c*refined[3],refined[1])
    eigenvalue_condition = plan.eigenvalue_condition
    condition = max(one(_SWFloat),coordinate_condition,parameter_condition,eigenvalue_condition)
    # Quad uses binary128 spacing, or the caller's BigFloat spacing if coarser.
    epsilon = precision===:quad ? max(_SWFloat(2)^(-112),_SWFloat(eps(BigFloat))) : _SWFloat(eps(Float64))
    # Include public output rounding, particularly overflow and underflow.
    rounded = Complex{_SWFloat}.(Complex{R}.(original[1:2]))
    finite = all(isfinite,rounded) && all(isfinite,refined)
    # Require relative value and slope changes <= 32*epsilon, a heuristic threshold.
    converged = !endpoint && finite && all(r -> r<=32epsilon,relative_changes)
    coefficient_error = max(plan.relative_change,plan.tail,plan.residual)*eigenvalue_condition
    # Take the largest error indicator: refinement change with a heuristic 4x
    # margin, output rounding, coefficient errors, or rounding amplified by conditioning.
    uncertainty = max(4relative_changes[1],relative(rounded[1]-refined[1],refined[1]),
                      coefficient_error,epsilon*condition)
    available = converged && !iszero(refined[1])
    digits = !available ? -1 : !isfinite(uncertainty) ? 0 :
        clamp(floor(Int,-log10(uncertainty)),0,precision===:quad ? 33 : 16)
    # Warn when conditioning can cost more than two decimal digits.
    flag = endpoint ? :unavailable : !finite ? :singular : !converged ? :poor : condition>100 ? :warning : :good
    diagnostics = (method=:continuation_refinement,converged,finite_flag=finite,conditioning_flag=flag,
        absolute_change_value=R(changes[1]),absolute_change_derivative=R(changes[2]),
        relative_change_value=R(relative_changes[1]),relative_change_derivative=R(relative_changes[2]),
        coordinate_condition=R(coordinate_condition),parameter_condition=R(parameter_condition),
        eigenvalue_condition=R(eigenvalue_condition),coefficient_tail=R(plan.tail))
    return (;digits,diagnostics)
end

function _complex_coordinate_accuracy(m,n,c,points,spheroid,precision,target,kind,normalize,normalization,diagnostics)
    _validate_complex_coordinates(m,n,c,points,spheroid,precision,target;kind,normalization)
    target===:angular && kind==2 && normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1"))
    angular_normalization = target===:angular && normalize
    # Match wave evaluation, then independently tighten arithmetic and steps.
    bits = _coordinate_precision(m,n,c,points)
    original = _with_swprecision(bits) do
        _complex_coordinate_data(m,n,c,points,spheroid,precision,target,kind,angular_normalization,normalization)
    end
    refined = _with_swprecision(bits+64) do
        _complex_coordinate_data(m,n,c,points,spheroid,precision,target,kind,angular_normalization,normalization;step_fraction=0.5)
    end
    records = _with_swprecision(bits+64) do
        [_coordinate_accuracy_record(a,b,refined.plan,z,precision,_coordinate_is_endpoint(z,target,spheroid))
         for (a,b,z) in zip(original.states,refined.states,points)]
    end
    digits = [record.digits for record in records]
    return diagnostics ? (;digits,diagnostics=[record.diagnostics for record in records]) : digits
end
