# Follow the eigenfunction, rather than the backend's sorted degree label.
# Native complex normalization fixes the square only, so its remaining sign
# must also be transported. Candidate labels always have the original parity.
struct _AngularModeEvaluator{F,G}
    sample::F
    eigenvalue::G
end
(evaluate::_AngularModeEvaluator)(args...) = evaluate.sample(args...)

function _angular_phase_evaluator(prefix, m, n, precision;endpoint_anchor=nothing,endpoint_coordinate=nothing)
    T = precision === :quad ? BigFloat : Float64
    points = T[0, 1//4, 1//2, 3//4]
    cache = Dict{Any,Any}()
    lambda_cache = Dict{Any,Any}()
    eigen = (c,degree) -> get!(lambda_cache,(c,degree)) do
        _call_complex_eigenvalue(prefix,m,degree,c;precision)
    end
    use_endpoint = Ref{Union{Nothing,Bool}}(endpoint_anchor)
    endpoint_point = Ref(endpoint_coordinate === nothing ? zero(T) : T(endpoint_coordinate))
    sample = (c, degree=n) -> get!(cache,(c,degree)) do
        r = _call_complex_smn_raw(prefix, m, degree, c, points; precision, normalize=true)
        lambda = eigen(c,degree)
        anchor = iseven(n-m) ? r.value[1] : r.derivative[1]
        profile = [r.value; r.derivative ./ (n+1)]
        if use_endpoint[] === nothing
            use_endpoint[] = !(isfinite(anchor) && !iszero(anchor))
        end
        if use_endpoint[]
            if iszero(endpoint_point[])
                distance = T(2)^(-ceil(Int,log2(32*(1+T(n)*(n+1)+abs(c)^2))))
                endpoint_point[] = 1-distance
            end
            endpoint = _call_complex_smn_raw(prefix,m,degree,c,[endpoint_point[]];precision,normalize=true)
            anchor = only(endpoint.value)
            append!(profile,[anchor,only(endpoint.derivative)/(n+1)])
        end
        isfinite(lambda) && isfinite(anchor) && !iszero(anchor) && all(isfinite, profile) ||
            error("Cannot continue angular phase: nonfinite samples or a zero phase anchor at c=$c")
        (; m, n=Int(degree), lambda, anchor, profile, sign=1,endpoint_anchor=use_endpoint[],endpoint_point=endpoint_point[])
    end
    return _AngularModeEvaluator(sample,eigen)
end

function _signed_angular_state(state, sign)
    return (; state..., anchor=sign*state.anchor, profile=sign.*state.profile, sign=sign*state.sign)
end

function _align_angular_state(reference, candidate)
    overlap = real(conj(reference.anchor / abs(reference.anchor)) *
                   (candidate.anchor / abs(candidate.anchor)))
    return _signed_angular_state(candidate, overlap < 0 ? -1 : 1)
end

_mode_eigenvalue_distance(a,b) = abs(a-b)/(1+max(abs(a),abs(b)))
_candidate_distance_bound(evaluate,reference,c,degree) = 0
_candidate_distance_bound(evaluate::_AngularModeEvaluator,reference,c,degree) =
    _mode_eigenvalue_distance(reference.lambda,evaluate.eigenvalue(c,degree))

function _angular_state_distance(a, b)
    anchor_error = abs(a.anchor-b.anchor) / max(abs(a.anchor), abs(b.anchor))
    scale = max(maximum(abs, a.profile), maximum(abs, b.profile))
    profile_error = maximum(abs, a.profile ./ scale .- b.profile ./ scale)
    eigenvalue_error = _mode_eigenvalue_distance(a.lambda,b.lambda)
    return max(anchor_error, profile_error, eigenvalue_error)
end

function _nearest_angular_state(evaluate, reference, c, branch_window)
    best = _align_angular_state(reference, evaluate(c, reference.n))
    score = _angular_state_distance(reference, best)
    runner_up = Inf
    lowest = reference.m + mod(reference.n-reference.m,2)
    for degree in max(lowest,reference.n-2branch_window):2:reference.n+2branch_window
        degree == reference.n && continue
        # Eigenvalue distance is a lower bound for the full matching score.
        # Exclude clearly separated modes before computing their angular data.
        bound = _candidate_distance_bound(evaluate,reference,c,degree)
        if bound > 2score
            runner_up = min(runner_up,bound)
            continue
        end
        candidate = _align_angular_state(reference, evaluate(c, degree))
        distance = _angular_state_distance(reference, candidate)
        if distance < score
            runner_up = score
            best,score = candidate,distance
        else
            runner_up = min(runner_up,distance)
        end
    end
    # Refine when two modes are similarly plausible. A small absolute change
    # alone is insufficient: nearby eigenvalues can exchange native labels.
    return (;best...,resolved=score < runner_up/2)
end

function _transport_angular_phase(evaluate, c0, state0, c1; branch_window=1)
    isfinite(c0) && isfinite(c1) || error("continuation path must be finite")
    c0 == c1 && return state0
    evaluations = Ref(0)
    function advance(a, initial, b, depth)
        midpoint = (a+b)/2
        if midpoint == a || midpoint == b
            last = _nearest_angular_state(evaluate, initial, b, branch_window)
            last.resolved && _angular_state_distance(initial, last) <= 0.2 && return last
        end
        if depth > 24 || midpoint == a || midpoint == b || evaluations[] >= 4096
            error("Complex mode continuation could not resolve the segment $a to $b. " *
                  "The path may approach a branch point or leave the available native mode window.")
        end
        # Limit parameter steps even when endpoint values happen to be close.
        if abs(b-a) > 0.5
            middle = advance(a, initial, midpoint, depth+1)
            return advance(midpoint, middle, b, depth+1)
        end
        middle = _nearest_angular_state(evaluate, initial, midpoint, branch_window)
        last = _nearest_angular_state(evaluate, middle, b, branch_window)
        evaluations[] += 2
        if middle.resolved && last.resolved && max(_angular_state_distance(initial, middle),
               _angular_state_distance(middle, last)) <= 0.2
            return last
        end
        middle = advance(a, initial, midpoint, depth+1)
        return advance(midpoint, middle, b, depth+1)
    end
    return advance(c0, state0, c1, 0)
end

const _complex_mode_cache_key = gensym(:complex_mode_cache)

function _complex_mode_state(prefix,m,n,c,precision)
    T = precision === :quad ? BigFloat : Float64
    parameter = Complex{T}(c)
    isfinite(parameter) || error("c is not finite at the requested precision")
    # Reuse the canonical path across coordinates, radial kinds, and diagnostics.
    # Task-local storage avoids shared mutable continuation state. Backend and
    # BigFloat context are part of the key; failed paths are never cached.
    cache = get!(task_local_storage(),_complex_mode_cache_key) do
        Dict{Any,Any}()
    end
    key = (prefix,m,n,parameter,precision,Base.precision(BigFloat),rounding(BigFloat),backend_library(;precision))
    haskey(cache,key) && return cache[key]
    evaluate = _angular_phase_evaluator(prefix,m,n,precision)
    state = _initial_angular_phase(evaluate,m,n,parameter)
    length(cache) >= 64 && empty!(cache)
    cache[key] = state
    return state
end

function _angular_mode_factor(m,n,state,normalize,T)
    (normalize || n == state.n) && return T(state.sign)
    # A backend relabeling must not change the requested integral norm.
    return T(state.sign*sqrt(_ferrers_norm2(m,n,BigFloat)/_ferrers_norm2(m,state.n,BigFloat)))
end

_radial_mode_factor(n,state) = isodd((state.n-n)÷2) ? -1 : 1

function _scale_mode_result!(result,factor)
    result.value .*= factor
    result.derivative .*= factor
    return result
end

function _initial_angular_phase(evaluate, m, n, c;branch_window=1)
    start = complex(real(c), zero(real(c)))
    state = evaluate(start)
    # On the real axis this is the Legendre/DLMF sign, before (-1)^m.
    negative = state.endpoint_anchor ? false : isodd((n-m)÷2)
    sign = signbit(real(state.anchor)) == negative ? 1 : -1
    state = _signed_angular_state(state, sign)
    return _transport_angular_phase(evaluate, start, state, c;branch_window)
end

function _call_complex_smn(prefix, m, n, c, eta; precision=:double, normalize=false)
    iszero(c) && return _spherical_smn_complex(m, n, eta; precision, normalize)
    T = precision === :quad ? BigFloat : Float64
    state = _require_angular_anchor(_complex_mode_state(prefix,m,n,c,precision))
    result = _call_complex_smn_raw(prefix, m, state.n, Complex{T}(c), eta; precision, normalize)
    return _scale_mode_result!(result,_angular_mode_factor(m,n,state,normalize,T))
end

function _require_angular_anchor(state)
    # A boundary sample can still identify eigenvalues/radial modes when the
    # angular backend loses its interior values. Do not return those spurious
    # zeros as successful angular evaluations or accuracy estimates.
    state.endpoint_anchor && error("Angular backend lost its nonzero origin anchor; interior angular values cannot be trusted")
    return state
end

function _call_complex_rmn(prefix,m,n,c,x;precision=:double,kind=1)
    state = _complex_mode_state(prefix,m,n,c,precision)
    result = _call_complex_rmn_raw(prefix,m,state.n,c,x;precision,kind)
    return _scale_mode_result!(result,_radial_mode_factor(n,state))
end

# Each stencil is continued from its center, so crossing a cut of the scalar
# reference path cannot turn a small finite difference into a sign jump.
function _angular_jacobian_evaluator(m, n, c, eta; spheroid, precision, normalize,kind=1)
    _validate_wave_arguments(m,n,c,eta,spheroid,precision,:angular;kind)
    T = precision === :quad ? BigFloat : Float64
    center = Complex{T}(c)
    prefix = spheroid === :prolate ? :cprolate : :coblate
    state = _require_angular_anchor(_complex_mode_state(prefix,m,n,center,precision))
    phase_evaluate = _angular_phase_evaluator(prefix, m, n, precision)
    cache = Dict{Complex{T}, Any}()
    return parameter -> get!(cache, Complex{T}(parameter)) do
        local_state = _transport_angular_phase(phase_evaluate, center, state, Complex{T}(parameter))
        if kind == 2
            return _angular_second_kind(m,n,Complex{T}(parameter),eta,spheroid,precision,
                                        normalize,false,false,false;mode=local_state)
        end
        result = _call_complex_smn_raw(prefix, m, local_state.n, Complex{T}(parameter), eta; precision, normalize)
        factor = _angular_mode_factor(m,n,local_state,normalize,T) * (isodd(m) ? -1 : 1)
        return _scale_mode_result!(result,factor)
    end
end

# Eigenvalue and radial finite differences follow the same local branch as the
# angular stencil, including when the scalar reference path has a cut nearby.
function _complex_local_evaluator(m,n,c,spheroid,precision;points=nothing,kind=1)
    _validate_wave_arguments(m,n,c,points === nothing ? [0] : points,spheroid,precision,
                             points === nothing ? :angular : :radial;kind)
    T = precision === :quad ? BigFloat : Float64
    center = Complex{T}(c)
    prefix = spheroid === :prolate ? :cprolate : :coblate
    state = _complex_mode_state(prefix,m,n,center,precision)
    evaluate = _angular_phase_evaluator(prefix,m,n,precision;
        endpoint_anchor=state.endpoint_anchor,endpoint_coordinate=state.endpoint_point)
    cache = Dict{Complex{T},Any}()
    return parameter -> get!(cache,Complex{T}(parameter)) do
        points !== nothing && _validate_radial_parameter(parameter)
        local_state = _transport_angular_phase(evaluate,center,state,Complex{T}(parameter))
        points === nothing && return local_state.lambda
        result = _call_complex_rmn_raw(prefix,m,local_state.n,parameter,points;precision,kind)
        return _scale_mode_result!(result,_radial_mode_factor(n,local_state))
    end
end
