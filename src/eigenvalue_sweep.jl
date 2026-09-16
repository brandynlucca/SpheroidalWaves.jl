"""
Track `lambda_mn(c)` over a real monotone `c` grid using continuation with branch lock.

This API is intended for stable parameter sweeps where accidental mode-hopping can
occur in challenging regimes. At each new `c` point, a predictor value is formed
from previous continuation states and the selected branch is corrected by choosing
the candidate closest to that predictor from a local degree window.

Keyword arguments:
- `spheroid`: `:prolate` or `:oblate`.
- `precision`: `:double` or `:quad`.
- `branch_lock`: if `true`, enables local branch selection around the previous
  degree index using predictor proximity.
- `branch_window`: nonnegative integer defining the half-width of the local degree
  search window around the previous selected degree.
- `use_jacobian_predictor`: if `true`, uses `jacobian_eigen` for first-step
  predictor initialization and falls back to a zeroth-order predictor when Jacobian
  diagnostics are not acceptable.

Returns a named tuple with fields:
- `c::Vector{T}`
- `lambda::Vector{T}`
- `selected_n::Vector{Int}`
- `switched_branch::Vector{Bool}`

`T` is `Float64` for double precision and `BigFloat` for quad precision.
Grid coordinates, predictors, and candidate comparisons retain this precision.
"""
function eigenvalue_sweep(m::Integer,
                          n::Integer,
                          c_grid::AbstractVector{<:Real};
                          spheroid::Symbol=:prolate,
                          precision::Symbol=:double,
                          branch_lock::Bool=true,
                          branch_window::Integer=1,
                          use_jacobian_predictor::Bool=true,
                          evaluator::Function=(mm, nn, cc) -> eigenvalue(mm, nn, cc; spheroid=spheroid, precision=precision))

    _validate_precision(precision)
    if spheroid != :prolate && spheroid != :oblate
        error("spheroid must be :prolate or :oblate, got :$spheroid")
    end
    if length(c_grid) == 0
        error("c_grid must be non-empty")
    end
    if branch_window < 0
        error("branch_window must be nonnegative, got $branch_window")
    end

    T = precision === :quad ? BigFloat : Float64
    cvals = T.(c_grid)
    all(isfinite, cvals) || error("c_grid must be finite")
    if length(cvals) > 1
        diffs = diff(cvals)
        if !all(isfinite.(diffs)) || all(d -> d == 0.0, diffs)
            error("c_grid must be strictly monotone and finite")
        end
        increasing = all(d -> d > 0.0, diffs)
        decreasing = all(d -> d < 0.0, diffs)
        if !(increasing || decreasing)
            error("c_grid must be strictly monotone (all increasing or all decreasing)")
        end
    end

    npts = length(cvals)
    lambdas = Vector{T}(undef, npts)
    selected_n = Vector{Int}(undef, npts)
    switched_branch = falses(npts)

    lambda0 = evaluator(m, n, cvals[1])
    if !(lambda0 isa Real) || !isfinite(T(lambda0))
        error("evaluator returned non-finite or non-real eigenvalue at first grid point")
    end
    lambdas[1] = T(lambda0)
    selected_n[1] = Int(n)

    if npts == 1
        return (c=cvals, lambda=lambdas, selected_n=selected_n, switched_branch=switched_branch)
    end

    for i in 2:npts
        ci = cvals[i]
        cprev = cvals[i - 1]
        dc = ci - cprev

        lambda_pred = if i == 2
            if use_jacobian_predictor
                try
                    jac = jacobian_eigen(m, selected_n[i - 1], cprev;
                                         spheroid=spheroid,
                                         precision=precision,
                                         with_metadata=true,
                                         adaptive=true)
                    dlambda = T(jac.derivative)
                    if isfinite(dlambda) && jac.metadata.suggested_action == :accept
                        lambdas[i - 1] + dlambda * dc
                    else
                        lambdas[i - 1]
                    end
                catch
                    lambdas[i - 1]
                end
            else
                lambdas[i - 1]
            end
        else
            slope = (lambdas[i - 1] - lambdas[i - 2]) / (cvals[i - 1] - cvals[i - 2])
            lambdas[i - 1] + slope * dc
        end

        n_center = selected_n[i - 1]
        n_min = max(Int(m), n_center - Int(branch_window))
        n_max = n_center + Int(branch_window)

        candidate_n_values = if branch_lock
            collect(n_min:n_max)
        else
            [Int(n)]
        end

        best_n = Int(n)
        best_lambda = T(NaN)
        best_score = T(Inf)

        for n_candidate in candidate_n_values
            lambda_candidate = evaluator(m, n_candidate, ci)
            if !(lambda_candidate isa Real)
                continue
            end
            lambda_value = T(lambda_candidate)
            if !isfinite(lambda_value)
                continue
            end
            score = abs(lambda_value - lambda_pred)
            if score < best_score
                best_score = score
                best_lambda = lambda_value
                best_n = n_candidate
            end
        end

        if !isfinite(best_score)
            error("could not find a finite continuation candidate at c = $ci")
        end

        lambdas[i] = best_lambda
        selected_n[i] = best_n
        switched_branch[i] = best_n != selected_n[i - 1]
    end

    return (c=cvals, lambda=lambdas, selected_n=selected_n, switched_branch=switched_branch)
end

"""
    eigenvalue_sweep(m, n, c_grid::AbstractVector{<:Complex};
                    spheroid=:prolate, precision=:double,
                    branch_lock=true, branch_window=1)

Continue an eigenvalue along the ordered, piecewise-linear complex path
`c_grid`. Initialize degree `n` by continuation from `real(first(c_grid))`
to the first point, as in scalar `eigenvalue`. Subsequent points continue from
the preceding state, with adaptive subdivision and angular-profile matching.
Repeated points and closed paths are allowed.

Only native degrees of the same parity can be selected. `branch_window` counts
neighbors in that parity class (one means degrees `n-2` and `n+2`). The returned
`selected_n` contains native labels; a label change need not be a discontinuity
of the continued eigenvalue. `switched_branch` marks those label changes.
`c` and `lambda` are complex vectors at the requested precision.

With `branch_lock=false`, evaluate the original native label independently at
each point. Complex continuation uses profiles, not the real sweep's Jacobian
predictor or custom eigenvalue-only evaluator. Unresolved paths raise an error.
Paths around branch points can return a different eigenvalue at their starting
coordinate; this function does not impose a globally single-valued branch.
"""
function eigenvalue_sweep(m::Integer,n::Integer,c_grid::AbstractVector{<:Complex};
        spheroid::Symbol=:prolate,precision::Symbol=:double,
        branch_lock::Bool=true,branch_window::Integer=1,
        use_jacobian_predictor::Bool=true,evaluator=nothing)
    _validate_precision(precision)
    spheroid in (:prolate,:oblate) || error("spheroid must be :prolate or :oblate")
    0 <= m <= n || error("require 0 <= m <= n")
    isempty(c_grid) && error("c_grid must be non-empty")
    branch_window >= 0 || error("branch_window must be nonnegative")
    evaluator === nothing || throw(ArgumentError("complex continuation requires native angular profiles; custom eigenvalue-only evaluators are unsupported"))
    T = precision === :quad ? BigFloat : Float64
    cvals = Complex{T}.(c_grid)
    all(isfinite,cvals) || error("c_grid must be finite at the requested precision")
    prefix = spheroid === :prolate ? :cprolate : :coblate
    lambdas = similar(cvals)
    selected_n = fill(Int(n),length(cvals))
    switched_branch = falses(length(cvals))
    if !branch_lock
        for i in eachindex(cvals)
            lambdas[i] = _call_complex_eigenvalue(prefix,m,n,cvals[i];precision)
        end
    else
        evaluate = _angular_phase_evaluator(prefix,m,n,precision)
        state = _initial_angular_phase(evaluate,m,n,first(cvals);branch_window)
        lambdas[1],selected_n[1] = state.lambda,state.n
        for i in 2:length(cvals)
            state = _transport_angular_phase(evaluate,cvals[i-1],state,cvals[i];branch_window)
            lambdas[i],selected_n[i] = state.lambda,state.n
            switched_branch[i] = selected_n[i] != selected_n[i-1]
        end
    end
    return (;c=cvals,lambda=lambdas,selected_n,switched_branch)
end

