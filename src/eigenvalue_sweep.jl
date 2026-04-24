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
- `c::Vector{Float64}`
- `lambda::Vector{Float64}`
- `selected_n::Vector{Int}`
- `switched_branch::Vector{Bool}`
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

    cvals = Float64.(c_grid)
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
    lambdas = Vector{Float64}(undef, npts)
    selected_n = Vector{Int}(undef, npts)
    switched_branch = falses(npts)

    lambda0 = evaluator(m, n, cvals[1])
    if !(lambda0 isa Real) || !isfinite(Float64(lambda0))
        error("evaluator returned non-finite or non-real eigenvalue at first grid point")
    end
    lambdas[1] = Float64(lambda0)
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
                    dlambda = Float64(jac.derivative)
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
        best_lambda = NaN
        best_score = Inf

        for n_candidate in candidate_n_values
            lambda_candidate = evaluator(m, n_candidate, ci)
            if !(lambda_candidate isa Real)
                continue
            end
            lambda_float = Float64(lambda_candidate)
            if !isfinite(lambda_float)
                continue
            end
            score = abs(lambda_float - lambda_pred)
            if score < best_score
                best_score = score
                best_lambda = lambda_float
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

