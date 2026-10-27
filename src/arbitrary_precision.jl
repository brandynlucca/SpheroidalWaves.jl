# Explicit arbitrary precision starts with separation eigenpairs. Native
# routines select the mode/sign only; every magnitude is recomputed in Julia.
function _validate_precision_bits(bits)
    bits isa Integer && !(bits isa Bool) && 2 <= bits <= div(typemax(Int), 4) ||
        throw(ArgumentError("precision_bits must be an integer of at least 2 bits"))
    return Int(bits)
end

_precision_output(x::Real, bits) = BigFloat(x, _swrounding(); precision = bits)
_precision_output(x::Complex, bits) =
    complex(_precision_output(real(x), bits), _precision_output(imag(x), bits))

_working_parameter(c::Real) = _SWFloat(c)
_working_parameter(c::Complex) = complex(_SWFloat(real(c)), _SWFloat(imag(c)))

# Compare after aligning the bilinearly normalized vectors' remaining sign.
# Eigenvalues and sensitivities use mixed absolute/relative accuracy; vectors
# use their orthonormal-basis norm, not the relative error of tiny entries.
function _arbitrary_plan_change(previous, result, derivatives)
    k = min(length(previous.v), length(result.v))
    phase = real(sum(conj.(previous.v[1:k]) .* result.v[1:k])) < 0 ? -1 : 1
    result.v .*= phase
    derivatives && (result.dv .*= phase)
    vector_error = sum(abs2, result.v[1:k]-previous.v[1:k]) +
                   sum(abs2, result.v[(k + 1):end]) +
                   sum(abs2, previous.v[(k + 1):end])
    change = max(sqrt(vector_error/sum(abs2, result.v)),
        abs(result.lambda-previous.lambda)/(1+abs(result.lambda)))
    if derivatives
        derivative_error = sum(abs2, result.dv[1:k]-previous.dv[1:k]) +
                           sum(abs2, result.dv[(k + 1):end]) +
                           sum(abs2, previous.dv[(k + 1):end])
        change = max(change,
            sqrt(derivative_error)/max(one(_SWFloat), sqrt(sum(abs2, result.dv))),
            abs(result.dlambda-previous.dlambda)/(1+abs(result.dlambda)))
    end
    return change
end

function _arbitrary_truncation_plan(m, n, c, spheroid, seed, tolerance,
        max_terms, derivatives)
    z = _working_parameter(c)
    bandwidth = ceil(Int, min(abs(z), _SWFloat(max_terms)))
    terms = min(max_terms, max(16, (n-m)÷2+9+bandwidth))
    previous = nothing
    while true
        # Leave guard bits for the linear solve and near-cancelling Rayleigh sum.
        result = _solve_coefficient_mode(m, n, z, terms, spheroid,
            _working_parameter(seed); derivatives, target = 16eps(_SWFloat))
        change = previous === nothing ?
                 (iszero(z) ? zero(_SWFloat) : _SWFloat(Inf)) :
                 _arbitrary_plan_change(previous, result, derivatives)
        if max(change, result.tail, result.residual, result.derivative_residual) <= tolerance
            return (; result..., terms, relative_change = change, tolerance,
                phase = Ref{Union{Nothing, Int}}(nothing), converged = true)
        end
        terms == max_terms &&
            error("Requested precision did not converge within max_terms=$max_terms")
        previous = result
        terms = min(max_terms, terms+max(16, terms÷2))
    end
end

const _arbitrary_coefficient_cache_key = gensym(:arbitrary_coefficient_plans)

function _arbitrary_coefficient_plan(m, n, c, spheroid, precision, bits,
        rtol, max_terms, derivatives)
    cache = get!(task_local_storage(), _arbitrary_coefficient_cache_key) do
        Dict{Any, Any}()
    end
    # Preserve the original input, including BigFloat bits beyond the ambient
    # precision. Native mode selection also depends on the ambient context.
    key = (m, n, typeof(c), c, spheroid, precision, bits, rtol, max_terms,
        derivatives, _swrounding(), Base.precision(BigFloat), rounding(BigFloat),
        backend_library(; precision))
    haskey(cache, key) && return cache[key]

    seed = if _use_small_parameter_expansion(c)
        _with_swprecision(bits+64) do
            z = _working_parameter(c)
            spherical = _SWFloat(n)*(n+1)
            z isa Real ? spherical : complex(spherical, zero(_SWFloat))
        end
    else
        eigenvalue(m, n, c; spheroid, precision)
    end
    isfinite(seed) || error("Native eigenvalue seed is not finite")
    # A tighter rtol must also increase working precision, rather than merely
    # making a fixed-precision solver fail its stopping test.
    workbits = _with_swprecision(max(192, bits+64)) do
        requested = rtol === nothing ? bits : max(bits, ceil(Int, -log2(_SWFloat(rtol))))
        max(192, requested+64)
    end
    previous = nothing
    for attempt in 1:4
        result, stable = _with_swprecision(workbits) do
            tolerance = ldexp(one(_SWFloat), -bits-8)
            rtol === nothing || (tolerance = min(tolerance, _SWFloat(rtol)/256))
            result = _arbitrary_truncation_plan(m, n, c, spheroid, seed, tolerance,
                max_terms, derivatives)
            # Guard against convergence to a different native-selected branch.
            if !_use_small_parameter_expansion(c)
                seed_tolerance = max(precision === :quad ? _SWFloat("1e-24") :
                                     _SWFloat("1e-9"), 32_SWFloat(eps(BigFloat)))
                abs(result.lambda-seed)/(1+abs(seed)) <= seed_tolerance ||
                    error("Refined eigenvalue disagrees with the selected native mode")
            end
            stable = previous !== nothing &&
                     _arbitrary_plan_change(previous, result, derivatives) <= tolerance
            result, stable
        end
        if stable
            length(cache) >= 32 && empty!(cache)
            cache[key] = result
            return result
        end
        previous = result
        workbits += max(64, bits÷2)
    end
    error("Requested precision could not be confirmed by working-precision refinement")
end

function _arbitrary_dmn(m, n, c; spheroid, precision, normalize, rtol,
        max_terms, derivatives, basis, precision_bits)
    bits = _validate_precision_bits(precision_bits)
    plan = _arbitrary_coefficient_plan(m, n, c, spheroid, precision, bits,
        rtol, max_terms, derivatives)
    return _with_swprecision(Base.precision(real(plan.lambda))) do
        phase = _use_small_parameter_expansion(c) ?
                (real(plan.v[(n-m)÷2+1]) < 0 ? -1 : 1) :
                _coefficient_phase(plan, precision)
        scale = phase*(normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(m, n, _SWFloat)))
        factors = basis === :orthonormal ? fill(scale, plan.terms) :
                  [scale/sqrt(_ferrers_norm2(m, l, _SWFloat)) for l in plan.degrees]
        # Explicit destinations preserve requested output bits without changing
        # Julia's process-wide defaults (including on Julia 1.10 and 1.11).
        output(x) = _precision_output(x, bits)
        (; coefficients = output.(plan.v .* factors),
            dcoefficients_dc = derivatives ? output.(plan.dv .* factors) : nothing,
            degrees = copy(plan.degrees), eigenvalue = output(plan.lambda),
            deigenvalue_dc = derivatives ? output(plan.dlambda) : nothing,
            terms = plan.terms, tail_norm = output(plan.tail),
            relative_change = output(plan.relative_change), converged = true)
    end
end
