# Native routines select the mode and sign; Julia recomputes each magnitude.
function _validate_precision_bits(bits)
    bits isa Integer && !(bits isa Bool) && 2 <= bits <= div(typemax(Int), 4) ||
        throw(ArgumentError("precision_bits must be an integer of at least 2 bits"))
    return Int(bits)
end

_round_to_bigfloat(x::Real, bits) = BigFloat(x, _swrounding(); precision = bits)
_round_to_bigfloat(x::Complex, bits) =
    complex(_round_to_bigfloat(real(x), bits), _round_to_bigfloat(imag(x), bits))

_SWFloat_complex(c::Complex) = complex(_SWFloat(real(c)), _SWFloat(imag(c)))

# Align signs, then compare vector norms and scaled eigenvalue sensitivities.
function _align_and_measure_eigenpair_change!(previous, result, derivatives)
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

function _refine_coefficient_truncation(m, n, c, spheroid, seed, tolerance,
        max_terms, derivatives)
    z = c isa Real ? _SWFloat(c) : _SWFloat_complex(c)
    eigenvalue_seed = seed isa Real ? _SWFloat(seed) : _SWFloat_complex(seed)
    bandwidth = ceil(Int, min(abs(z), _SWFloat(max_terms)))
    terms = min(max_terms, max(16, (n-m)÷2+9+bandwidth))
    previous = nothing
    while true
        # Target requested accuracy above the matrix sum's rounding floor.
        result = _solve_coefficient_mode(m, n, z, terms, spheroid,
            eigenvalue_seed; derivatives, target = tolerance/16)
        change = previous === nothing ?
                 (iszero(z) ? zero(_SWFloat) : _SWFloat(Inf)) :
                 _align_and_measure_eigenpair_change!(previous, result, derivatives)
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

const _coefficient_eigenpair_cache_key = gensym(:coefficient_eigenpairs)

function _solve_coefficient_eigenpair_at_precision(m, n, c, spheroid, precision, bits,
        rtol, max_terms, derivatives)
    cache = get!(task_local_storage(), _coefficient_eigenpair_cache_key) do
        Dict{Any, Any}()
    end
    # Cache original inputs and the arithmetic settings used for native seeds.
    key = (m, n, typeof(c), c, spheroid, precision, bits, rtol, max_terms,
        derivatives, _swrounding(), Base.precision(BigFloat), rounding(BigFloat),
        backend_library(; precision))
    haskey(cache, key) && return cache[key]

    seed = if _use_small_parameter_expansion(c)
        _with_swprecision(bits+64) do
            spherical = _SWFloat(n)*(n+1)
            c isa Real ? spherical : complex(spherical, zero(_SWFloat))
        end
    else
        # Native bandwidths must be nonnegative; Julia retains signed sensitivities.
        eigenvalue(m, n, c isa Real ? abs(c) : c; spheroid, precision)
    end
    isfinite(seed) || error("Native eigenvalue seed is not finite")
    # Tighter tolerances also need more working bits.
    workbits = _with_swprecision(max(192, bits+64)) do
        requested = rtol === nothing ? bits : max(bits, ceil(Int, -log2(_SWFloat(rtol))))
        max(192, requested+64)
    end
    previous = nothing
    for attempt in 1:4
        result, stable = _with_swprecision(workbits) do
            tolerance = ldexp(one(_SWFloat), -bits-8)
            rtol === nothing || (tolerance = min(tolerance, _SWFloat(rtol)/256))
            result = _refine_coefficient_truncation(m, n, c, spheroid, seed, tolerance,
                max_terms, derivatives)
            # Guard against convergence to a different native-selected branch.
            if !_use_small_parameter_expansion(c)
                seed_tolerance = max(precision === :quad ? _SWFloat("1e-24") :
                                     _SWFloat("1e-9"), 32_SWFloat(eps(BigFloat)))
                abs(result.lambda-seed)/(1+abs(seed)) <= seed_tolerance ||
                    error("Refined eigenvalue disagrees with the selected native mode")
            end
            stable = previous !== nothing &&
                     _align_and_measure_eigenpair_change!(previous, result, derivatives) <= tolerance
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

function _dmn_at_precision(m, n, c; spheroid, precision, normalize, rtol,
        max_terms, derivatives, basis, precision_bits)
    bits = _validate_precision_bits(precision_bits)
    eigenpair = _solve_coefficient_eigenpair_at_precision(m, n, c, spheroid, precision, bits,
        rtol, max_terms, derivatives)
    return _with_swprecision(Base.precision(real(eigenpair.lambda))) do
        phase = _use_small_parameter_expansion(c) ?
                (real(eigenpair.v[(n-m)÷2+1]) < 0 ? -1 : 1) :
                _coefficient_phase(c isa Real ? (; eigenpair..., c = abs(eigenpair.c)) : eigenpair,
            precision)
        scale = phase*(normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(m, n, _SWFloat)))
        factors = basis === :orthonormal ? fill(scale, eigenpair.terms) :
                  [scale/sqrt(_ferrers_norm2(m, l, _SWFloat)) for l in eigenpair.degrees]
        # Explicit destinations preserve requested bits without changing global settings.
        round_to_bigfloat(x) = _round_to_bigfloat(x, bits)
        (; coefficients = round_to_bigfloat.(eigenpair.v .* factors),
            dcoefficients_dc = derivatives ? round_to_bigfloat.(eigenpair.dv .* factors) : nothing,
            degrees = copy(eigenpair.degrees), eigenvalue = round_to_bigfloat(eigenpair.lambda),
            deigenvalue_dc = derivatives ? round_to_bigfloat(eigenpair.dlambda) : nothing,
            terms = eigenpair.terms, tail_norm = round_to_bigfloat(eigenpair.tail),
            relative_change = round_to_bigfloat(eigenpair.relative_change), converged = true)
    end
end
