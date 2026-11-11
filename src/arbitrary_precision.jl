# _solve_coefficient_eigenpair_at_precision selects MPFR modes before refinement.
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
        max_terms, derivatives; min_terms = 0)
    0<=min_terms<=max_terms || throw(ArgumentError("min_terms must be between zero and max_terms"))
    z = c isa Real ? _SWFloat(c) : _SWFloat_complex(c)
    eigenvalue_seed = seed isa Real ? _SWFloat(seed) : _SWFloat_complex(seed)
    bandwidth = ceil(Int, min(abs(z), _SWFloat(max_terms)))
    terms = min(max_terms, max(16, (n-m)÷2+9+bandwidth, min_terms))
    previous = if terms == max_terms && terms>(n-m)÷2+1 && !iszero(z)
        _solve_coefficient_mode(m, n, z, terms-1, spheroid,
            eigenvalue_seed; derivatives, target = tolerance/16)
    else
        nothing
    end
    while true
        # Target requested accuracy above the matrix sum's rounding floor.
        result = _solve_coefficient_mode(m, n, z, terms, spheroid,
            eigenvalue_seed; derivatives, target = tolerance/16)
        change = previous === nothing ?
                 (iszero(z) ? zero(_SWFloat) : _SWFloat(Inf)) :
                 _align_and_measure_eigenpair_change!(previous, result, derivatives)
        if max(change, result.tail, result.residual, result.derivative_residual) <= tolerance
            return (; result..., terms, relative_change = change, tolerance,
                phase = Ref{Union{Nothing, Int}}(nothing),
                legendre_factors = Ref{Union{Nothing, Tuple{Int, Base.MPFR.MPFRRoundingMode,
                    Vector{Tuple{_SWFloat, _SWFloat}}}}}(nothing), converged = true)
        end
        terms == max_terms &&
            error("Requested precision did not converge within max_terms=$max_terms")
        previous = result
        terms = min(max_terms, terms+max(16, terms÷2))
    end
end

const _coefficient_eigenpair_cache_key = gensym(:coefficient_eigenpairs)
const _coefficient_eigenpair_cache_limit = 16*1024^2

function _coefficient_eigenpair_bytes(result)
    stored = get(result, :cache_bytes, nothing)
    stored===nothing ? _coefficient_eigenpair_bytes(result, nothing) : stored
end

function _coefficient_eigenpair_bytes(result, key)
    try
        # Reserve both tuple storage and array slots, including Julia-version-dependent BigFloat storage.
        element_bytes = Base.Checked.checked_add(
            Base.Checked.checked_mul(2, Base.summarysize(real(result.lambda))),
            2sizeof(Tuple{_SWFloat, _SWFloat}))
        count = max(0, last(result.degrees)-result.m)
        factor_bytes = Base.Checked.checked_add(128, Base.Checked.checked_mul(count, element_bytes))
        retained = Base.summarysize((key, result))
        Base.Checked.checked_add(retained, Base.Checked.checked_add(factor_bytes, 1024))
    catch exception
        exception isa OverflowError || rethrow()
        typemax(Int)
    end
end

function _cached_coefficient_eigenpair(cache, key)
    result = get(cache, key, nothing)
    result===nothing || (result.last_used[] = time_ns())
    return result
end

function _cache_coefficient_eigenpair!(cache, key, result; native_label = false)
    limit = _coefficient_eigenpair_cache_limit
    result = (; result..., native_label, last_used = Ref(time_ns()))
    required = _coefficient_eigenpair_bytes(result, key)
    result = (; result..., cache_bytes = required)
    required>limit && return result
    retained = sum(_coefficient_eigenpair_bytes, values(cache); init = 0)
    while !isempty(cache) && (length(cache)>=32 || retained>limit-required)
        victim = first(cache)
        for candidate in cache
            if (candidate.second.native_label && !victim.second.native_label) ||
               (candidate.second.native_label==victim.second.native_label &&
                candidate.second.last_used[]<victim.second.last_used[])
                victim = candidate
            end
        end
        native_label && !victim.second.native_label && return result
        removed = pop!(cache, victim.first)
        retained -= _coefficient_eigenpair_bytes(removed)
    end
    cache[key] = result
    return result
end

function _solve_coefficient_eigenpair_at_precision(m, n, c, spheroid, precision, bits,
        rtol, max_terms, derivatives; native_label = false, min_terms = 0, max_refinements = 4)
    0<=min_terms<=max_terms || throw(ArgumentError("min_terms must be between zero and max_terms"))
    cache = get!(task_local_storage(), _coefficient_eigenpair_cache_key) do
        Dict{Any, Any}()
    end
    # Cache original inputs, requested precision, rounding, and the MPFR library.
    key = (m, n, typeof(c), c, spheroid, precision, bits, rtol, max_terms, min_terms, native_label,
        derivatives, _swrounding(), Base.precision(BigFloat), rounding(BigFloat),
        backend_library(; precision = :mpfr), max_refinements)
    cached = _cached_coefficient_eigenpair(cache, key)
    cached===nothing || return cached

    mode_state = if c isa Complex && !native_label
        _with_swprecision(max(192, bits+96)) do
            evaluate = _mpfr_angular_phase_evaluator(m, n, spheroid, precision, bits, max_terms)
            _initial_angular_phase(evaluate, m, n, _SWFloat_complex(c))
        end
    else
        nothing
    end
    seed = if mode_state !== nothing
        mode_state.lambda
    elseif _use_small_parameter_expansion(c)
        _with_swprecision(bits+64) do
            spherical = _SWFloat(n)*(n+1)
            c isa Real ? spherical : complex(spherical, zero(_SWFloat))
        end
    else
        _with_swprecision(bits+64) do
            _mpfr_eigenvalue(m, n, c, spheroid, bits+64)
        end
    end
    isfinite(seed) || error("MPFR eigenvalue seed is not finite")
    # Tighter tolerances also need more working bits.
    workbits = _with_swprecision(max(192, bits+64)) do
        requested = rtol === nothing ? bits : max(bits, ceil(Int, -log2(_SWFloat(rtol))))
        max(192, requested+64)
    end
    previous = nothing
    for attempt in 1:max_refinements
        result, stable = _with_swprecision(workbits) do
            tolerance = ldexp(one(_SWFloat), -bits-8)
            rtol === nothing || (tolerance = min(tolerance, _SWFloat(rtol)/256))
            minimum_terms = mode_state === nothing ? min_terms : max(min_terms, (mode_state.n-m)÷2+1)
            result = _refine_coefficient_truncation(m, n, c, spheroid, seed, tolerance,
                max_terms, derivatives; min_terms = minimum_terms)
            # _refine_coefficient_truncation must retain the selected MPFR mode.
            if !_use_small_parameter_expansion(c)
                seed_tolerance = ldexp(one(_SWFloat), -min(bits, 64))
                abs(result.lambda-seed)/(1+abs(seed)) <= seed_tolerance ||
                    error("Refined eigenvalue disagrees with the selected MPFR mode")
            end
            stable = previous !== nothing &&
                     _align_and_measure_eigenpair_change!(previous, result, derivatives) <= tolerance
            (; result..., mode_state), stable
        end
        if stable
            return _cache_coefficient_eigenpair!(cache, key, result; native_label)
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
        phase = _coefficient_phase_at_precision(eigenpair)
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
