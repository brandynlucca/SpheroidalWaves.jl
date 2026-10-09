function _legendre_quadrature(count, T)
    nodes, weights = zeros(T, count), zeros(T, count)
    for i in 1:cld(count, 2)
        z = cospi((T(i)-T(1)/4)/(count+T(1)/2))
        for iteration in 1:50
            p0, p1 = one(T), z
            for k in 2:count
                p0, p1 = p1, ((2k-1)*z*p1-(k-1)*p0)/k
            end
            derivative = count*(z*p1-p0)/(z^2-1)
            delta = p1/derivative
            z -= delta
            abs(delta) <= 4eps(T)*max(one(T), abs(z)) && break
            iteration == 50 && error("Legendre quadrature did not converge")
        end
        # Re-evaluate the derivative at the converged node for the weight.
        p0, p1 = one(T), z
        for k in 2:count
            p0, p1 = p1, ((2k-1)*z*p1-(k-1)*p0)/k
        end
        derivative = count*(z*p1-p0)/(z^2-1)
        w = 2/((1-z)*(1+z)*derivative^2)
        nodes[i], nodes[count + 1 - i] = -z, z
        weights[i]=weights[count + 1 - i]=w
    end
    return nodes, weights
end

function _ferrers_norm2(m, n, T)
    value = T(2)/(2T(n)+1)
    for k in (n - m + 1):(n + m)
        value *= k
    end
    return value
end

# Multiplication by x in the orthonormal Ferrers basis has entries a_l.
# H = diag(l(l+1)) + sigma*c^2*x^2 uses the package's separation constant.
# For complex c this is symmetric under transpose, not Hermitian.
function _coefficient_matrix(m, n, c, terms, spheroid)
    degrees = collect((m + mod(n - m, 2)):2:(m + mod(n - m, 2) + 2terms - 2))
    a(l) = l < m ? zero(_SWFloat) :
           sqrt((_SWFloat(l+1)^2-_SWFloat(m)^2) /
                ((2_SWFloat(l)+1)*(2_SWFloat(l)+3)))
    b = [a(l)^2+a(l-1)^2 for l in degrees]
    e = [a(l)*a(l+1) for l in degrees[1:(end - 1)]]
    sigma = spheroid === :prolate ? 1 : -1
    diagonal = [_SWFloat(l)*(l+1) for l in degrees] .+ sigma*c^2 .* b
    off = sigma*c^2 .* e
    return degrees, Tridiagonal(off, diagonal, copy(off)),
    Tridiagonal(2sigma*c .* e, 2sigma*c .* b, 2sigma*c .* e)
end

function _solve_coefficient_mode(m, n, c, terms, spheroid, seed)
    degrees, H, dH = _coefficient_matrix(m, n, c, terms, spheroid)
    v = zeros(eltype(H), terms)
    v[min(terms, (n-m)÷2+1)] = 1
    target = eps(_SWFloat)^(2//3)
    shift = target^(1//2)*(1+abs(seed))
    shifted = Tridiagonal(copy(H.dl), H.d .- (seed+shift), copy(H.du))
    lambda = seed
    residual = _SWFloat(Inf)
    max_iterations = max(16, cld(Base.precision(_SWFloat), 32))
    for iteration in 1:max_iterations
        w = iszero(c) ? v : shifted \ v
        norm2 = sum(w .* w)
        abs(norm2) > sqrt(eps(_SWFloat))*sum(abs2, w) ||
            error("Coefficient eigenvector is nearly self-orthogonal; the mode is unresolved")
        v = w ./ sqrt(norm2)
        product = H*v
        lambda = sum(v .* product)
        residual = maximum(abs, product-lambda*v)/((1+abs(lambda))*maximum(abs, v))
        residual <= target && break
        iteration == max_iterations && error("Coefficient eigenvector did not converge")
    end
    # Hellmann-Feynman with the analytic bilinear normalization v^T v = 1.
    dlambda = sum(v .* (dH*v))
    rhs = dlambda*v-dH*v
    pivot = argmax(abs.(v))
    dv = zero(v)
    # Pin one nonzero component, solve the two tridiagonal blocks, then impose
    # v^T dv = 0. This avoids an ill-conditioned shifted inverse for sensitivities.
    for indices in (1:(pivot - 1), (pivot + 1):terms)
        isempty(indices) && continue
        lo, hi = first(indices), last(indices)
        dv[indices] = Tridiagonal(H.dl[lo:(hi - 1)], H.d[indices] .- lambda, H.du[lo:(hi - 1)]) \
                      rhs[indices]
    end
    dv .-= v*sum(v .* dv)
    derivative_residual = maximum(abs, (H*dv-lambda*dv)-rhs) /
                          max(one(_SWFloat), maximum(abs, rhs), (1+abs(lambda))*maximum(abs, dv))
    tail = max(sqrt(sum(abs2, v[max(1, end - 3):end]))/sqrt(sum(abs2, v)),
        sqrt(sum(abs2, dv[max(1, end - 3):end]))/max(one(_SWFloat), sqrt(sum(abs2, dv))))
    iszero(c) && (tail = zero(_SWFloat))
    eigenvalue_condition = sum(abs2, v)/abs(sum(v .* v))
    return (; m, n, c, spheroid, degrees, v, dv, lambda, dlambda,
        residual, derivative_residual, tail, eigenvalue_condition)
end

# Bounded task-local storage. Include arithmetic settings and the backend used
# to select the eigenvalue branch. Cached vectors are never returned to callers.
function _coefficient_plan(
        m, n, c; spheroid = :prolate, precision = :double, rtol = nothing,
        max_terms = 512, eigenvalue_seed = nothing, min_terms = 0)
    _validate_wave_arguments(m, n, c, [0], spheroid, precision, :angular)
    tolerance = rtol === nothing ?
                (precision === :quad ? _SWFloat("1e-29") : _SWFloat("1e-14")) :
                _SWFloat(rtol)
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("rtol must be finite and positive"))
    minimum_terms = max((n-m)÷2+1, min_terms)
    max_terms >= minimum_terms ||
        throw(ArgumentError("max_terms must include the requested degree"))
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? _input_float(T, c) : _input_float(Complex{T}, c)
    cache = get!(task_local_storage(), :SpheroidalWaves_coefficient_plans) do
        Dict{Any, Any}()
    end
    key = (m, n, typeof(parameter), parameter, spheroid, precision, tolerance,
        max_terms, eigenvalue_seed, min_terms, Base.precision(_SWFloat),
        rounding(_SWFloat), backend_library(; precision))
    haskey(cache, key) && return cache[key]
    seed = eigenvalue_seed === nothing ? eigenvalue(m, n, parameter; spheroid, precision) :
           eigenvalue_seed
    plan = _with_swprecision(max(192, Base.precision(_SWFloat))) do
        z = parameter isa Real ? _SWFloat(parameter) : Complex{_SWFloat}(parameter)
        terms = min(max_terms, max(16, minimum_terms+ceil(Int, abs(z))+8))
        previous = nothing
        while true
            result = _solve_coefficient_mode(m, n, z, terms, spheroid, seed)
            change = _SWFloat(Inf)
            if previous !== nothing
                k = length(previous.v)
                sign = real(sum(conj.(previous.v) .* result.v[1:k])) < 0 ? -1 : 1
                result.v .*= sign
                result.dv .*= sign
                change = max(
                    sqrt(sum(abs2, result.v[1:k]-previous.v))/sqrt(sum(abs2, result.v)),
                    sqrt(sum(abs2, result.dv[1:k]-previous.dv))/max(one(_SWFloat), sqrt(sum(abs2, result.dv))),
                    abs(result.dlambda-previous.dlambda)/max(one(_SWFloat), abs(result.dlambda)))
            elseif iszero(z)
                change = zero(_SWFloat)
            end
            converged = max(change, result.tail, result.residual, result.derivative_residual) <=
                        tolerance
            if converged || terms == max_terms
                # Confirm that refinement has retained the backend-selected mode.
                agreement = abs(result.lambda-seed)/(1+abs(seed))
                (!converged ||
                 agreement <=
                 (precision === :quad ? _SWFloat("1e-24") : _SWFloat("1e-9"))) ||
                    error("Coefficient eigenvalue disagrees with the selected backend mode")
                return (; result..., terms, relative_change = change, converged,
                    tolerance, phase = Ref{Union{Nothing, Int}}(nothing))
            end
            previous = result
            terms = min(max_terms, terms+max(16, terms÷2))
        end
    end
    if plan.converged
        length(cache) >= 32 && empty!(cache)
        cache[key] = plan
    end
    return plan
end

# Near zero, the perturbation sigma*c^2*x^2 is much smaller than the
# spherical eigenvalue gaps. Select that isolated mode without native seeds
# or angular anchors, which can be unreliable for tiny nonzero parameters.
_use_small_parameter_expansion(c) = abs(c) <= 1//10000

function _small_parameter_plan(m, n, c, spheroid, precision)
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? _input_float(T, c) : _input_float(Complex{T}, c)
    # Extra guard bits preserve relative accuracy in the O(c^2) coefficients,
    # even when the absolute eigenproblem residual is already tiny.
    return _with_swprecision(max(320, Base.precision(_SWFloat))) do
        z = parameter isa Real ? _SWFloat(parameter) : Complex{_SWFloat}(parameter)
        # Refine the spherical seed before the cross-expansion agreement check.
        seed = _solve_coefficient_mode(
            m, n, z, max(16, (n-m)÷2+9), spheroid, _SWFloat(n)*(n+1)).lambda
        plan = _coefficient_plan(m, n, parameter; spheroid, precision,
            eigenvalue_seed = seed, max_terms = max(64, (n-m)÷2+32))
        plan.converged || error("Small-parameter coefficient expansion did not converge")
        return plan
    end
end

function _small_parameter_smn(m, n, c, points, spheroid, precision, normalize)
    plan = _small_parameter_plan(m, n, c, spheroid, precision)
    # Match the native convention (without Condon–Shortley); the public
    # boundary supplies that phase. The dominant coefficient stays nonzero.
    phase = real(plan.v[(n - m) ÷ 2 + 1]) < 0 ? -1 : 1
    scale = (isodd(m) ? -phase : phase) *
            (normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(m, n, _SWFloat)))
    return _evaluate_coefficient_vector(plan, scale .* plan.v, points)
end

# Evaluate a coefficient vector and its coordinate derivative in one recurrence.
# Factoring (1-x^2)^(m/2) gives explicit endpoint limits without Inf-Inf sums.
function _evaluate_coefficient_vector(
        plan, coefficients, points; second_derivative = false, regular_factor = false)
    values = zeros(eltype(coefficients), length(points))
    derivatives = similar(values)
    seconds = similar(values)
    m = plan.m
    for (i, point) in enumerate(points)
        x = _SWFloat(point)
        p = inv(sqrt(_SWFloat("2")))
        for k in 1:m
            p *= -sqrt(_SWFloat(2k+1)/(2k))
        end
        previous, dp, dprevious = zero(p), zero(p), zero(p)
        ddp, ddprevious = zero(p), zero(p)
        f, df = zero(eltype(coefficients)), zero(eltype(coefficients))
        ddf = zero(f)
        index = 1
        for l in m:last(plan.degrees)
            if l == plan.degrees[index]
                f += coefficients[index]*p
                df += coefficients[index]*dp
                second_derivative && (ddf += coefficients[index]*ddp)
                index += 1
                index > length(coefficients) && break
            end
            k = _SWFloat(l+1)
            a = sqrt((4k^2-1)/((k-m)*(k+m)))
            b = l == m ? zero(k) : sqrt((2k+1)*(k-1-m)*(k-1+m)/((2k-3)*(k-m)*(k+m)))
            second_derivative && ((ddprevious, ddp) = (ddp, a*(2dp+x*ddp)-b*ddprevious))
            previous, p, dprevious, dp = p, a*x*p-b*previous, dp, a*(p+x*dp)-b*dprevious
        end
        u = (1-x)*(1+x)
        if regular_factor
            values[i], derivatives[i], seconds[i] = f, df, ddf
        elseif iszero(u)
            values[i] = m == 0 ? f : zero(f)
            derivatives[i] = m == 0 ? df :
                             m == 1 ?
                             (iszero(f) ? zero(f) :
                              _directed_infinity(
                -x*f, coefficients isa
                      AbstractVector{<:Real}, _SWFloat)) :
                             m == 2 ? -2x*f : zero(f)
            if second_derivative
                seconds[i] = m==0 ? ddf :
                             m==2 ? -2*f-4x*df :
                             m==4 ? 8f :
                             m>4 ? zero(f) :
                             _directed_infinity(
                    (m==1 ? -1 : 1)*f, coefficients isa AbstractVector{<:Real}, _SWFloat)
            end
        else
            factor = u^(m//2)
            values[i] = factor*f
            derivatives[i] = factor*(df-m*x/u*f)
            second_derivative &&
                (seconds[i] = factor*(ddf-2m*x/u*df+(m*(m-2)*x^2/u^2-m/u)*f))
        end
    end
    result = (; value = values, derivative = derivatives)
    return second_derivative ? (; result..., second_derivative = seconds) : result
end

function _coefficient_phase(plan, precision; mode = nothing)
    # Match the existing angular convention only once per prepared plan. Its
    # sign is locally constant, so it multiplies both v and dv without a derivative.
    mode === nothing && plan.phase[] !== nothing && return plan.phase[]
    endpoint = 1-_SWFloat(2)^(-ceil(Int, log2(32*(1+_SWFloat(plan.n)*(plan.n+1)+abs(plan.c)^2))))
    # Oblate modes at large real c concentrate near the endpoints; their origin
    # sample can be lost even while the eigenvalue and endpoint values are sound.
    points = plan.c isa Real && plan.spheroid === :oblate && abs(plan.c)>plan.n+1 ?
             (endpoint, _SWFloat("0")) : (_SWFloat("0"), endpoint)
    for x in points
        # These anchors are strictly interior. Use the native/continued
        # evaluator directly, avoiding a recursive public smn dispatch cycle.
        reference = if mode !== nothing
            prefix = plan.spheroid === :prolate ? :cprolate : :coblate
            result = _call_complex_smn_raw(
                prefix, plan.m, mode.n, plan.c, [x]; precision, normalize = true)
            _scale_mode_result!(result, mode.sign)
        elseif plan.c isa Complex && iszero(real(plan.c))
            prefix = plan.spheroid === :prolate ? :oblate : :psms
            _call_real_smn(
                prefix, plan.m, plan.n, abs(imag(plan.c)), [x]; precision, normalize = true)
        elseif plan.c isa Real
            prefix = plan.spheroid === :prolate ? :psms : :oblate
            _call_real_smn(prefix, plan.m, plan.n, plan.c, [x]; precision, normalize = true)
        else
            prefix = plan.spheroid === :prolate ? :cprolate : :coblate
            _call_complex_smn(
                prefix, plan.m, plan.n, plan.c, [x]; precision, normalize = true)
        end
        reference = _angular_phase!(reference, plan.m)
        calculated = _evaluate_coefficient_vector(plan, plan.v, [x])
        name = iszero(x) && isodd(plan.n-plan.m) ? :derivative : :value
        a, b = only(getproperty(reference, name)), only(getproperty(calculated, name))
        isfinite(a) && !iszero(a) && isfinite(b) && !iszero(b) || continue
        phase = abs(a-b) <= abs(a+b) ? 1 : -1
        mode === nothing && (plan.phase[] = phase)
        return phase
    end
    error("Angular coefficient phase cannot be resolved")
end

# Angular values can be exponentially smaller than their Legendre
# summands. Fixed native precision then loses digits even when lambda is sound.
_use_angular_expansion(n, c, spheroid) = c isa Real && abs(c)>n+1

_angular_precision(c) = max(Base.precision(_SWFloat)+64, 256+ceil(Int, 5abs(c)))

function _guarded_angular_plan(m, n, c, spheroid, precision; eigenvalue_seed = nothing)
    # Form the exponential in the working type, before it can underflow in
    # Float64 for a large ordinary (or complex-double) bandwidth.
    tolerance = min(_SWFloat("1e-40"), exp(-2_SWFloat(abs(c)))*_SWFloat("1e-40"))
    max_terms = max(512, (n-m)÷2+2ceil(Int, abs(c))+32)
    plan = _coefficient_plan(
        m, n, c; spheroid, precision, rtol = tolerance, max_terms, eigenvalue_seed)
    plan.converged || error("Angular coefficient expansion did not converge")
    return plan
end

function _coefficient_smn(m, n, c, points, spheroid, precision, normalize,
        scaled, logderivative, second_derivative)
    T = precision === :quad ? BigFloat : Float64
    parameter, coordinates = _input_float(T, c), _SWFloat.(points)
    # Allow for O(exp(-abs(c))) interior values, with guard bits left over for
    # the eigenvector solve and cancellation in values and derivatives.
    bits = _angular_precision(parameter)
    wide = _with_swprecision(bits) do
        # The truncation error must shrink along with the smallest interior
        # values; a fixed absolute coefficient tolerance is insufficient.
        plan = _guarded_angular_plan(m, n, parameter, spheroid, precision)
        scale = _coefficient_phase(plan, precision)*(normalize ? one(_SWFloat) :
                                                     sqrt(_ferrers_norm2(m, n, _SWFloat)))
        _evaluate_coefficient_vector(plan, scale .* plan.v, coordinates; second_derivative)
    end
    output = map(v -> scaled ? _decimal_scaled(v, T) : T.(v), wide)
    if logderivative
        output = (; output...,
            logderivative = T[_wave_logderivative(v, d, x, m, spheroid, :angular, 1)
                              for (v, d, x) in zip(wide.value, wide.derivative, points)])
    end
    return output
end

# Compatibility for internal callers of the original coefficient interface.
function _angular_coefficients(m::Integer, n::Integer, c::Union{Real, Complex};
        spheroid::Symbol = :prolate, precision::Symbol = :double, normalize::Bool = false,
        rtol = nothing, max_terms::Integer = 512)
    (; dmn(m, n, c; spheroid, precision, normalize, rtol, max_terms)...,
        quadrature_points = 0)
end

function _with_coefficient_expansion(
        f, m, n, c, spheroid, precision, normalize, rtol, max_terms)
    _validate_wave_arguments(m, n, c, [0], spheroid, precision, :angular)
    return _with_swprecision(max(320, _angular_precision(c))) do
        small = _use_small_parameter_expansion(c)
        seed = small ? _small_parameter_plan(m, n, c, spheroid, precision).lambda : nothing
        plan = _coefficient_plan(
            m, n, c; spheroid, precision, rtol, max_terms, eigenvalue_seed = seed)
        plan.converged || error("Coefficient expansion did not converge within max_terms")
        phase = small ? (real(plan.v[(n - m) ÷ 2 + 1]) < 0 ? -1 : 1) :
                _coefficient_phase(plan, precision)
        scale = phase*(normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(m, n, _SWFloat)))
        f(plan, scale)
    end
end

"""
    dmn(m, n, c; spheroid=:prolate, precision=:double,
                           normalize=false, rtol=nothing, max_terms=512)

Return the Ferrers expansion of `smn(m,n,c,x; kind=1, normalize, ...)`:
`S(x) = sum(coefficients[i] * P_degrees[i]^m(x))`, with Condon–Shortley phase
included in the Ferrers basis. Only degrees of the parity of `n-m` occur.
In DLMF 30.8.1 notation, the returned coefficient at degree `l=n+2k` is
`(-1)^k * a[n,k]^m`, not `a[n,k]^m` itself.

The named tuple contains `degrees`, `coefficients`, `dcoefficients_dc`,
`eigenvalue`, `deigenvalue_dc`, `terms`, `tail_norm`, `relative_change`, and
`converged`. Derivatives hold the chosen normalization fixed. Returned arrays
are independent copies. `normalize=true` uses the same unit bilinear norm as
`smn`. Complex coefficients are analytic continuations, not Hermitian-normalized.

Supports integers `0 ≤ m ≤ n`, real or complex `c`, both geometries, and `c=0`.
`precision=:double` returns Float64 values, `:quad` returns BigFloat values,
with complex counterparts for complex `c`. `rtol` controls coefficient-vector
convergence, not relative accuracy of every tiny coefficient. An unresolved
mode or exhausted `max_terms` raises an error.
"""
function dmn(m::Integer, n::Integer, c::Union{Real, Complex};
        spheroid::Symbol = :prolate, precision::Symbol = :double, normalize::Bool = false,
        rtol = nothing, max_terms::Integer = 512)
    T = precision === :quad ? BigFloat : Float64
    V = c isa Real ? T : Complex{T}
    return _with_coefficient_expansion(
        m, n, c, spheroid, precision, normalize, rtol, max_terms) do plan, scale
        factors = [scale/sqrt(_ferrers_norm2(m, l, _SWFloat)) for l in plan.degrees]
        (; coefficients = V.(plan.v .* factors),
            dcoefficients_dc = V.(plan.dv .* factors),
            degrees = copy(plan.degrees), eigenvalue = V(plan.lambda), deigenvalue_dc = V(plan.dlambda),
            terms = plan.terms, tail_norm = T(plan.tail),
            relative_change = T(plan.relative_change), converged = true)
    end
end

"""
    amn(m, n, c; spheroid=:prolate, precision=:double,
                  rtol=nothing, max_terms=512)

Return the dimensionless DLMF 30.11.4 radial factor `A_n^m(σ*c^2)` in
Meixner–Schäfke normalization, where `σ=1` for prolate and `σ=-1` for oblate.
This function accepts **signed** integer order `-n ≤ m ≤ n`. For negative order
`-p`, it returns `sum(d_l * (l+p)!/(l-p)!) * (n-p)!/(n+p)!`, where `d_l` are
the unnormalized angular coefficients of order `p`. Positive order includes
the finite extension below degree `m`, so it is generally not `sum(d_l)`.

`A_n^m(0)=1`. Both real and complex parameters are supported. The result is
Float64 or BigFloat according to `precision`, or its complex counterpart.
An unresolved coefficient expansion or singular finite extension raises an error.
See [`dmn`](@ref) for truncation controls.
"""
function amn(m::Integer, n::Integer, c::Union{Real, Complex};
        spheroid::Symbol = :prolate, precision::Symbol = :double, rtol = nothing, max_terms::Integer = 512)
    T = precision === :quad ? BigFloat : Float64
    V = c isa Real ? T : Complex{T}
    p = abs(m)
    return _with_coefficient_expansion(
        p, n, c, spheroid, precision, false, rtol, max_terms) do plan, scale
        sums = _coefficient_joining_sums(plan, scale; extend = m>=0)
        value = m < 0 ?
                sums.weighted/prod(_SWFloat(k) for k in (n - p + 1):(n + p); init = one(_SWFloat)) :
                sums.joining
        V(value)
    end
end

"""
    kmn(m, n, c; spheroid=:prolate, precision=:double,
                   normalize=false, rtol=nothing, max_terms=512)

Return the dimensionless connection coefficient `K_n^m(γ)` in DLMF 30.11.8:
`S_n^{m(1)}(z,γ) = K_n^m(γ) * Ps_n^m(z,γ^2)`. Here `Ps` is the exterior
Legendre-branch continuation, not the Ferrers value at a real radial coordinate.
For prolate `γ=c`. For oblate `γ=im*c`, and the corresponding radial coordinate
is `z=-im*x`. The oblate result is complex even for real `c`.

`normalize=true` gives the multiplier for the unit-normalized angular function.
Supports `0 ≤ m ≤ n`, real or complex `c`, and both precisions. At `c=0`, the
factor is 1 for `(m,n)=(0,0)` and 0 for other modes before optional normalization.
These limits do not define `rmn` at zero parameter. Formulas 30.11.10–11 use
origin values or derivatives according to parity and avoid division by angular
zeros at an arbitrary matching coordinate. Unresolved sums raise an error.
"""
function kmn(m::Integer, n::Integer, c::Union{Real, Complex};
        spheroid::Symbol = :prolate, precision::Symbol = :double, normalize::Bool = false,
        rtol = nothing, max_terms::Integer = 512)
    T = precision === :quad ? BigFloat : Float64
    V = c isa Real && spheroid === :prolate ? T : Complex{T}
    return _with_coefficient_expansion(
        m, n, c, spheroid, precision, normalize, rtol, max_terms) do plan, scale
        V(_joining_factor_value(plan, scale))
    end
end

function _joining_factor_value(plan, scale)
    m, n = plan.m, plan.n
    factors = [sqrt(prod(_SWFloat(k) for k in (l - m + 1):(l + m); init = one(_SWFloat))*(2l+1)/2)
               for l in plan.degrees]
    weights = factors .* plan.v
    denominator = sum(weights)
    abs(denominator)>sqrt(eps(_SWFloat))*sum(abs, weights) ||
        error("Joining factor denominator is singular or numerically unresolved")
    origin = _evaluate_coefficient_vector(plan, scale .* plan.v, [zero(_SWFloat)])
    parity = mod(n-m, 2)
    anchor = only(parity == 0 ? origin.value : origin.derivative)
    isfinite(anchor) && !iszero(anchor) ||
        error("Joining factor origin normalization is unresolved")
    degree = m+parity
    gamma = plan.spheroid === :prolate ? plan.c : im*plan.c
    odd_factorial = prod(_SWFloat(k) for k in 1:2:(2degree + 1); init = one(_SWFloat))
    (-1)^(m+(degree-n)÷2)*gamma^degree*first(weights)/(odd_factorial*denominator*anchor)
end
