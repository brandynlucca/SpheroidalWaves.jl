function _bisect_wave(evaluate, a, b, fa, fb, tolerance; maxiter = 512)
    iszero(fa) && return a
    iszero(fb) && return b
    for _ in 1:maxiter
        midpoint=(a+b)/2
        (midpoint == a || midpoint == b || b-a <= tolerance*max(one(a), abs(midpoint))) &&
            return midpoint
        fm=evaluate([midpoint])[1]
        isfinite(fm) || error("Nonfinite value while locating a coordinate zero")
        iszero(fm) && return midpoint
        if signbit(fm) == signbit(fa)
            a, fa=midpoint, fm
        else
            b, fb=midpoint, fm
        end
    end
    error("Coordinate zero did not converge")
end

function _coordinate_zeros(
        evaluate, a, b, tolerance, initial, max_points; expected = nothing, angular = false, maxiter = 512)
    count=initial
    previous=nothing
    while count <= max_points
        grid = angular ? [(a+b)/2-(b-a)/2*cospi(typeof(a)(k)/count) for k in 0:count] :
               collect(range(a, b; length = count+1))
        values=evaluate(grid)
        all(isfinite, values) || error("Nonfinite values in coordinate zero search")
        all(iszero, values) &&
            throw(DomainError((a, b), "function is identically zero on the sampled interval"))
        roots=typeof(a)[]
        for i in 1:count
            if iszero(values[i])
                push!(roots, grid[i])
            elseif !iszero(values[i + 1]) && signbit(values[i]) != signbit(values[i + 1])
                push!(roots, _bisect_wave(
                    evaluate, grid[i], grid[i + 1], values[i], values[i + 1], tolerance; maxiter))
            end
        end
        iszero(values[end]) && push!(roots, b)
        count_ok = expected === nothing || length(roots) == expected
        if count_ok && previous !== nothing && length(roots) == length(previous) &&
           all(abs.(roots-previous) .<= 4tolerance .* max.(one(a), abs.(roots)))
            return roots
        end
        previous=roots
        count*=2
    end
    error("Coordinate zero search did not stabilize within max_points=$max_points")
end

"""
    angular_zeros(m, n, c::Real; stationary=false, spheroid=:prolate,
                  precision=:double, rtol=nothing, max_points=16384)

Locate real angular zeros strictly inside (-1,1), or zeros of the first
derivative with `stationary=true`. A Chebyshev grid is refined until the roots
stabilize. Function-zero searches additionally require exactly n-m roots.
An unresolved search throws an error. Constant functions have no isolated
stationary points and are rejected for a stationary-point search.
`precision_bits=p` returns p-bit roots after precision refinement; `rtol` defaults
to `2^(-p-8)`, and `max_terms` bounds the coefficient expansion.
"""
function angular_zeros(m::Integer, n::Integer, c::Real; stationary::Bool = false,
        spheroid::Symbol = :prolate, precision::Symbol = :double, rtol = nothing, max_points::Integer = 16384,
        precision_bits = nothing, max_terms::Integer = 512)
    _validate_wave_arguments(m, n, c, [0], spheroid, precision, :angular)
    stationary && m == n == 0 && iszero(c) &&
        throw(DomainError(c, "constant function has no isolated stationary points"))
    precision_bits === nothing || return _coordinate_zeros_at_precision(
        m, n, c, (-1, 1), :angular, spheroid, precision, 1, stationary,
        rtol, max_points, precision_bits, max_terms)
    T=precision === :quad ? BigFloat : Float64
    tolerance=rtol === nothing ? T(10)^(precision === :quad ? -28 : -12) :
              _input_float(T, rtol)
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("rtol must be finite and positive"))
    margin=precision === :quad ? T(2)^(-108) : 16eps(T)
    !stationary && n == m && return T[]
    evaluate(x) = getproperty(
        smn(m, n, c, x; spheroid, precision, normalize = true, scaled = true),
        stationary ? :derivative : :value).mantissa
    initial=max(32, 4(n-m+1), ceil(Int, 4abs(c)))
    max_points >= 2initial ||
        throw(ArgumentError("max_points must allow two search grids (at least $(2initial))"))
    return _coordinate_zeros(
        evaluate, -one(T)+margin, one(T)-margin, tolerance, initial, max_points;
        expected = stationary ? nothing : n-m, angular = true)
end

"""
    radial_zeros(m, n, c::Real, (a,b); kind=1, stationary=false,
                 spheroid=:prolate, precision=:double, rtol=nothing, max_points=16384)

Locate zeros of real radial kinds 1 or 2 (or of their coordinate derivatives)
in a finite closed interval. Require 1<a<b for prolate, 0<=a<b for oblate.
Two refined grids must agree before returning roots; bisection refines each
sign change. Grid stability is a numerical check, not a proof of completeness.
Use a larger `max_points` for long intervals or highly oscillatory solutions.
`precision_bits=p` returns p-bit roots after precision refinement; `rtol` defaults
to `2^(-p-8)`, and `max_terms` bounds the coefficient expansion.
"""
function radial_zeros(m::Integer, n::Integer, c::Real, interval::Tuple{<:Real, <:Real};
        kind::Integer = 1, stationary::Bool = false, spheroid::Symbol = :prolate,
        precision::Symbol = :double, rtol = nothing, max_points::Integer = 16384,
        precision_bits = nothing, max_terms::Integer = 512)
    kind in (1, 2) ||
        throw(ArgumentError("real zero searches support radial kinds 1 and 2"))
    _validate_wave_arguments(m, n, c, collect(interval), spheroid, precision, :radial; kind)
    precision_bits === nothing || return _coordinate_zeros_at_precision(
        m, n, c, interval, :radial, spheroid, precision, kind, stationary,
        rtol, max_points, precision_bits, max_terms)
    T=precision === :quad ? BigFloat : Float64
    a, b=_input_float.(T, interval)
    a < b && (spheroid !== :prolate || a > 1) ||
        throw(ArgumentError("require a nonsingular interval with a < b"))
    tolerance=rtol === nothing ? T(10)^(precision === :quad ? -28 : -12) :
              _input_float(T, rtol)
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("rtol must be finite and positive"))
    evaluate(x) = real.(getproperty(
        rmn(m, n, c, x; spheroid, precision, kind, scaled = true),
        stationary ? :derivative : :value).mantissa)
    initial=max(32, 4(n-m+1), ceil(Int, 16abs(c)*(b-a)/T(pi)))
    max_points >= 2initial ||
        throw(ArgumentError("max_points must allow two search grids (at least $(2initial))"))
    return _coordinate_zeros(evaluate, a, b, tolerance, initial, max_points)
end

function _coordinate_zeros_at_precision(m, n, c, interval, target, spheroid, precision,
        kind, stationary, rtol, max_points, bits, max_terms)
    bits = _validate_precision_bits(bits)
    _validate_coefficient_options(m, n, c, spheroid, precision, rtol, max_terms, :orthonormal)
    workbits = max(bits+64, _precision_input_bits(c)+32,
        maximum(_precision_input_bits, interval)+32)
    if rtol!==nothing
        workbits = _with_swprecision(workbits) do
            tolerance = _SWFloat(rtol)
            isfinite(tolerance) && tolerance>0 || throw(ArgumentError("rtol must be finite and positive"))
            max(workbits, ceil(Int, -log2(tolerance))+64)
        end
    end
    previous = nothing
    for attempt in 1:4
        roots, stable = _with_swprecision(workbits) do
            parameter = _SWFloat(c)
            a, b = _SWFloat.(interval)
            a < b && (target===:angular || spheroid!==:prolate || a>1) ||
                throw(ArgumentError("require a nonsingular interval with a < b"))
            tolerance = rtol===nothing ? ldexp(one(_SWFloat), -bits-8) : _SWFloat(rtol)
            isfinite(tolerance) && tolerance>0 || throw(ArgumentError("rtol must be finite and positive"))
            if target===:angular
                margin = ldexp(one(_SWFloat), -workbits+8)
                a += margin
                b -= margin
            end
            initial = max(32, 4(n-m+1), ceil(Int, target===:angular ?
                4abs(parameter) : 16abs(parameter)*(b-a)/_SWFloat(pi)))
            max_points>=2initial || throw(ArgumentError("max_points must allow two search grids (at least $(2initial))"))
            !stationary && target===:angular && n==m && return (_SWFloat[], true)
            pair = _solve_coefficient_eigenpair_at_precision(m, n, c, spheroid, precision,
                workbits, nothing, max_terms, target===:radial)
            evaluate(points) = _with_swprecision(Base.precision(real(pair.lambda))) do
                values = if target===:angular
                    _angular_values_at_precision(pair, points, 1, true, 1;
                        rtol = ldexp(one(_SWFloat), -workbits), max_terms)
                else
                    data = _radial_states_at_precision(pair, points, kind, :standard,
                        ldexp(one(_SWFloat), -workbits), 1//attempt)
                    (; value = [s[1] for s in data.states], derivative = [s[2] for s in data.states])
                end
                real.(getproperty(values, stationary ? :derivative : :value))
            end
            found = _coordinate_zeros(evaluate, a, b, tolerance, initial, max_points;
                angular = target===:angular, expected = target===:angular && !stationary ? n-m : nothing,
                maxiter = workbits+512)
            stable = previous!==nothing && length(found)==length(previous) &&
                all(abs(x-y)<=4tolerance*max(one(x), abs(x)) for (x, y) in zip(found, previous))
            found, stable
        end
        stable && return _round_to_bigfloat.(roots, bits)
        previous = roots
        workbits += max(64, bits÷2)
    end
    error("Coordinate zeros did not stabilize under precision refinement")
end
