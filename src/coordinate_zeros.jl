function _bisect_wave(evaluate,a,b,fa,fb,tolerance)
    iszero(fa) && return a
    iszero(fb) && return b
    for _ in 1:512
        midpoint=(a+b)/2
        (midpoint == a || midpoint == b || b-a <= tolerance*max(one(a),abs(midpoint))) && return midpoint
        fm=evaluate([midpoint])[1]
        isfinite(fm) || error("Nonfinite value while locating a coordinate zero")
        iszero(fm) && return midpoint
        if signbit(fm) == signbit(fa)
            a,fa=midpoint,fm
        else
            b,fb=midpoint,fm
        end
    end
    error("Coordinate zero did not converge")
end

function _coordinate_zeros(evaluate,a,b,tolerance,initial,max_points;expected=nothing,angular=false)
    count=initial
    previous=nothing
    while count <= max_points
        grid = angular ? [(a+b)/2-(b-a)/2*cospi(typeof(a)(k)/count) for k in 0:count] : collect(range(a,b;length=count+1))
        values=evaluate(grid)
        all(isfinite,values) || error("Nonfinite values in coordinate zero search")
        all(iszero,values) && throw(DomainError((a,b),"function is identically zero on the sampled interval"))
        roots=typeof(a)[]
        for i in 1:count
            if iszero(values[i])
                push!(roots,grid[i])
            elseif !iszero(values[i+1]) && signbit(values[i]) != signbit(values[i+1])
                push!(roots,_bisect_wave(evaluate,grid[i],grid[i+1],values[i],values[i+1],tolerance))
            end
        end
        iszero(values[end]) && push!(roots,b)
        count_ok = expected === nothing || length(roots) == expected
        if count_ok && previous !== nothing && length(roots) == length(previous) &&
                all(abs.(roots-previous) .<= 4tolerance.*max.(one(a),abs.(roots)))
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
"""
function angular_zeros(m::Integer,n::Integer,c::Real;stationary::Bool=false,
        spheroid::Symbol=:prolate,precision::Symbol=:double,rtol=nothing,max_points::Integer=16384)
    _validate_wave_arguments(m,n,c,[0],spheroid,precision,:angular)
    stationary && m == n == 0 && iszero(c) && throw(DomainError(c,"constant function has no isolated stationary points"))
    T=precision === :quad ? BigFloat : Float64
    tolerance=rtol === nothing ? T(10)^(precision === :quad ? -28 : -12) : T(rtol)
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("rtol must be finite and positive"))
    margin=precision === :quad ? T(2)^(-108) : 16eps(T)
    !stationary && n == m && return T[]
    evaluate(x) = getproperty(smn(m,n,c,x;spheroid,precision,normalize=true,scaled=true),stationary ? :derivative : :value).mantissa
    initial=max(32,4(n-m+1),ceil(Int,4abs(c)))
    max_points >= 2initial || throw(ArgumentError("max_points must allow two search grids (at least $(2initial))"))
    return _coordinate_zeros(evaluate,-one(T)+margin,one(T)-margin,tolerance,initial,max_points;
                             expected=stationary ? nothing : n-m,angular=true)
end

"""
    radial_zeros(m, n, c::Real, (a,b); kind=1, stationary=false,
                 spheroid=:prolate, precision=:double, rtol=nothing, max_points=16384)

Locate zeros of real radial kinds 1 or 2 (or of their coordinate derivatives)
in a finite closed interval. Require 1<a<b for prolate, 0<=a<b for oblate.
Two refined grids must agree before returning roots; bisection refines each
sign change. Grid stability is a numerical check, not a proof of completeness.
Use a larger `max_points` for long intervals or highly oscillatory solutions.
"""
function radial_zeros(m::Integer,n::Integer,c::Real,interval::Tuple{<:Real,<:Real};
        kind::Integer=1,stationary::Bool=false,spheroid::Symbol=:prolate,
        precision::Symbol=:double,rtol=nothing,max_points::Integer=16384)
    kind in (1,2) || throw(ArgumentError("real zero searches support radial kinds 1 and 2"))
    _validate_wave_arguments(m,n,c,collect(interval),spheroid,precision,:radial;kind)
    T=precision === :quad ? BigFloat : Float64
    a,b=T.(interval)
    a < b && (spheroid !== :prolate || a > 1) || throw(ArgumentError("require a nonsingular interval with a < b"))
    tolerance=rtol === nothing ? T(10)^(precision === :quad ? -28 : -12) : T(rtol)
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("rtol must be finite and positive"))
    evaluate(x) = real.(getproperty(rmn(m,n,c,x;spheroid,precision,kind,scaled=true),stationary ? :derivative : :value).mantissa)
    initial=max(32,4(n-m+1),ceil(Int,16abs(c)*(b-a)/T(pi)))
    max_points >= 2initial || throw(ArgumentError("max_points must allow two search grids (at least $(2initial))"))
    return _coordinate_zeros(evaluate,a,b,tolerance,initial,max_points)
end
