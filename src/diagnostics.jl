# Coefficients of A*y'' + B*y' + C*y = 0, in the package's lambda convention.
function _wave_equation(m, c, x, lambda, spheroid, target)
    sigma = spheroid === :prolate ? 1 : -1
    if target === :angular
        a = (1-x)*(1+x)
        return a, -2x, lambda-sigma*c^2*x^2-m^2/a
    end
    a = spheroid === :prolate ? (x-1)*(x+1) : x^2+1
    return a, 2x, c^2*x^2-lambda-sigma*m^2/a
end

function _wave_second_derivative(result, m, n, c, points, spheroid, precision, target; option=0)
    lambda = eigenvalue(m,n,c;spheroid,precision)
    T = precision === :quad ? BigFloat : Float64
    second = similar(result.value)
    for (i, point) in enumerate(points)
        x = eltype(result.value) <: Union{BigFloat,Complex{BigFloat}} ? BigFloat(point) : T(point)
        singular = target === :angular ? abs(x) == 1 : spheroid === :prolate && x == 1
        if singular
            second[i] = _wave_endpoint_second(result.value[i],result.derivative[i],m,n,c,x,
                                               lambda,spheroid,precision,target,option)
        elseif (target === :angular && 1-abs(x) < 1//65536) ||
               (target === :radial && spheroid === :prolate && option == 1 && x-1 < 1//65536)
            second[i] = _regular_near_endpoint_second(m,n,c,BigFloat(point),lambda,spheroid,precision,target,option)
        else
            a,b,d = _wave_equation(m,c,x,lambda,spheroid,target)
            second[i] = -(b*result.derivative[i]+d*result.value[i])/a
        end
    end
    return (;result...,second_derivative=second)
end

function _wave_endpoint_second(value,derivative,m,n,c,x,lambda,spheroid,precision,target,option)
    T = precision === :quad ? BigFloat : Float64
    target === :radial && option != 1 && return T(NaN)
    m > 4 && return zero(value)
    amplitude,polynomial = _regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option)
    result = if m == 0
        2get(polynomial,3,zero(amplitude))*amplitude
    elseif m == 2
        (target === :angular ? 1 : -1)*2amplitude*(2get(polynomial,2,zero(amplitude))-1)
    elseif m == 4
        8amplitude
    else
        return _directed_infinity((m == 1 ? -1 : 1)*amplitude,c isa Real,T)
    end
    return c isa Real ? real(result) : result
end

function _directed_infinity(direction,real_result,T)
    component(z) = iszero(z) ? zero(T) : copysign(T(Inf),z)
    return real_result ? component(real(direction)) : complex(component(real(direction)),component(imag(direction)))
end

function _regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option; radius=0)
    # Recover the regular factor at the endpoint from a nonsingular binary
    # anchor. The same series applies radially with t=1-x < 0.
    q = (spheroid === :prolate ? 1 : -1)*c^2
    L = lambda isa Real ? BigFloat(lambda) : Complex{BigFloat}(lambda)
    Q = q isa Real ? BigFloat(q) : Complex{BigFloat}(q)
    distance = BigFloat(2)^(-ceil(Int,log2(max(big"256",32*(1+abs(L)+abs(Q)+m*(m+1))))))
    point = target === :angular ? sign(x)*(1-distance) : 1+distance
    r = _scaled_native_values(m,n,c,[point],spheroid,precision,target,option)
    polynomial = _angular_endpoint_coefficients(m,L,Q,max(distance,radius))
    t = target === :angular ? distance : -distance
    shape = _angular_endpoint_polynomial(polynomial,t).value
    factor = distance*(target === :angular ? 2-distance : 2+distance)
    amplitude = only(r.value)/(factor^(m//2)*shape)
    return amplitude,polynomial
end

function _regular_near_endpoint_second(m,n,c,x,lambda,spheroid,precision,target,option)
    amplitude,polynomial = _regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option;radius=abs(1-abs(x)))
    coordinate = abs(x)
    t = 1-coordinate
    shape = _angular_endpoint_polynomial(polynomial,t)
    u = target === :angular ? (1-coordinate)*(1+coordinate) : (coordinate-1)*(coordinate+1)
    sigma = target === :angular ? 1 : -1
    result = amplitude*u^(m//2)*(shape.second_derivative+2sigma*m*coordinate/u*shape.derivative+
                    (m*(m-2)*coordinate^2/u^2-sigma*m/u)*shape.value)
    return c isa Real ? real(result) : result
end

function _fix_regular_endpoints!(result,m,n,c,points,spheroid,precision,target,option)
    iszero(c) && return result # Existing analytic spherical values already use limits.
    target === :radial && (spheroid !== :prolate || option != 1) && return result
    indices = target === :angular ? findall(x -> abs(x)==1,points) : findall(isone,points)
    isempty(indices) && return result
    lambda=eigenvalue(m,n,c;spheroid,precision)
    T=precision === :quad ? BigFloat : Float64
    for i in indices
        x=points[i]
        if m > 2
            result.value[i]=0
            result.derivative[i]=0
            continue
        end
        amplitude,p=_regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option)
        c isa Real && (amplitude=real(amplitude))
        result.value[i]=m == 0 ? amplitude : zero(amplitude)
        direction=target === :angular ? -sign(x) : 1
        result.derivative[i] = if m == 0
            -get(p,2,zero(amplitude))*amplitude*(target === :angular ? sign(x) : 1)
        elseif m == 1
            _directed_infinity(direction*amplitude,c isa Real,T)
        else
            2direction*amplitude
        end
    end
    return result
end

# Exact polynomial differentiation weights. They are independent of the ODE
# and of the native coordinate derivatives being checked.
function _difference_weights(offsets)
    first, second = Rational{BigInt}[], Rational{BigInt}[]
    for j in eachindex(offsets)
        polynomial = Rational{BigInt}[1]
        denominator = big(1)
        for k in eachindex(offsets)
            k == j && continue
            next = zeros(Rational{BigInt},length(polynomial)+1)
            next[1:end-1] .-= offsets[k].*polynomial
            next[2:end] .+= polynomial
            polynomial = next
            denominator *= offsets[j]-offsets[k]
        end
        push!(first,polynomial[2]/denominator)
        push!(second,2polynomial[3]/denominator)
    end
    return first,second
end

# Internal validation helper; not part of the public API.
#     _spheroidal_residual(m, n, c, x; target=:angular, spheroid=:prolate,
#                         precision=:double, kind=1, normalize=false, h=nothing)
#
# Check the differential equation using first and second coordinate derivatives
# estimated independently from nine function values. Return vectors `residual`,
# `relative_residual` (divided by the sum of absolute ODE terms), `derivative`,
# `second_derivative`, `derivative_change`, `second_derivative_change`, and `step`.
# The changes compare steps h and h/2; the reported derivatives use h/2.
# These are consistency diagnostics, not accuracy bounds. No ODE-derived second
# derivative enters the check. Endpoints singular in the ODE are rejected; oblate
# x=0 uses a forward stencil. Scalar coordinates return length-one vectors.
function _spheroidal_residual(m::Integer,n::Integer,c::Union{Real,Complex},points::AbstractVector{<:Real};
        target::Symbol=:angular,spheroid::Symbol=:prolate,precision::Symbol=:double,
        kind::Integer=1,normalize::Bool=false,h=nothing)
    target in (:angular,:radial) || throw(ArgumentError("target must be :angular or :radial"))
    _validate_wave_arguments(m,n,c,points,spheroid,precision,target;kind)
    T = precision === :quad ? BigFloat : Float64
    xs = T.(points)
    singular = target === :angular ? any(x -> abs(x)==1,xs) : spheroid === :prolate && any(isone,xs)
    singular && throw(DomainError(points,"residual stencils require nonsingular coordinates"))
    h !== nothing && !(isfinite(h) && h > 0) && throw(ArgumentError("h must be finite and positive"))
    evaluate(z) = target === :angular ? smn(m,n,c,z;spheroid,precision,normalize,kind) :
                                       rmn(m,n,c,z;spheroid,precision,kind)
    lambda = eigenvalue(m,n,c;spheroid,precision)
    center = evaluate(xs)
    residual, first, second = similar(center.value),similar(center.value),similar(center.value)
    relative, first_change, second_change, steps = (zeros(T,length(xs)) for _ in 1:4)
    for (i,x) in enumerate(xs)
        epsilon = precision === :quad ? T(2)^(-112) : eps(T)
        step = h === nothing ? epsilon^(T(1)/10)*max(one(T),abs(x))/4 : T(h)
        if target === :angular
            step = min(step,(1-abs(x))/32)
        elseif spheroid === :prolate
            step = min(step,(x-1)/32)
        end
        offsets = target === :radial && spheroid === :oblate && x < 4step ? collect(0:8) : collect(-4:4)
        w1,w2 = _difference_weights(offsets)
        function estimate(s)
            grid = x .+ s.*offsets
            length(unique(grid)) == 9 || throw(ArgumentError("finite difference step is below coordinate resolution"))
            values = evaluate(grid).value .- center.value[i]
            return sum(T.(w1).*values)/s,sum(T.(w2).*values)/s^2
        end
        coarse1,coarse2 = estimate(step)
        first[i],second[i] = estimate(step/2)
        first_change[i],second_change[i] = abs(first[i]-coarse1),abs(second[i]-coarse2)
        a,b,d = _wave_equation(m,c,x,lambda,spheroid,target)
        terms = (a*second[i],b*first[i],d*center.value[i])
        residual[i] = sum(terms)
        scale = sum(abs,terms)
        relative[i] = iszero(scale) ? zero(T) : abs(residual[i])/scale
        steps[i] = step/2
    end
    return (;residual,relative_residual=relative,derivative=first,second_derivative=second,
            derivative_change=first_change,second_derivative_change=second_change,step=steps)
end

_spheroidal_residual(m::Integer,n::Integer,c::Union{Real,Complex},x::Real;kwargs...) =
    _spheroidal_residual(m,n,c,[x];kwargs...)
