# Ferrers (on-the-cut) Qs on -1 < x < 1. DLMF 30.5.2 and 30.5.4 fix
# its origin data. The joining sums include the finite extension of the
# Legendre coefficients below degree m: DLMF 30.8(ii), 30.11.4.
function _qs_initial_data(m,n,c,spheroid,precision;mode=nothing)
    rtol = min(big"1e-40",exp(-2BigFloat(abs(c)))*big"1e-40")
    max_terms = max(512,(n-m)÷2+2ceil(Int,abs(c))+32)
    plan = _coefficient_plan(m,n,c;spheroid,precision,rtol,max_terms,
                             eigenvalue_seed=mode === nothing ? nothing : mode.lambda)
    plan.converged || error("Qs coefficient expansion did not converge")
    scale = _coefficient_phase(plan,precision;mode)*sqrt(_ferrers_norm2(m,n,BigFloat))
    d = [scale*v/sqrt(_ferrers_norm2(m,l,BigFloat)) for (v,l) in zip(plan.v,plan.degrees)]
    dc = [scale*v/sqrt(_ferrers_norm2(m,l,BigFloat)) for (v,l) in zip(plan.dv,plan.degrees)]
    g = (spheroid === :prolate ? 1 : -1)*plan.c^2
    dg = 2*(spheroid === :prolate ? 1 : -1)*plan.c
    # Weighted sum = (n+m)!/(n-m)! * A^(-m).
    weights = [prod(BigFloat(k) for k in l-m+1:l+m;init=big"1") for l in plan.degrees]
    weighted = sum(d.*weights)
    joining = sum(d)
    dweighted,djoining = sum(dc.*weights),sum(dc)
    if m > 0 && !iszero(g)
        # Work with d_l=(-1)^k a_{n,k}^m, so off-diagonal signs reverse.
        previous,current = zero(g),one(g)
        dprevious,dcurrent = zero(g),zero(g)
        lower = typeof(g)[]
        dlower = typeof(g)[]
        for l in -m+mod(n-m,2):2:first(plan.degrees)-2
            push!(lower,current)
            push!(dlower,dcurrent)
            a = g*(l-m-1)*(l-m)/((2l-3)*(2l-1))
            shift = g*(1-2BigFloat(l*(l+1)-1+m^2)/((2l-1)*(2l+3)))
            b = l*(l+1)+shift-plan.lambda
            cc = g*(l+m+1)*(l+m+2)/((2l+3)*(2l+5))
            numerator = a*previous+b*current
            bound = abs(a*previous)+(abs(l*(l+1))+abs(shift)+abs(plan.lambda))*abs(current)
            abs(numerator) > sqrt(eps(BigFloat))*bound ||
                throw(DomainError(c,"Qs finite coefficient extension is singular or numerically unresolved (DLMF 30.8(ii))"))
            # Differentiate the finite extension, including its normalization.
            dnumerator = (dg/g*a)*previous+a*dprevious+
                         (dg/g*shift-plan.dlambda)*current+b*dcurrent
            dprevious,dcurrent = dcurrent,-dnumerator/cc+numerator*(dg/g)/cc
            previous,current = current,-numerator/cc
        end
        iszero(current) && throw(DomainError(c,"Qs is undefined: its finite coefficient extension is singular (DLMF 30.8(ii))"))
        joining += sum(lower)*(first(d)/current)
        djoining += sum(dlower)*(first(d)/current)+
                    sum(lower)*(first(dc)/current-first(d)*dcurrent/current^2)
    end
    # For c=0 the extension vanishes and the factorial Wronskian is exact.
    constant = joining*weighted
    dconstant = djoining*weighted+joining*dweighted
    isfinite(constant) && !iszero(constant) ||
        error("Qs normalization is singular or numerically unresolved")
    tail = max(sum(abs,d[max(1,end-3):end])/abs(joining),
               sum(abs,(d.*weights)[max(1,end-3):end])/abs(weighted))
    tail <= big"1e-34" || error("Qs joining sums did not converge to the required precision")
    origin = _evaluate_coefficient_vector(plan,scale.*plan.v,[big"0"])
    dorigin = _evaluate_coefficient_vector(plan,scale.*plan.dv,[big"0"])
    if iseven(n-m)
        y,dy = zero(constant),constant/only(origin.value)
        z,dz = zero(constant),(dconstant-dy*only(dorigin.value))/only(origin.value)
    else
        y,dy = -constant/only(origin.derivative),zero(constant)
        z,dz = (-dconstant-y*only(dorigin.derivative))/only(origin.derivative),zero(constant)
    end
    isfinite(y) && isfinite(dy) || error("Qs origin normalization is numerically unresolved")
    # P ~ a*(1-x^2)^(m/2) at +1. This also fixes the direction of Q's
    # logarithmic (m=0) or algebraic (m>0) endpoint divergence.
    a = (-1)^m*weighted/(BigFloat(2)^m*factorial(big(m)))
    amplitude = constant/(m == 0 ? a : 2m*a)
    return (;lambda=plan.lambda,dlambda=plan.dlambda,g,dg,y,dy,z,dz,amplitude,constant,dconstant,plan)
end

# Taylor propagation of the polynomial-coefficient equation
# (1-x^2)^2 y'' - 2x(1-x^2)y' + ((lambda-g*x^2)(1-x^2)-m^2)y=0.
# Steps remain inside one quarter of the distance to its nearest singularity.
# A tail check controls both the value and its derivative at every step.
function _qs_step(m,lambda,g,x,y,dy,h;sensitivity=nothing)
    u = (1-x)*(1+x)
    a = (u^2,-4x*u,6x^2-2,4x,one(x))
    b = (-2x*u,6x^2-2,6x,big"2")
    d = ((lambda-g*x^2)*u-m^2,-2x*lambda-2g*x+4g*x^3,-lambda-g+6g*x^2,4g*x,g)
    coefficients = [y,dy]
    if sensitivity !== nothing
        z,dz,dlambda,dg = sensitivity
        tangents = [z,dz]
        dc = ((dlambda-dg*x^2)*u,-2x*dlambda-2dg*x+4dg*x^3,
              -dlambda-dg+6dg*x^2,4dg*x,dg)
    end
    tolerance = big"1e-45"
    for k in 0:254
        total = zero(y)
        for j in 1:min(4,k)
            total += a[j+1]*(k-j+2)*(k-j+1)*coefficients[k-j+3]
        end
        for j in 0:min(3,k)
            total += b[j+1]*(k-j+1)*coefficients[k-j+2]
        end
        for j in 0:min(4,k)
            total += d[j+1]*coefficients[k-j+1]
        end
        push!(coefficients,-total/(a[1]*(k+1)*(k+2)))
        if sensitivity !== nothing
            tangent = zero(y)
            for j in 1:min(4,k)
                tangent += a[j+1]*(k-j+2)*(k-j+1)*tangents[k-j+3]
            end
            for j in 0:min(3,k)
                tangent += b[j+1]*(k-j+1)*tangents[k-j+2]
            end
            for j in 0:min(4,k)
                tangent += d[j+1]*tangents[k-j+1]+dc[j+1]*coefficients[k-j+1]
            end
            push!(tangents,-tangent/(a[1]*(k+1)*(k+2)))
        end
        if k >= 30 && k % 8 == 6
            order = length(coefficients)-1
            value,derivative = coefficients[end],zero(y)
            for j in order-1:-1:0
                derivative = derivative*h+value
                value = value*h+coefficients[j+1]
            end
            tail = sum(abs(coefficients[j+1]*h^j) for j in order-7:order)
            dtail = sum(abs(j*coefficients[j+1]*h^(j-1)) for j in order-7:order)
            if tail <= tolerance*max(abs(value),abs(y),abs(h*dy)) &&
               dtail <= tolerance*max(abs(derivative),abs(dy),abs(y/h))
                sensitivity === nothing && return value,derivative
                v,dv = tangents[end],zero(y)
                for j in order-1:-1:0
                    dv = dv*h+v
                    v = v*h+tangents[j+1]
                end
                ztail = sum(abs(tangents[j+1]*h^j) for j in order-7:order)
                dztail = sum(abs(j*tangents[j+1]*h^(j-1)) for j in order-7:order)
                if ztail <= tolerance*max(abs(v),abs(z),abs(h*dz)) &&
                   dztail <= tolerance*max(abs(dv),abs(dz),abs(z/h))
                    return value,derivative,v,dv
                end
            end
        end
    end
    error("Qs Taylor propagation did not converge")
end

function _qs_values(m,n,c,points,spheroid,precision;sensitivity=false,mode=nothing)
    data = _qs_initial_data(m,n,c,spheroid,precision;mode)
    values = fill(zero(data.y),length(points))
    derivatives,seconds = similar(values),similar(values)
    tangents,dtangents = similar(values),similar(values)
    x,y,dy = big"0",data.y,data.dy
    z,dz = data.z,data.dz
    step_limit = inv(1+sqrt(abs(data.lambda))+abs(c)+m)
    for i in sortperm(abs.(points))
        target = abs(BigFloat(points[i]))
        parity = points[i] < 0 && iseven(n-m) ? -1 : 1
        dparity = points[i] < 0 ? -parity : parity
        if target == 1
            divergent = _directed_infinity(data.amplitude,c isa Real,BigFloat)
            values[i],derivatives[i],seconds[i] = parity*divergent,dparity*divergent,parity*divergent
            # Qs is singular at the endpoint itself. Do not infer a parameter
            # derivative by subtracting infinities or just its leading amplitude.
            tangents[i] = dtangents[i] = oftype(data.y,NaN)
            continue
        end
        steps = 0
        while x < target
            h = min(target-x,(1-x)/4,step_limit)
            x+h > x || error("Qs coordinate propagation exhausted working precision")
            if sensitivity
                y,dy,z,dz = _qs_step(m,data.lambda,data.g,x,y,dy,h;
                                     sensitivity=(z,dz,data.dlambda,data.dg))
            else
                y,dy = _qs_step(m,data.lambda,data.g,x,y,dy,h)
            end
            x += h
            steps += 1
            steps <= 100000 || error("Qs coordinate propagation exceeded its step limit")
        end
        u = (1-target)*(1+target)
        values[i],derivatives[i] = parity*y,dparity*dy
        seconds[i] = parity*(2target*dy-(data.lambda-data.g*target^2-m^2/u)*y)/u
        tangents[i],dtangents[i] = parity*z,dparity*dz
    end
    result = (;value=values,derivative=derivatives,second_derivative=seconds)
    return sensitivity ? (;result...,dvalue_dc=tangents,dderivative_dc=dtangents,plan=data.plan) : result
end

_qs_working_bits(c) = max(256,Base.precision(BigFloat)+64,256+ceil(Int,5abs(c)))

function _angular_second_kind(m,n,c,points,spheroid,precision,normalize,scaled,logderivative,second_derivative;mode=nothing)
    # Unit-integral normalization is not defined for the singular family.
    normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1; Qs uses the DLMF second-kind normalization"))
    R = precision === :quad ? BigFloat : Float64
    T = c isa Real ? R : Complex{R}
    parameter = c isa Real ? R(c) : Complex{R}(c)
    coordinates = BigFloat.(points)
    # Extra guard digits protect coefficient sums and propagation, including
    # oblate solutions concentrated near an endpoint. Quad inputs never pass
    # through Float64. This is still a double/quad API, not arbitrary precision.
    bits = _qs_working_bits(parameter)
    result = setprecision(BigFloat,bits) do
        _qs_values(m,n,parameter,coordinates,spheroid,precision;mode)
    end
    selected = second_derivative ? result : (;result.value,result.derivative)
    # Complex{BigFloat}(z) is an identity conversion when z already has that
    # type. Round its components explicitly to avoid leaking guard precision.
    output = map(selected) do v
        scaled ? _decimal_scaled(v,T) :
            c isa Real ? R.(v) : complex.(R.(real.(v)),R.(imag.(v)))
    end
    if logderivative
        ratio = [abs(x)==1 ? T(copysign(Inf,x)) : iszero(v) ? T(NaN) : T(d/v)
                 for (x,v,d) in zip(coordinates,result.value,result.derivative)]
        output = (;output...,logderivative=ratio)
    end
    return output
end
