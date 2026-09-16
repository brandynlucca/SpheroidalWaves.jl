# Differentiate DLMF 30.11.3, including the joining sum in its denominator.
# No parameter differences or native radial values enter this calculation.
function _radial_spherical_bessel(top,z,kind)
    b = zeros(kind>=3 ? typeof(complex(z)) : typeof(z),top+2)
    if kind >= 3
        # DLMF 10.49.6-7 at degrees zero and one. Preserve the outgoing or
        # incoming exponential directly, including when it is very small.
        phase = kind==3 ? im : -im
        wave = exp(phase*z)/z
        b[1],b[2] = -phase*wave,-wave*(1+phase/z)
        for l in 1:top
            b[l+2] = (2l+1)/z*b[l+1]-b[l]
        end
    elseif kind == 2
        b[1],b[2] = -cos(z)/z,-cos(z)/z^2-sin(z)/z
        for l in 1:top
            b[l+2] = (2l+1)/z*b[l+1]-b[l]
        end
    elseif abs(z)>top+32
        b[1],b[2] = sin(z)/z,(sin(z)/z-cos(z))/z
        for l in 1:top
            b[l+2] = (2l+1)/z*b[l+1]-b[l]
        end
    else
        # Miller recurrence selects the recessive solution at large degree.
        # Normalize at whichever of j0,j1 has the larger magnitude.
        start = top+32+ceil(Int,abs(z))+cld(Base.precision(BigFloat),4)
        next,current = zero(z),one(z)
        for l in start:-1:1
            previous = (2l+1)/z*current-next
            l-1<=top+1 && (b[l]=previous)
            next,current = current,previous
        end
        j0,j1 = sin(z)/z,(sin(z)/z-cos(z))/z
        scale = abs(j0)>=abs(j1) ? j0/b[1] : j1/b[2]
        b .*= scale
    end
    return b
end

function _radial_expansion_data(plan,x,kind)
    m,n,c = plan.m,plan.n,plan.c
    sigma = plan.spheroid === :prolate ? 1 : -1
    z = c*x
    b = _radial_spherical_bessel(last(plan.degrees),z,kind)
    factors = [sqrt(prod(BigFloat(k) for k in l-m+1:l+m;init=big"1")*(2l+1)/2) for l in plan.degrees]
    weights,dweights = factors.*plan.v,factors.*plan.dv
    denominator,ddenominator = sum(weights),sum(dweights)
    abs(denominator)>eps(BigFloat)^(1//2)*sum(abs,weights) ||
        error("Radial joining sum is singular or numerically unresolved")
    sums,tail = zeros(eltype(b),4),zeros(BigFloat,4)
    for (i,l) in enumerate(plan.degrees)
        w,dw = (-1)^((l-n)÷2)*weights[i],(-1)^((l-n)÷2)*dweights[i]
        v = b[l+1]
        d = l/z*v-b[l+2]
        dd = (l*(l+1)/z^2-1)*v-2d/z
        terms = (w*v,w*c*d,dw*v+w*x*d,dw*c*d+w*(d+c*x*dd))
        for k in 1:4
            sums[k] += terms[k]
            i>length(weights)-8 && (tail[k]+=abs(terms[k]))
        end
    end
    # A fixed absolute floor would hide truncation error in decaying waves.
    relative_tail = maximum(iszero(s) ? (iszero(t) ? zero(t) : oftype(t,Inf)) :
                            t/abs(s) for (t,s) in zip(tail,sums))
    f,df = sums[1]/denominator,sums[2]/denominator
    tangent = (sums[3]-ddenominator*f)/denominator
    dtangent = (sums[4]-ddenominator*df)/denominator
    if sigma==1 && x==1 && m>0
        slope(v) = m==1 ? _directed_infinity(v,c isa Real,BigFloat) : m==2 ? 2v : zero(v)
        return (;state=(zero(f),slope(f),zero(tangent),slope(tangent)),relative_tail)
    end
    factor = (1-sigma/x^2)^(m//2)
    slope = sigma*m/(x*(x^2-sigma))
    # For m=0 at the regular prolate boundary this is exactly zero, not 0/0.
    m==0 && (slope=zero(x))
    return (;state=(factor*f,factor*(df+slope*f),factor*tangent,factor*(dtangent+slope*tangent)),relative_tail)
end

# Polynomial form of the radial equation and its parameter derivative.
function _radial_sensitivity_step(plan,x,state,h)
    y,dy,z,dz = state
    sigma = plan.spheroid === :prolate ? 1 : -1
    u,q,dq = x^2-sigma,plan.c^2,2plan.c
    lambda,dlambda = plan.lambda,plan.dlambda
    a = (u^2,4x*u,6x^2-2sigma,4x,one(x))
    b = (2x*u,6x^2-2sigma,6x,big"2")
    d = ((q*x^2-lambda)*u-sigma*plan.m^2,4q*x^3-2*(sigma*q+lambda)*x,
         6q*x^2-sigma*q-lambda,4q*x,q)
    dc = ((dq*x^2-dlambda)*u,4dq*x^3-2*(sigma*dq+dlambda)*x,
          6dq*x^2-sigma*dq-dlambda,4dq*x,dq)
    values,tangents = [y,dy],[z,dz]
    for k in 0:254
        rhs,drhs = zero(y),zero(z)
        for j in 1:min(4,k)
            factor = a[j+1]*(k-j+2)*(k-j+1)
            rhs += factor*values[k-j+3]
            drhs += factor*tangents[k-j+3]
        end
        for j in 0:min(3,k)
            factor = b[j+1]*(k-j+1)
            rhs += factor*values[k-j+2]
            drhs += factor*tangents[k-j+2]
        end
        for j in 0:min(4,k)
            rhs += d[j+1]*values[k-j+1]
            drhs += d[j+1]*tangents[k-j+1]+dc[j+1]*values[k-j+1]
        end
        push!(values,-rhs/(a[1]*(k+1)*(k+2)))
        push!(tangents,-drhs/(a[1]*(k+1)*(k+2)))
        if k>=30 && k%8==6
            order = length(values)-1
            output = map((values,tangents)) do coefficients
                v,dv = last(coefficients),zero(y)
                for j in order-1:-1:0
                    dv = dv*h+v
                    v = v*h+coefficients[j+1]
                end
                tail = sum(abs(coefficients[j+1]*h^j) for j in order-7:order)
                dtail = sum(abs(j*coefficients[j+1]*h^(j-1)) for j in order-7:order)
                good = tail<=big"1e-42"*max(abs(v),abs(coefficients[1]),abs(h*coefficients[2])) &&
                       dtail<=big"1e-42"*max(abs(dv),abs(coefficients[2]),abs(coefficients[1]/h))
                (;v,dv,good)
            end
            all(r->r.good,output) && return (output[1].v,output[1].dv,output[2].v,output[2].dv)
        end
    end
    error("Differentiated radial equation did not converge")
end

function _radial_expansion_batch(plan,points,kind)
    V = kind>=3 ? typeof(complex(plan.c)) : typeof(plan.c)
    states = Vector{NTuple{4,V}}(undef,length(points))
    tail = big"0"
    propagate = kind>=2 || plan.spheroid===:oblate
    interior = findall(x->x<2 && !(plan.spheroid===:prolate && x==1 && kind>=2),points)
    if propagate && !isempty(interior)
        initial = _radial_expansion_data(plan,big"2",kind)
        tail = max(tail,initial.relative_tail)
        x,state = big"2",initial.state
        for i in sort(interior;by=i->points[i],rev=true)
            target = BigFloat(points[i])
            steps = 0
            while x>target
                radius = plan.spheroid===:prolate ? x-1 : sqrt(x^2+1)
                h = -min(x-target,radius/4,inv(1+abs(plan.c)+sqrt(abs(plan.lambda))))
                x+h<x || error("Radial sensitivity propagation exhausted coordinate precision")
                state = _radial_sensitivity_step(plan,x,state,h)
                x += h
                steps += 1
                steps<=100000 || error("Radial sensitivity propagation exceeded its step limit")
            end
            if plan.spheroid===:oblate && target==0 && kind==1
                # The regular solution has exact parity at the oblate origin.
                states[i] = isodd(plan.n-plan.m) ? (zero(state[1]),state[2],zero(state[3]),state[4]) :
                                                   (state[1],zero(state[2]),state[3],zero(state[4]))
            else
                states[i] = state
            end
        end
    end
    for (i,x) in enumerate(points)
        if plan.spheroid===:prolate && x==1 && kind>=2
            states[i] = ntuple(_->V(NaN),4)
        elseif !(propagate && i in interior)
            result = _radial_expansion_data(plan,BigFloat(x),kind)
            states[i] = result.state
            tail = max(tail,result.relative_tail)
        end
    end
    return (;states,tail,propagated=propagate && !isempty(interior))
end

_radial_needs_analytic(c,points,kind) = c isa Complex &&
    (iszero(real(c)) || kind>=3)

function _radial_analytic_data(m,n,c,points,spheroid,precision,kind;eigenvalue_seed=nothing)
    spheroid===:prolate && c isa Complex && any(isone,points) &&
        throw(DomainError(points,"complex prolate radial evaluation requires x > 1"))
    T = precision===:quad ? BigFloat : Float64
    parameter = c isa Real ? T(c) : Complex{T}(c)
    coordinates = BigFloat.(points)
    # Direct Hankel waves do not lose exp(2*abs(imag(c*x))) through an R1 +/- iR2
    # subtraction. Keep guards for the coefficient solve and degree recurrence.
    cancellation = kind>=3 ? zero(real(parameter)) :
                   4abs(imag(parameter))*max(2,maximum(coordinates))
    bits = max(Base.precision(BigFloat)+64,256+2*(n+m)+ceil(Int,8abs(parameter)+cancellation))
    result,plan = setprecision(BigFloat,bits) do
        # Neumann terms converge geometrically at x=2, even when the angular
        # coefficients themselves have already become very small.
        minimum = max(kind==1 ? 32 : 96,(n-m)÷2+16+ceil(Int,abs(parameter)))
        for attempt in 1:4
            max_terms = max(512,2minimum)
            rtol = min(big"1e-42",exp(-2BigFloat(abs(parameter)))*big"1e-42")
            plan = _coefficient_plan(m,n,parameter;spheroid,precision,rtol,min_terms=minimum,max_terms,eigenvalue_seed)
            plan.converged || error("Radial sensitivity coefficients did not converge")
            result = _radial_expansion_batch(plan,coordinates,kind)
            result.tail<=big"1e-34" && return result,plan
            minimum *= 2
        end
        error("Differentiated radial expansion did not converge")
    end
    return result,plan
end

function _radial_analytic_values(m,n,c,points,spheroid,precision,kind;second_derivative=false,eigenvalue_seed=nothing)
    result,plan = _radial_analytic_data(m,n,c,points,spheroid,precision,kind;eigenvalue_seed)
    values = (;value=[s[1] for s in result.states],derivative=[s[2] for s in result.states])
    second_derivative || return values
    seconds = setprecision(BigFloat,Base.precision(real(plan.c))) do
        [begin
            a,b,d=_wave_equation(m,plan.c,BigFloat(x),plan.lambda,spheroid,:radial)
            -(b*s[2]+d*s[1])/a
         end for (x,s) in zip(points,result.states)]
    end
    return (;values...,second_derivative=seconds)
end

function _radial_analytic_wave(m,n,c,points,spheroid,precision,kind,scaled,logderivative,second_derivative)
    result=_radial_analytic_values(m,n,c,points,spheroid,precision,kind;second_derivative)
    R=precision===:quad ? BigFloat : Float64
    output=map(v->scaled ? _decimal_scaled(v,Complex{R}) : complex.(R.(real.(v)),R.(imag.(v))),result)
    if logderivative
        ratios=[iszero(v) ? complex(R(NaN),zero(R)) : complex(R(real(d/v)),R(imag(d/v))) for (v,d) in zip(result.value,result.derivative)]
        output=(;output...,logderivative=ratios)
    end
    return output
end

function _radial_analytic_jacobian(m,n,c,points,spheroid,precision,kind,with_metadata)
    result,plan=_radial_analytic_data(m,n,c,points,spheroid,precision,kind)
    T=precision===:quad ? BigFloat : Float64
    # Match rmn's complex output type even for real first/second-kind results.
    rounded(z) = complex(T(real(z)),T(imag(z)))
    value = [rounded(s[3]) for s in result.states]
    derivative = [rounded(s[4]) for s in result.states]
    metadata(v) = (;_coefficient_derivative_metadata(plan,v)...,
        method=result.propagated ? :differentiated_equation : :differentiated_expansion,
        radial_series_tail=result.tail)
    mv,md = metadata(value),metadata(derivative)
    if c isa Real
        output = (;dvalue_dc=value,dderivative_dc=derivative)
        return with_metadata ? (;output...,metadata_value=mv,metadata_derivative=md) : output
    end
    output = (;dvalue_dcreal=value,dvalue_dcimag=im.*value,
               dderivative_dcreal=derivative,dderivative_dcimag=im.*derivative)
    return with_metadata ? (;output...,metadata_value_dcreal=mv,metadata_value_dcimag=mv,
                            metadata_derivative_dcreal=md,metadata_derivative_dcimag=md) : output
end
