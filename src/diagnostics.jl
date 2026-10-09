# Complex coordinates continue the real-domain normalization on fixed cuts.
function _validate_complex_coordinates(m,n,c,points,spheroid,precision,target;kind=1,normalization=:standard)
    _validate_wave_arguments(m,n,c,[target===:angular ? 0 : 2],spheroid,precision,target;kind,normalization)
    isempty(points) && throw(ArgumentError("evaluation points must not be empty"))
    all(isfinite,points) || throw(ArgumentError("evaluation points must be finite"))
    valid(z) = target===:angular ? (!isreal(z) || abs(real(z))<=1) : spheroid===:prolate ?
        (!isreal(z) || real(z)>=1 || real(z)<-1) : (!iszero(real(z)) || abs(imag(z))<1)
    all(valid,points) || throw(DomainError(points,"coordinates lie on a branch cut or an unsupported singularity"))
    return nothing
end

_coordinate_is_endpoint(z,target,spheroid) = isreal(z) &&
    (target===:angular ? abs(real(z))==1 : real(z)==(spheroid===:prolate ? 1 : 0))

function _coordinate_precision(m,n,c,points)
    parameter = Complex{_SWFloat}(c)
    radius = maximum(z -> abs(Complex{_SWFloat}(z)),points)
    growth = ceil(Int,32abs(parameter)*(2+radius)+8(n+m+1)*log2(2+radius))
    max(384,_angular_precision(parameter),_static_precision(m,n,parameter))+growth
end

function _coordinate_initial_data(m,n,c,spheroid,precision,target,kind,normalize,normalization;mode=nothing)
    tolerance = sqrt(eps(_SWFloat))
    if target===:angular
        if kind==2
            data = _qs_initial_data(m,n,c,spheroid,precision;mode,rtol=tolerance)
            plan,state = data.plan,(data.y,data.dy,data.z,data.dz)
        else
            plan = _coefficient_plan(m,n,c;spheroid,precision,rtol=tolerance,
                max_terms=max(512,(n-m)÷2+2ceil(Int,abs(c))+32),eigenvalue_seed=mode===nothing ? nothing : mode.lambda)
            plan.converged || error("Complex-coordinate angular expansion did not converge")
            scale = _coefficient_phase(plan,precision;mode)*(normalize ? one(_SWFloat) : sqrt(_ferrers_norm2(m,n,_SWFloat)))
            value = _evaluate_coefficient_vector(plan,scale.*plan.v,[0])
            tangent = _evaluate_coefficient_vector(plan,scale.*plan.dv,[0])
            state = (only(value.value),only(value.derivative),only(tangent.value),only(tangent.derivative))
        end
        sigma = spheroid===:prolate ? 1 : -1
        equation = (;plan...,spheroid=:prolate,q=sigma*plan.c^2,dq=2sigma*plan.c)
        return (;plan,equation,state=complex.(state),anchor=complex(zero(_SWFloat)))
    end
    anchor = _SWFloat(2)
    if iszero(c)
        plan = _sensitivity_plan(m,n,c,spheroid,precision)
        sigma = spheroid===:prolate ? 1 : -1
        y,dy = kind==1 ? _static_regular_radial(m,n,anchor,sigma) : _static_irregular_radial(m,n,anchor,sigma)
        state = (y,dy,zero(y),zero(y))
    else
        data,plan = _radial_analytic_data(m,n,c,[anchor],spheroid,precision,kind;
            eigenvalue_seed=mode===nothing ? nothing : mode.lambda,rtol=tolerance)
        state = only(data.states)
        if normalization===:static
            power = kind==1 ? -n : n+1
            y,dy,t,dt = state
            factor = plan.c^power
            state = (factor*y,factor*dy,factor*(t+power*y/plan.c),factor*(dt+power*dy/plan.c))
        end
    end
    return (;plan,equation=plan,state=complex.(state),anchor=complex(anchor))
end

function _coordinate_path(z,target,spheroid)
    target===:angular && return [z]
    spheroid===:oblate && return real(z)<0 ? [zero(z),z] : [z]
    height = iszero(imag(z)) && real(z)<-1 ? one(real(z)) : imag(z)
    return [complex(oftype(real(z),2),height),complex(real(z),height),z]
end

function _coordinate_continuation(equation,state,anchor,path;max_steps=100000,step_fraction=1,step_count=Ref(0))
    x = anchor
    tolerance = sqrt(eps(_SWFloat))
    singularity = equation.spheroid===:prolate ? one(x) : complex(zero(real(x)),one(real(x)))
    step_limit = inv(1+abs(equation.c)+sqrt(abs(equation.lambda)))
    steps = 0
    for destination in path
        while x!=destination
            delta = destination-x
            radius = min(abs(x-singularity),abs(x+singularity))
            length = min(abs(delta),step_fraction*min(radius/4,step_limit))
            h = delta*(length/abs(delta))
            x+h!=x || error("Complex-coordinate continuation exhausted coordinate precision")
            state = _radial_sensitivity_step(equation,x,state,h;rtol=tolerance)
            x = length==abs(delta) ? destination : x+h
            steps += 1
            step_count[] += 1
            steps<=max_steps || error("Complex-coordinate continuation exceeded its step limit")
        end
    end
    return state
end

function _coordinate_batch(initial,points,target,spheroid;step_fraction=1)
    paths = [_coordinate_path(z,target,spheroid) for z in points]
    states = Vector{typeof(initial.state)}(undef,length(points))
    step_count = Ref(0)
    function continue_paths(indices,depth,anchor,state)
        groups = Dict{Tuple{Int,_SWFloat},Vector{Int}}()
        for i in indices
            delta = paths[i][depth]-anchor
            direction = iszero(real(delta)) ? (0,_SWFloat(sign(imag(delta)))) :
                (real(delta)>0 ? 1 : -1,_with_swprecision(() -> imag(delta)/abs(real(delta)),2Base.precision(_SWFloat)+16))
            push!(get!(() -> Int[],groups,direction),i)
        end
        for indices in values(groups)
            sort!(indices;by=i -> (abs(paths[i][depth]-anchor),real(paths[i][depth]),imag(paths[i][depth])))
            current,position = state,anchor
            first = 1
            while first<=length(indices)
                destination = paths[indices[first]][depth]
                last = first
                while last<length(indices) && paths[indices[last+1]][depth]==destination
                    last += 1
                end
                current = _coordinate_continuation(initial.equation,current,position,[destination];step_fraction,step_count)
                position = destination
                children = Int[]
                for j in first:last
                    i = indices[j]
                    depth==length(paths[i]) ? (states[i]=current) : push!(children,i)
                end
                isempty(children) || continue_paths(children,depth+1,position,current)
                first = last+1
            end
        end
    end
    continue_paths(eachindex(points),1,initial.anchor,initial.state)
    return (;states,steps=step_count[])
end

function _complex_coordinate_data(m,n,c,points,spheroid,precision,target,kind,normalize,normalization;mode=nothing,step_fraction=1)
    return _with_swprecision(_coordinate_precision(m,n,c,points)) do
        R = precision===:quad ? BigFloat : Float64
        parameter = c isa Real ? _input_float(R,c) : _input_float(Complex{R},c)
        initial = _coordinate_initial_data(m,n,parameter,spheroid,precision,target,kind,normalize,normalization;mode)
        interior = findall(z -> !_coordinate_is_endpoint(z,target,spheroid),points)
        batch = _coordinate_batch(initial,Complex{_SWFloat}.(points[interior]),target,spheroid;step_fraction)
        states = Vector{typeof(initial.state)}(undef,length(points))
        states[interior] = batch.states
        for (i,point) in enumerate(points)
            if _coordinate_is_endpoint(point,target,spheroid)
                options = target===:angular ? (;normalize) : (;normalization)
                f,jac = target===:angular ? (smn,jacobian_smn) : (rmn,jacobian_rmn)
                value = f(m,n,c,real(point);spheroid,precision,kind,options...)
                tangent = jac(m,n,c,[real(point)];spheroid,precision,kind,options...)
                v,d = c isa Real ? (tangent.dvalue_dc,tangent.dderivative_dc) : (tangent.dvalue_dcreal,tangent.dderivative_dcreal)
                states[i] = Complex{_SWFloat}.((only(value.value),only(value.derivative),only(v),only(d)))
            end
        end
        return (;states,initial.plan,batch.steps)
    end
end

function _complex_coordinate_wave(m,n,c,points,spheroid,precision,target,kind,normalize,
                                   scaled,logderivative,second_derivative,derivatives,normalization)
    _validate_complex_coordinates(m,n,c,points,spheroid,precision,target;kind,normalization)
    kind==2 && normalize && throw(ArgumentError("normalize=true is only defined for angular kind=1"))
    derivatives in 1:4 || throw(ArgumentError("derivatives must be in 1:4"))
    order = max(derivatives,second_derivative ? 2 : 1)
    R = precision===:quad ? BigFloat : Float64
    T = Complex{R}
    fields = (:value,:derivative,:second_derivative,:third_derivative,:fourth_derivative)[1:order+1]
    values = _with_swprecision(_coordinate_precision(m,n,c,points)) do
        data = _complex_coordinate_data(m,n,c,points,spheroid,precision,target,kind,normalize,normalization)
        arrays = [zeros(Complex{_SWFloat},length(points)) for _ in fields]
        for (i,point) in enumerate(points)
            if _coordinate_is_endpoint(point,target,spheroid)
                f = target===:angular ? smn : rmn
                options = target===:angular ? (;normalize) : (;normalization)
                endpoint = f(m,n,c,real(point);spheroid,precision,kind,derivatives=order,options...)
                for (array,field) in zip(arrays,fields)
                    array[i] = only(getproperty(endpoint,field))
                end
                continue
            end
            z = Complex{_SWFloat}(point)
            y,dy = data.states[i][1:2]
            arrays[1][i],arrays[2][i] = y,dy
            order==1 && continue
            a,b,d = _wave_equation(m,data.plan.c,z,data.plan.lambda,spheroid,target)
            ddy = -(b*dy+d*y)/a
            arrays[3][i] = ddy
            if order>=3
                sigma = spheroid===:prolate ? 1 : -1
                bprime = target===:angular ? -2 : 2
                q = (target===:angular ? -sigma : 1)*data.plan.c^2
                mass = (target===:angular ? 1 : sigma)*m^2
                dprime = 2q*z+mass*b/a^2
                dsecond = 2q+mass*(bprime/a^2-2b^2/a^3)
                third = -(2b*ddy+(bprime+d)*dy+dprime*y)/a
                arrays[4][i] = third
                order==4 && (arrays[5][i] = -(3b*third+(3bprime+d)*ddy+2dprime*dy+dsecond*y)/a)
            end
        end
        NamedTuple{fields}(Tuple(arrays))
    end
    output = map(v -> scaled ? _decimal_scaled(v,T) : T.(v),values)
    if logderivative
        ratio = T[isreal(z) ? _wave_logderivative(v,d,real(z),m,spheroid,target,kind) :
            iszero(v) ? T(NaN) : d/v for (z,v,d) in zip(points,values.value,values.derivative)]
        output = (;output...,logderivative=ratio)
    end
    return output
end

# Coefficients of A*y'' + B*y' + C*y = 0, in the package's lambda convention.
function _wave_logderivative(value,derivative,x,m,spheroid,target,kind)
    if target === :angular && abs(x)==1 && (m>0 || kind==2)
        return oftype(value,copysign(Inf,kind==1 ? -x : x))
    elseif target === :radial && kind==1 &&
           ((spheroid === :prolate && x==1 && m>0) || (spheroid === :oblate && x==0 && iszero(value) && !iszero(derivative)))
        return oftype(value,Inf)
    end
    return iszero(value) ? oftype(value,NaN) : derivative/value
end

function _coordinate_derivatives(m,n,c,points,spheroid,precision,target,kind,normalize,
                                 scaled,logderivative,order;normalization=:standard)
    result = target === :angular ? smn(m,n,c,points;spheroid,precision,kind,normalize,scaled,
                                       logderivative,second_derivative=true) :
                                   rmn(m,n,c,points;spheroid,precision,kind,scaled,normalization,
                                       logderivative,second_derivative=true)
    R = precision === :quad ? BigFloat : Float64
    T = target === :angular && c isa Real ? R : Complex{R}
    third,fourth = _with_swprecision(max(320,_angular_precision(c))) do
        swfloat(z) = z isa Real ? _SWFloat(z) : Complex{_SWFloat}(z)
        wave_value(v,i) = scaled ? swfloat(v.mantissa[i])*_SWFloat(10)^v.exponent[i] : swfloat(v[i])
        parameter = swfloat(c)
        lambda = swfloat(eigenvalue(m,n,c;spheroid,precision))
        thirds,fourths = swfloat.(zeros(T,length(points))),swfloat.(zeros(T,length(points)))
        sigma = spheroid === :prolate ? 1 : -1
        for (i,point) in enumerate(points)
            x = _SWFloat(point)
            singular = target === :angular ? abs(x)==1 : spheroid === :prolate && x==1
            if singular
                if kind==1
                    option = target === :angular ? Int(normalize) : 1
                    amplitude,p = _regular_endpoint_factor(m,n,parameter,x,lambda,spheroid,precision,target,option;normalization)
                    thirds[i] = _endpoint_derivative(m,3,amplitude,p,x,target,c isa Real,_SWFloat)
                    fourths[i] = _endpoint_derivative(m,4,amplitude,p,x,target,c isa Real,_SWFloat)
                else
                    fourths[i] = wave_value(result.second_derivative,i)
                    thirds[i] = sign(x)*fourths[i]
                end
                continue
            end
            near = target === :angular ? 1-abs(x)<1//65536 :
                                         spheroid === :prolate && x-1<1//65536
            if kind==1 && near
                thirds[i],fourths[i] = _endpoint_derivatives(m,n,parameter,x,lambda,
                                                            spheroid,precision,target,normalize;normalization)
                continue
            end
            y,dy,ddy = wave_value(result.value,i),wave_value(result.derivative,i),wave_value(result.second_derivative,i)
            a,b,d = _wave_equation(m,parameter,x,lambda,spheroid,target)
            bprime = target === :angular ? -2 : 2
            q = (target === :angular ? -sigma : 1)*parameter^2
            mass = (target === :angular ? 1 : sigma)*m^2
            dprime = 2q*x+mass*b/a^2
            dsecond = 2q+mass*(bprime/a^2-2b^2/a^3)
            thirds[i] = -(2b*ddy+(bprime+d)*dy+dprime*y)/a
            fourths[i] = -(3b*thirds[i]+(3bprime+d)*ddy+2dprime*dy+dsecond*y)/a
        end
        thirds,fourths
    end
    convert_values(v) = scaled ? _decimal_scaled(v,T) : T.(v)
    output = (;result...,third_derivative=convert_values(third))
    return order==4 ? (;output...,fourth_derivative=convert_values(fourth)) : output
end

# Differentiate the regular endpoint factor without subtracting singular ODE terms.
function _endpoint_derivatives(m,n,c,x,lambda,spheroid,precision,target,normalize;normalization=:standard)
    option = target === :angular ? Int(normalize) : 1
    amplitude,coefficients = _regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option;
                                                      radius=abs(1-abs(x)),normalization)
    swfloat(z) = z isa Real ? _SWFloat(z) : Complex{_SWFloat}(z)
    polynomial = fill(zero(swfloat(amplitude)),5)
    t = 1-abs(x)
    for coefficient in reverse(coefficients)
        for j in 4:-1:1
            polynomial[j+1] = polynomial[j+1]*t+j*polynomial[j]
        end
        polynomial[1] = polynomial[1]*t+swfloat(coefficient)
    end
    sign = target === :angular ? -1 : 1
    u = sign*(abs(x)-1)*(abs(x)+1)
    du,ddu = 2sign*abs(x),2sign
    exponent = m//2
    factors = [u^exponent]
    for k in 0:3
        previous = k==0 ? zero(u) : factors[k]
        push!(factors,((exponent-k)*du*factors[k+1]+(exponent*k-k*(k-1)/2)*ddu*previous)/u)
    end
    direction = target === :angular && x<0 ? -1 : 1
    result = ntuple(2) do i
        k = i+2
        swfloat(amplitude)*direction^k*sum(binomial(k,j)*factors[j+1]*(-1)^(k-j)*polynomial[k-j+1] for j in 0:k)
    end
    return c isa Real ? real.(result) : result
end

function _wave_equation(m, c, x, lambda, spheroid, target)
    sigma = spheroid === :prolate ? 1 : -1
    if target === :angular
        a = (1-x)*(1+x)
        return a, -2x, lambda-sigma*c^2*x^2-m^2/a
    end
    a = spheroid === :prolate ? (x-1)*(x+1) : x^2+1
    return a, 2x, c^2*x^2-lambda-sigma*m^2/a
end

function _wave_second_derivative(result, m, n, c, points, spheroid, precision, target; option=0,normalization=:standard)
    lambda = eigenvalue(m,n,c;spheroid,precision)
    T = precision === :quad ? BigFloat : Float64
    second = similar(result.value)
    for (i, point) in enumerate(points)
        x = _input_bigfloat(point)
        singular = target === :angular ? abs(x) == 1 : spheroid === :prolate && x == 1
        if singular
            second[i] = _wave_endpoint_second(result.value[i],result.derivative[i],m,n,c,x,
                                               lambda,spheroid,precision,target,option;normalization)
        elseif (target === :angular && 1-abs(x) < 1//65536) ||
               (target === :radial && spheroid === :prolate && option == 1 && x-1 < 1//65536)
            second[i] = _regular_near_endpoint_second(m,n,c,x,lambda,spheroid,precision,target,option;normalization)
        else
            a,b,d = _wave_equation(m,c,x,lambda,spheroid,target)
            second[i] = -(b*result.derivative[i]+d*result.value[i])/a
        end
    end
    return (;result...,second_derivative=second)
end

function _wave_endpoint_second(value,derivative,m,n,c,x,lambda,spheroid,precision,target,option;normalization=:standard)
    T = precision === :quad ? BigFloat : Float64
    target === :radial && option != 1 && return T(NaN)
    m > 4 && return zero(value)
    amplitude,polynomial = _regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option;normalization)
    return _endpoint_derivative(m,2,amplitude,polynomial,x,target,c isa Real,T)
end

function _endpoint_derivative(m,order,amplitude,coefficients,x,target,real_result,T)
    m > 2order && return zero(amplitude)
    direction = target === :angular ? -sign(x) : 1
    if isodd(m)
        leading = amplitude*direction^order*prod(m-2j for j in 0:order-1)/2^order
        return _directed_infinity(leading,real_result,T)
    end
    power = m÷2
    sigma = target === :angular ? -1 : 1
    coefficient = sum(binomial(power,j)*2^(power-j)*sigma^j*(-sigma)^(order-power-j)*
        get(coefficients,order-power-j+1,zero(amplitude)) for j in 0:min(power,order-power))
    value = amplitude*direction^order*factorial(order)*coefficient
    return real_result ? real(value) : value
end

function _directed_infinity(direction,real_result,T)
    component(z) = iszero(z) ? zero(T) : copysign(T(Inf),z)
    return real_result ? component(real(direction)) : complex(component(real(direction)),component(imag(direction)))
end

function _regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option; radius=0,normalization=:standard)
    # Recover the regular factor at the endpoint from a nonsingular binary
    # anchor. The same series applies radially with t=1-x < 0.
    q = (spheroid === :prolate ? 1 : -1)*c^2
    L = lambda isa Real ? BigFloat(lambda) : Complex{BigFloat}(lambda)
    Q = q isa Real ? BigFloat(q) : Complex{BigFloat}(q)
    distance = BigFloat(2)^(-ceil(Int,log2(max(big"256",32*(1+abs(L)+abs(Q)+m*(m+1))))))
    point = target === :angular ? sign(x)*(1-distance) : 1+distance
    r = target === :radial && normalization === :static ? _static_radial_values(m,n,c,[point],spheroid,precision,1) :
        _scaled_native_values(m,n,c,[point],spheroid,precision,target,option)
    polynomial = _angular_endpoint_coefficients(m,L,Q,max(distance,radius))
    t = target === :angular ? distance : -distance
    shape = _angular_endpoint_polynomial(polynomial,t).value
    factor = distance*(target === :angular ? 2-distance : 2+distance)
    amplitude = only(r.value)/(factor^(BigFloat(m)/2)*shape)
    return amplitude,polynomial
end

function _regular_near_endpoint_second(m,n,c,x,lambda,spheroid,precision,target,option;normalization=:standard)
    amplitude,polynomial = _regular_endpoint_factor(m,n,c,x,lambda,spheroid,precision,target,option;radius=abs(1-abs(x)),normalization)
    coordinate = abs(x)
    t = 1-coordinate
    shape = _angular_endpoint_polynomial(polynomial,t)
    u = target === :angular ? (1-coordinate)*(1+coordinate) : (coordinate-1)*(coordinate+1)
    sigma = target === :angular ? 1 : -1
    result = amplitude*u^(BigFloat(m)/2)*(shape.second_derivative+2sigma*m*coordinate/u*shape.derivative+
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
    xs = _input_float.(T,points)
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
        step = h === nothing ? epsilon^(T(1)/10)*max(one(T),abs(x))/4 : _input_float(T,h)
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
