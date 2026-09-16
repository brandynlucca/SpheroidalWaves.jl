# F_c f(x) = integral(exp(-i*c*x*t)*f(t), t=-1..1).
# With orthonormal Legendre coefficients v, evaluate this identity at x=0
# for even modes, or differentiate it at x=0 for odd modes. Orthogonality
# leaves just v[1] in the integral. No oscillatory quadrature is needed here.
function _integral_eigenvalue_data(n,c,precision; sensitivity=false)
    bits = _angular_working_bits(c)
    previous = nothing
    tolerance = precision===:quad ? big"1e-34" : big"1e-17"
    for attempt in 1:4
        result = setprecision(BigFloat,bits) do
            # Avoid an unresolved native eigenvalue at extremely small c.
            # The perturbation has norm <= c^2, so the spherical eigenvalue
            # identifies the same isolated mode well inside the seed tolerance.
            seed = c<big"1e-20" ? BigFloat(n)*(n+1) : nothing
            plan = _guarded_angular_plan(0,n,c,:prolate,precision;eigenvalue_seed=seed)
            origin = _evaluate_coefficient_vector(plan,plan.v,[0])
            ratio = if iseven(n)
                sqrt(big"2.0")*first(plan.v)/only(origin.value)
            else
                BigFloat(c)*sqrt(big"2.0"/3)*first(plan.v)/only(origin.derivative)
            end
            amplitude = abs(ratio)
            concentration = BigFloat(c)*amplitude^2/(2BigFloat(pi))
            complement = 1-concentration
            if sensitivity
                endpoint = only(_evaluate_coefficient_vector(plan,plan.v,[1]).value)
                # Differentiate the sinc kernel: its quadratic form is
                # |integral(exp(i*c*t)*psi(t),t=-1..1)|^2/pi.
                # The Fourier identity at x=1 then gives Lambda'/Lambda.
                dlog = 2endpoint^2/BigFloat(c)
                dconcentration = amplitude^2*endpoint^2/BigFloat(pi)
                damplitude = if n==0 && c<1//2
                    # Avoid subtracting endpoint^2 - 1/2 = O(c^2).
                    dorigin = _evaluate_coefficient_vector(plan,plan.dv,[0])
                    amplitude*(first(plan.dv)/first(plan.v)-only(dorigin.value)/only(origin.value))
                else
                    amplitude*(endpoint^2-1//2)/BigFloat(c)
                end
                (;amplitude,concentration,complement,dlog,dconcentration,damplitude)
            else
                (;amplitude,concentration,complement)
            end
        end
        valid = all(isfinite,values(result)) && result.amplitude>0 &&
                0<result.concentration<1 && result.complement>0
        if valid && previous !== nothing &&
                all(abs(a-b)<=tolerance*abs(a) for (a,b) in zip(values(result),values(previous)))
            return result
        end
        previous = result
        bits += 128
    end
    error("Integral-operator eigenvalue did not converge, including its distance from zero and one")
end

function _integral_parameter(m,n,c,spheroid,precision,operator,form)
    m==0 && n>=0 || throw(ArgumentError("integral operators require m=0 and n>=0"))
    spheroid===:prolate || throw(ArgumentError("integral operators require spheroid=:prolate"))
    c isa Real && isfinite(c) && c>=0 ||
        throw(DomainError(c,"integral operators require finite real c >= 0"))
    form in (:value,:log,:complement) || throw(ArgumentError("form must be :value, :log or :complement"))
    operator===:fourier && form!==:value &&
        throw(ArgumentError("the Fourier operator supports only form=:value"))
    T = precision===:quad ? BigFloat : Float64
    parameter = T(c)
    isfinite(parameter) || throw(DomainError(c,"c is not finite at the requested precision"))
    !iszero(c) && iszero(parameter) &&
        throw(DomainError(c,"c rounds to zero at the requested precision; use precision=:quad"))
    return parameter
end

function _integral_eigenvalue(m,n,c,spheroid,precision,operator,form)
    parameter = _integral_parameter(m,n,c,spheroid,precision,operator,form)
    T = typeof(parameter)
    if iszero(parameter)
        operator===:fourier && return complex(n==0 ? T(2) : zero(T))
        return form===:complement ? one(T) : form===:log ? T(-Inf) : zero(T)
    end
    result = _integral_eigenvalue_data(n,parameter,precision)
    if operator===:fourier
        phase = (1,-im,-1,im)[mod(n,4)+1]
        return Complex{T}(phase)*T(result.amplitude)
    end
    form===:value && return T(result.concentration)
    form===:complement && return T(result.complement)
    # Retain a tiny negative logarithm even when concentration rounds to one.
    logarithm = setprecision(BigFloat,Base.precision(result.concentration)) do
        result.concentration>1//2 ? log1p(-result.complement) : log(result.concentration)
    end
    return T(logarithm)
end

function _integral_eigen_jacobian(m,n,c,spheroid,precision,operator,form,
                                  h,with_metadata,adaptive,rtol,atol)
    parameter = _integral_parameter(m,n,c,spheroid,precision,operator,form)
    T = typeof(parameter)
    if h !== nothing
        step = _resolve_jacobian_step(parameter,h,precision)
        iszero(parameter) && form===:log &&
            throw(DomainError(parameter,"log concentration is singular at c=0; omit h for its derivative limit"))
        !iszero(parameter) && step>=parameter &&
            throw(ArgumentError("h must be smaller than c for a centered integral-operator difference"))
        f(x) = _integral_eigenvalue(m,n,x,spheroid,precision,operator,form)
        calc(s) = iszero(parameter) ? (-3f(parameter)+4f(parameter+s)-f(parameter+2s))/(2s) :
                                     (f(parameter+s)-f(parameter-s))/(2s)
        derivative,metadata = _finite_difference_with_metadata(calc,step;precision,adaptive,rtol,atol)
        return with_metadata ? (;derivative,metadata) : derivative
    end
    if iszero(parameter)
        derivative = if operator===:fourier
            n==1 ? complex(zero(T),-T(2)/3) : complex(zero(T))
        elseif form===:log
            T(Inf)
        else
            (form===:complement ? -1 : 1)*(n==0 ? T(2)/T(pi) : zero(T))
        end
    else
        result = _integral_eigenvalue_data(n,parameter,precision;sensitivity=true)
        derivative = operator===:fourier ? Complex{T}((1,-im,-1,im)[mod(n,4)+1])*T(result.damplitude) :
                     form===:log ? T(result.dlog) :
                     T(form===:complement ? -result.dconcentration : result.dconcentration)
    end
    metadata = (method=iszero(parameter) ? :right_limit : :integral_identity,
                step_used=nothing,relative_change_when_halving_step=nothing,
                finite_flag=isfinite(derivative),
                conditioning_flag=isfinite(derivative) ? :good : :singular,
                suggested_action=isfinite(derivative) ? :accept : :singular_limit)
    return with_metadata ? (;derivative,metadata) : derivative
end

# Invert log(Lambda/(1-Lambda)), which remains sensitive in either tail.
# The user's target and returned residual retain the requested form.
function _find_integral_bandwidth(m,n,target,bracket,spheroid,precision,form,
                                   atol,rtol,maxiter,use_jacobian)
    T = precision===:quad ? BigFloat : Float64
    a = _integral_parameter(m,n,bracket[1],spheroid,precision,:concentration,form)
    b = _integral_parameter(m,n,bracket[2],spheroid,precision,:concentration,form)
    a<b || throw(ArgumentError("bracket must satisfy 0 <= c_lo < c_hi"))
    y = T(target)
    !iszero(target) && iszero(y) && throw(DomainError(target,"target rounds to zero; use precision=:quad or a logarithmic target"))
    isfinite(target) && !isfinite(y) && throw(DomainError(target,"target is not finite at the requested precision; use precision=:quad"))
    form===:complement && isone(y) && !isone(target) &&
        throw(DomainError(target,"complement rounds to one; specify concentration directly or use precision=:quad"))
    valid = form===:value ? 0<=y<1 : form===:complement ? 0<y<=1 : y<0
    valid && !isnan(y) || throw(DomainError(target,
        "concentration targets require 0 <= value < 1, 0 < complement <= 1, or log < 0; use complement for targets near one"))
    zero_target = form===:value ? iszero(y) : form===:complement ? isone(y) : y==-Inf
    if zero_target
        iszero(a) || throw(ArgumentError("zero concentration requires a bracket containing c=0"))
        return (converged=true,c=a,residual=zero(T),iterations=0,bracket=(a,b),method=:endpoint)
    end
    isfinite(y) || throw(DomainError(target,"target must be finite except log=-Inf"))
    t = BigFloat(y)
    odds = form===:value ? log(t)-log1p(-t) :
           form===:complement ? log1p(-t)-log(t) : t-log(-expm1(t))
    function evaluate(c)
        if iszero(c)
            value = form===:value ? zero(T) : form===:complement ? one(T) : T(-Inf)
            return (f=BigFloat(-Inf),d=BigFloat(Inf),residual=value-y)
        end
        result = _integral_eigenvalue_data(n,c,precision;sensitivity=use_jacobian)
        setprecision(BigFloat,Base.precision(result.concentration)) do
            logvalue = result.concentration>1//2 ? log1p(-result.complement) : log(result.concentration)
            value = form===:value ? result.concentration : form===:complement ? result.complement : logvalue
            (f=logvalue-log(result.complement)-odds,
             d=use_jacobian ? result.dlog/result.complement : BigFloat(NaN),
             residual=T(value-BigFloat(y)))
        end
    end
    left,right = evaluate(a),evaluate(b)
    for (c,r) in ((a,left),(b,right))
        (iszero(r.f) || (abs(r.f)<=rtol && (!use_jacobian || abs(r.f/r.d)<=atol+rtol*abs(c)))) &&
            return (converged=true,c,residual=r.residual,iterations=0,bracket=(a,b),method=:endpoint)
    end
    left.f<=0<=right.f || throw(ArgumentError("bracket endpoints must straddle the requested concentration"))
    c = (a+b)/2
    method = :bisection
    for iteration in 1:maxiter
        result = evaluate(c)
        xtol = T(atol)+T(rtol)*abs(c)
        # A small absolute eigenvalue residual alone is insufficient in a tail.
        if iszero(result.f) || (abs(result.f)<=rtol &&
           (b-a<=xtol || (use_jacobian && abs(result.f/result.d)<=xtol)))
            return (converged=true,c,residual=result.residual,iterations=iteration,bracket=(a,b),method)
        end
        if result.f<0
            a,left = c,result
        else
            b,right = c,result
        end
        iteration==maxiter && return (converged=false,c,residual=result.residual,
                                      iterations=iteration,bracket=(a,b),method=:maxiter)
        midpoint = (a+b)/2
        candidate = use_jacobian ? T(c-result.f/result.d) : midpoint
        # Periodic bisection guarantees progress even for poor Newton steps.
        if iteration%3!=0 && isfinite(candidate) && a<candidate<b && candidate!=c
            c,method = candidate,use_jacobian ? :newton : :bisection
        else
            c,method = midpoint,:bisection
        end
        (c==a || c==b) && return (converged=b-a<=xtol && abs((c==a ? left : right).f)<=rtol,c,
            residual=(c==a ? left : right).residual,iterations=iteration,bracket=(a,b),method=:bisection)
    end
end
