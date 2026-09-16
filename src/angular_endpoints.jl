# The regular angular solution is S(x) = (1-x^2)^(m/2) F(1-x).
# Evaluate its analytic factor F in BigFloat near the endpoint, keeping the
# requested distance from one instead of rounding x to native binary128 first.
_near_angular_endpoint(x) = 0 < 1-abs(x) < 1//65536

function _angular_endpoint_plan(m, n, c, points, spheroid; precision=:quad)
    indices = findall(_near_angular_endpoint, points)
    isempty(indices) && return nothing
    parameter = c isa Real ? BigFloat(c) : Complex{BigFloat}(c)
    # This plan is attached to a raw native evaluation. Its degree has already
    # been selected by the caller; do not continue that label a second time.
    native_lambda = parameter isa Complex ?
        _call_complex_eigenvalue(spheroid === :prolate ? :cprolate : :coblate,m,n,parameter;precision) :
        eigenvalue(m,n,parameter;spheroid,precision)
    lambda = native_lambda isa Real ? BigFloat(native_lambda) : Complex{BigFloat}(native_lambda)
    q = (spheroid === :prolate ? 1 : -1)*parameter^2
    # An exactly representable anchor sufficiently close to one keeps F away
    # from its zeros and avoids cancellation when setting the normalization.
    scale = max(big"256",32*(1+abs(lambda)+abs(q)+BigFloat(m)*(m+1)))
    anchor_distance = BigFloat(2)^(-ceil(Int,log2(scale)))
    anchor = 1-anchor_distance
    anchor == 1 && error("Angular endpoint anchor exceeds the native coordinate resolution")
    native_points = BigFloat.(points)
    native_points[indices] .= anchor
    radius = max(anchor_distance,maximum(i -> 1-abs(BigFloat(points[i])),indices))
    coefficients = _angular_endpoint_coefficients(m,lambda,q,radius)
    (;indices,points=BigFloat.(points),native_points,coefficients,anchor_distance)
end

function _angular_endpoint_coefficients(m,lambda,q,radius)
    coefficients = [one(lambda+q)]
    term_sum, derivative_sum = one(lambda+q),zero(lambda+q)
    power, small = one(radius),0
    for k in 0:1023
        current = coefficients[k+1]
        previous = k >= 1 ? coefficients[k] : zero(current)
        previous2 = k >= 2 ? coefficients[k-1] : zero(current)
        next = ((BigFloat(k)*(k+2m+1)-lambda+BigFloat(m)*(m+1)+q)*current -
                2q*previous+q*previous2)/(BigFloat(2)*(k+1)*(k+m+1))
        push!(coefficients,next)
        derivative_term = (k+1)*next*power
        power *= radius
        term = next*power
        term_sum += term
        derivative_sum += derivative_term
        tiny = abs(term) <= eps(BigFloat)*abs(term_sum) &&
               abs(derivative_term) <= eps(BigFloat)*abs(derivative_sum)
        small = tiny ? small+1 : 0
        small >= 3 && break
        k == 1023 && error("Angular endpoint series did not converge")
    end
    return coefficients
end

function _angular_endpoint_polynomial(coefficients, distance)
    value, derivative, second = last(coefficients),zero(last(coefficients)),zero(last(coefficients))
    for k in length(coefficients)-1:-1:1
        second = second*distance+2derivative
        derivative = derivative*distance+value
        value = value*distance+coefficients[k]
    end
    (;value,derivative,second_derivative=second)
end

function _angular_endpoint_reconstruct!(value, derivative, plan, m, n)
    plan === nothing && return nothing
    d0 = plan.anchor_distance
    anchor_shape = _angular_endpoint_polynomial(plan.coefficients,d0).value
    anchor_factor = (d0*(2-d0))^(m//2)
    for i in plan.indices
        x = abs(plan.points[i])
        distance = 1-x
        u = distance*(1+x)
        shape = _angular_endpoint_polynomial(plan.coefficients,distance)
        amplitude = value[i]/(anchor_factor*anchor_shape)
        factor = u^(m//2)
        value[i] = amplitude*factor*shape.value
        derivative[i] = -amplitude*factor*(shape.derivative+m*x/u*shape.value)
        if plan.points[i] < 0
            value[i] *= isodd(n-m) ? -1 : 1
            derivative[i] *= isodd(n-m) ? 1 : -1
        end
    end
    return nothing
end
