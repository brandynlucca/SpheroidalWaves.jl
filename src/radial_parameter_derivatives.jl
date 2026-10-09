# Differentiate DLMF 30.11.3, including the joining sum in its denominator.
# No parameter differences or native radial values enter this calculation.
function _static_precision(m, n, parameter)
    # Inverse iteration resolves small coefficients at roughly one third of
    # the arithmetic precision. Static sensitivities cancel another c^2.
    guards = iszero(parameter) ? 0 : 6max(0, -exponent(_input_bigfloat(abs(parameter))))
    max(320, Base.precision(_SWFloat), 320+2(n+m)+guards)
end

function _static_radial_values(
        m, n, c, points, spheroid, precision, kind; eigenvalue_seed = nothing)
    T = precision === :quad ? BigFloat : Float64
    parameter = c isa Real ? _input_float(T, c) : _input_float(Complex{T}, c)
    return _with_swprecision(_static_precision(m, n, parameter)) do
        if !iszero(parameter)
            data, plan = _radial_analytic_data(
                m, n, parameter, points, spheroid, precision, kind; eigenvalue_seed)
            z = parameter isa Real ? _SWFloat(parameter) : Complex{_SWFloat}(parameter)
            factor = z^(kind==1 ? -n : n+1)
            value = [factor*s[1] for s in data.states]
            derivative = [factor*s[2] for s in data.states]
            if spheroid === :prolate && kind==1 && m==1 && any(isone, points)
                regular = _radial_expansion_data(plan, one(_SWFloat), 1; regular_factor = true).state[1]
                derivative[findall(isone, points)] .= _directed_infinity(
                    factor*regular, c isa
                                    Real, _SWFloat)
            end
            return (; value, derivative)
        end
        sigma = spheroid === :prolate ? 1 : -1
        value, derivative = (zeros(_SWFloat, length(points)) for _ in 1:2)
        if kind==1
            for (i, x) in enumerate(points)
                value[i], derivative[i] = _static_regular_radial(m, n, _SWFloat(x), sigma)
            end
        else
            anchor = _SWFloat(2)
            y, dy = _static_irregular_radial(m, n, anchor, sigma)
            state = (y, dy, zero(y), zero(y))
            plan = (;
                m, c = zero(y), lambda = _SWFloat(n)*(n+1), dlambda = zero(y), spheroid)
            for i in sortperm(points; rev = true)
                x = _SWFloat(points[i])
                if spheroid === :prolate && x==1
                    value[i] = derivative[i] = _SWFloat(NaN)
                elseif x>=2
                    value[i], derivative[i] = _static_irregular_radial(m, n, x, sigma)
                else
                    while anchor>x
                        radius = sigma==1 ? anchor-1 : sqrt(anchor^2+1)
                        step = -min(anchor-x, radius/4, inv(_SWFloat(n+1)))
                        state = _radial_sensitivity_step(plan, anchor, state, step)
                        anchor += step
                    end
                    value[i], derivative[i] = state[1], state[2]
                end
            end
        end
        (; value, derivative)
    end
end

# The growing Laplace solution has leading term x^n/(2n+1)!!.
function _static_regular_radial(m, n, x, sigma)
    degree = n-m
    coefficients = zeros(typeof(x), degree+1)
    term = inv(prod(typeof(x)(k) for k in 1:2:(2n + 1); init = one(x)))
    for k in 0:(degree ÷ 2)
        coefficients[degree - 2k + 1] = term
        k<degree÷2 && (term *= -sigma*typeof(x)(n-k)*(degree-2k)*(degree-2k-1) /
                 (typeof(x)(k+1)*(2n-2k)*(2n-2k-1)))
    end
    polynomial, derivative = last(coefficients), zero(x)
    for coefficient in reverse(coefficients[1:(end - 1)])
        derivative = derivative*x+polynomial
        polynomial = polynomial*x+coefficient
    end
    if sigma==1 && x==1
        slope = m==0 ? derivative :
                m==1 ? _directed_infinity(polynomial, true, typeof(x)) :
                m==2 ? 2polynomial : zero(x)
        return m==0 ? polynomial : zero(x), slope
    end
    u = x^2-sigma
    factor = u^(m//2)
    slope = m==0 ? zero(x) : m*x/u
    return factor*polynomial, factor*(derivative+slope*polynomial)
end

# DLMF 14.3.7, normalized to -(2n-1)!!/x^(n+1) at infinity.
function _static_irregular_radial(m, n, x, sigma; max_terms = 4096)
    a, b, d = (typeof(x)(n+m+1)/2, typeof(x)(n+m+2)/2, typeof(x)(n)+3//2)
    term, total, derivative = one(x), one(x), zero(x)
    for k in 1:max_terms
        term *= (a+k-1)*(b+k-1)*sigma/((d+k-1)*k*x^2)
        total += term
        derivative -= 2k*term/x
        if abs(term)<=eps(typeof(x))*abs(total) &&
           abs(2k*term/x)<=eps(typeof(x))*abs(derivative)
            u = x^2-sigma
            factor = -prod(typeof(x)(j) for j in 1:2:(2n - 1); init = one(x))*u^(m//2)*x^(-n-m-1)
            return factor*total, factor*(derivative+(m*x/u-(n+m+1)/x)*total)
        end
    end
    error("Static radial series did not converge")
end

function _radial_spherical_bessel(top, z, kind)
    b = zeros(kind>=3 ? typeof(complex(z)) : typeof(z), top+2)
    if kind >= 3
        # DLMF 10.49.6-7 at degrees zero and one. Preserve the outgoing or
        # incoming exponential directly, including when it is very small.
        phase = kind==3 ? im : -im
        wave = exp(phase*z)/z
        b[1], b[2] = -phase*wave, -wave*(1+phase/z)
        for l in 1:top
            b[l + 2] = (2l+1)/z*b[l + 1]-b[l]
        end
    elseif kind == 2
        b[1], b[2] = -cos(z)/z, -cos(z)/z^2-sin(z)/z
        for l in 1:top
            b[l + 2] = (2l+1)/z*b[l + 1]-b[l]
        end
    elseif abs(z)>top+32
        b[1], b[2] = sin(z)/z, (sin(z)/z-cos(z))/z
        for l in 1:top
            b[l + 2] = (2l+1)/z*b[l + 1]-b[l]
        end
    else
        # Miller recurrence selects the recessive solution at large degree.
        # Normalize at whichever of j0,j1 has the larger magnitude.
        start = top+32+ceil(Int, abs(z))+cld(Base.precision(_SWFloat), 4)
        next, current = zero(z), one(z)
        for l in start:-1:1
            previous = (2l+1)/z*current-next
            l-1<=top+1 && (b[l]=previous)
            next, current = current, previous
        end
        j0, j1 = sin(z)/z, (sin(z)/z-cos(z))/z
        scale = abs(j0)>=abs(j1) ? j0/b[1] : j1/b[2]
        b .*= scale
    end
    return b
end

function _radial_expansion_data(plan, x, kind; regular_factor = false)
    m, n, c = plan.m, plan.n, plan.c
    sigma = plan.spheroid === :prolate ? 1 : -1
    z = c*x
    b = _radial_spherical_bessel(last(plan.degrees), z, kind)
    factors = [sqrt(prod(_SWFloat(k) for k in (l - m + 1):(l + m); init = _SWFloat("1"))*(2l+1)/2)
               for l in plan.degrees]
    weights, dweights = factors .* plan.v, factors .* plan.dv
    denominator, ddenominator = sum(weights), sum(dweights)
    abs(denominator)>eps(_SWFloat)^(1//2)*sum(abs, weights) ||
        error("Radial joining sum is singular or numerically unresolved")
    sums, tail = zeros(eltype(b), 4), zeros(_SWFloat, 4)
    for (i, l) in enumerate(plan.degrees)
        w, dw = (-1)^((l-n)÷2)*weights[i], (-1)^((l-n)÷2)*dweights[i]
        v = b[l + 1]
        d = l/z*v-b[l + 2]
        dd = (l*(l+1)/z^2-1)*v-2d/z
        terms = (w*v, w*c*d, dw*v+w*x*d, dw*c*d+w*(d+c*x*dd))
        for k in 1:4
            sums[k] += terms[k]
            i>length(weights)-8 && (tail[k]+=abs(terms[k]))
        end
    end
    # A fixed absolute floor would hide truncation error in decaying waves.
    relative_tail = maximum(iszero(s) ? (iszero(t) ? zero(t) : oftype(t, Inf)) :
                            t/abs(s) for (t, s) in zip(tail, sums))
    f, df = sums[1]/denominator, sums[2]/denominator
    tangent = (sums[3]-ddenominator*f)/denominator
    dtangent = (sums[4]-ddenominator*df)/denominator
    regular_factor && return (; state = (f, df, tangent, dtangent), relative_tail)
    if sigma==1 && x==1 && m>0
        slope(v) = m==1 ? _directed_infinity(v, c isa Real, _SWFloat) : m==2 ? 2v : zero(v)
        return (; state = (zero(f), slope(f), zero(tangent), slope(tangent)), relative_tail)
    end
    factor = (1-sigma/x^2)^(m//2)
    slope = sigma*m/(x*(x^2-sigma))
    # For m=0 at the regular prolate boundary this is exactly zero, not 0/0.
    m==0 && (slope=zero(x))
    return (;
        state = (
            factor*f, factor*(df+slope*f), factor*tangent, factor*(dtangent+slope*tangent)),
        relative_tail)
end

# Polynomial form of the radial equation and its parameter derivative.
# A two-entry state propagates the value and coordinate slope. A four-entry
# state also propagates their parameter sensitivities.
function _radial_sensitivity_step(plan, x, state, h; rtol = _SWFloat("1e-42"),
        max_terms = max(256, cld(Base.precision(_SWFloat), 2)+64))
    y, dy = state[1:2]
    sensitivity = length(state) == 4
    sigma = plan.spheroid === :prolate ? 1 : -1
    u, q, dq = x^2-sigma, get(plan, :q, plan.c^2), get(plan, :dq, 2plan.c)
    lambda, dlambda = plan.lambda, plan.dlambda
    a = (u^2, 4x*u, 6x^2-2sigma, 4x, one(x))
    b = (2x*u, 6x^2-2sigma, 6x, _SWFloat("2"))
    d = ((q*x^2-lambda)*u-sigma*plan.m^2, 4q*x^3-2*(sigma*q+lambda)*x,
        6q*x^2-sigma*q-lambda, 4q*x, q)
    dc = sensitivity ?
         ((dq*x^2-dlambda)*u, 4dq*x^3-2*(sigma*dq+dlambda)*x,
        6dq*x^2-sigma*dq-dlambda, 4dq*x, dq) : nothing
    values = [y, dy]
    tangents = sensitivity ? [state[3], state[4]] : nothing
    for k in 0:(max_terms - 2)
        rhs = zero(y)
        drhs = sensitivity ? zero(state[3]) : nothing
        for j in 1:min(4, k)
            factor = a[j + 1]*(k-j+2)*(k-j+1)
            rhs += factor*values[k - j + 3]
            sensitivity && (drhs += factor*tangents[k - j + 3])
        end
        for j in 0:min(3, k)
            factor = b[j + 1]*(k-j+1)
            rhs += factor*values[k - j + 2]
            sensitivity && (drhs += factor*tangents[k - j + 2])
        end
        for j in 0:min(4, k)
            rhs += d[j + 1]*values[k - j + 1]
            sensitivity &&
                (drhs += d[j + 1]*tangents[k - j + 1]+dc[j + 1]*values[k - j + 1])
        end
        push!(values, -rhs/(a[1]*(k+1)*(k+2)))
        sensitivity && push!(tangents, -drhs/(a[1]*(k+1)*(k+2)))
        if k>=30 && k%8==6
            order = length(values)-1
            radius = abs(h)
            first_power = radius^(order-8)
            output = map(sensitivity ? (values, tangents) : (values,)) do coefficients
                v, dv = last(coefficients), zero(y)
                for j in (order - 1):-1:0
                    dv = dv*h+v
                    v = v*h+coefficients[j + 1]
                end
                # Consecutive real powers avoid repeating complex exponentiation
                # for every term in the value and derivative tail bounds.
                tail, dtail, power = zero(radius), zero(radius), first_power
                for j in (order - 7):order
                    magnitude = abs(coefficients[j + 1])*power
                    dtail += j*magnitude
                    tail += radius*magnitude
                    power *= radius
                end
                good = tail<=rtol*max(abs(v), abs(coefficients[1]), abs(h*coefficients[2])) &&
                       dtail<=rtol*max(abs(dv), abs(coefficients[2]), abs(coefficients[1]/h))
                (; v, dv, good)
            end
            if all(r->r.good, output)
                return sensitivity ?
                       (output[1].v, output[1].dv, output[2].v, output[2].dv) :
                       (output[1].v, output[1].dv)
            end
        end
    end
    error("Differentiated radial equation did not converge")
end

function _radial_expansion_batch(plan, points, kind)
    V = kind>=3 ? typeof(complex(plan.c)) : typeof(plan.c)
    states = Vector{NTuple{4, V}}(undef, length(points))
    tail = _SWFloat("0")
    propagate = kind>=2 || plan.spheroid===:oblate
    interior = findall(x->x<2 && !(plan.spheroid===:prolate && x==1 && kind>=2), points)
    if propagate && !isempty(interior)
        initial = _radial_expansion_data(plan, _SWFloat("2"), kind)
        tail = max(tail, initial.relative_tail)
        x, state = _SWFloat("2"), initial.state
        for i in sort(interior; by = i->points[i], rev = true)
            target = _SWFloat(points[i])
            steps = 0
            while x>target
                radius = plan.spheroid===:prolate ? x-1 : sqrt(x^2+1)
                h = -min(x-target, radius/4, inv(1+abs(plan.c)+sqrt(abs(plan.lambda))))
                x+h<x ||
                    error("Radial sensitivity propagation exhausted coordinate precision")
                state = _radial_sensitivity_step(plan, x, state, h)
                x += h
                steps += 1
                steps<=100000 ||
                    error("Radial sensitivity propagation exceeded its step limit")
            end
            if plan.spheroid===:oblate && target==0 && kind==1
                # The regular solution has exact parity at the oblate origin.
                states[i] = isodd(plan.n-plan.m) ?
                            (zero(state[1]), state[2], zero(state[3]), state[4]) :
                            (state[1], zero(state[2]), state[3], zero(state[4]))
            else
                states[i] = state
            end
        end
    end
    for (i, x) in enumerate(points)
        if plan.spheroid===:prolate && x==1 && kind>=2
            states[i] = ntuple(_->V(NaN), 4)
        elseif !(propagate && i in interior)
            result = _radial_expansion_data(plan, _SWFloat(x), kind)
            states[i] = result.state
            tail = max(tail, result.relative_tail)
        end
    end
    return (; states, tail, propagated = propagate && !isempty(interior))
end

_radial_needs_analytic(c, points, kind) = c isa Complex &&
                                          (iszero(real(c)) || kind>=3)

function _radial_analytic_data(
        m, n, c, points, spheroid, precision, kind; eigenvalue_seed = nothing,
        coefficient_plan = _coefficient_plan, expansion_batch = _radial_expansion_batch, rtol = nothing)
    T = precision===:quad ? BigFloat : Float64
    parameter = c isa Real ? _input_float(T, c) : _input_float(Complex{T}, c)
    coordinates = _SWFloat.(points)
    # Direct Hankel waves do not lose exp(2*abs(imag(c*x))) through an R1 +/- iR2
    # subtraction. Keep guards for the coefficient solve and degree recurrence.
    cancellation = kind>=3 ? zero(real(parameter)) :
                   4abs(imag(parameter))*max(2, maximum(coordinates))
    bits = max(Base.precision(_SWFloat)+64, 256+2*(n+m)+ceil(Int, 8abs(parameter)+cancellation))
    result, plan = _with_swprecision(bits) do
        # Neumann terms converge geometrically at x=2, even when the angular
        # coefficients themselves have already become very small.
        minimum = max(kind==1 ? 32 : 96, (n-m)÷2+16+ceil(Int, abs(parameter)))
        for attempt in 1:4
            max_terms = max(512, 2minimum)
            tolerance = min(
                _SWFloat("1e-42"), exp(-2_SWFloat(abs(parameter)))*_SWFloat("1e-42"),
                rtol===nothing ? one(_SWFloat) : rtol)
            plan = coefficient_plan(m, n, parameter; spheroid, precision, rtol = tolerance,
                min_terms = minimum, max_terms, eigenvalue_seed)
            plan.converged || error("Radial sensitivity coefficients did not converge")
            result = expansion_batch(plan, coordinates, kind)
            result.tail<=(rtol===nothing ? _SWFloat("1e-34") : rtol) && return result, plan
            minimum *= 2
        end
        error("Differentiated radial expansion did not converge")
    end
    return result, plan
end

function _radial_analytic_values(m, n, c, points, spheroid, precision, kind;
        second_derivative = false, eigenvalue_seed = nothing)
    result, plan = _radial_analytic_data(
        m, n, c, points, spheroid, precision, kind; eigenvalue_seed)
    values = (;
        value = [s[1] for s in result.states], derivative = [s[2] for s in result.states])
    second_derivative || return values
    seconds = _with_swprecision(Base.precision(real(plan.c))) do
        [begin
             if spheroid === :prolate && x==1
                 _wave_endpoint_second(s[1], s[2], m, n, c, x, plan.lambda,
                     spheroid, precision, :radial, kind)
             elseif spheroid === :prolate && kind==1 && x-1<1//65536
                 _regular_near_endpoint_second(
                     m, n, c, _SWFloat(x), plan.lambda, spheroid, precision, :radial, kind)
             else
                 a, b, d=_wave_equation(
                     m, plan.c, _SWFloat(x), plan.lambda, spheroid, :radial)
                 -(b*s[2]+d*s[1])/a
             end
         end
         for (x, s) in zip(points, result.states)]
    end
    return (; values..., second_derivative = seconds)
end

function _radial_analytic_wave(m, n, c, points, spheroid, precision, kind, scaled,
        logderivative, second_derivative; normalization = :standard)
    result = if normalization === :static
        values = _static_radial_values(m, n, c, points, spheroid, precision, kind)
        second_derivative ?
        _wave_second_derivative(values, m, n, c, points, spheroid, precision,
            :radial; option = kind, normalization) : values
    else
        _radial_analytic_values(
            m, n, c, points, spheroid, precision, kind; second_derivative)
    end
    R=precision===:quad ? BigFloat : Float64
    output=map(
        v->scaled ? _decimal_scaled(v, Complex{R}) :
           complex.(R.(real.(v)), R.(imag.(v))), result)
    if logderivative
        ratios=Complex{R}[_wave_logderivative(v, d, x, m, spheroid, :radial, kind)
                          for (v, d, x) in zip(result.value, result.derivative, points)]
        output=(; output..., logderivative = ratios)
    end
    return output
end

function _radial_analytic_jacobian(m, n, c, points, spheroid, precision, kind,
        diagnostics; normalization = :standard)
    result, plan = if normalization === :static
        _with_swprecision(_static_precision(m, n, c)) do
            if iszero(c)
                plan = _sensitivity_plan(m, n, c, spheroid, precision)
                states = [ntuple(
                              _ -> spheroid===:prolate && x==1 && kind==2 ? _SWFloat(NaN) :
                                   zero(_SWFloat),
                              4) for x in points]
                return (; states, tail = zero(_SWFloat), propagated = false), plan
            end
            data, plan = _radial_analytic_data(m, n, c, points, spheroid, precision, kind)
            power = kind==1 ? -n : n+1
            states = map(zip(points, data.states)) do (x, state)
                f, df, fc, dfc = state
                factor = plan.c^power
                v = power==0 ? fc : factor*(fc+power*f/plan.c)
                d = power==0 ? dfc : factor*(dfc+power*df/plan.c)
                if spheroid === :prolate && x==1 && m==1 && kind==1
                    regular = _radial_expansion_data(plan, one(_SWFloat), 1; regular_factor = true).state
                    direction = factor*(regular[3]+power*regular[1]/plan.c)
                    d = _directed_infinity(direction, c isa Real, _SWFloat)
                end
                (factor*f, factor*df, v, d)
            end
            (; data..., states), plan
        end
    else
        _radial_analytic_data(m, n, c, points, spheroid, precision, kind)
    end
    T=precision===:quad ? BigFloat : Float64
    # Match rmn's complex output type even for real first/second-kind results.
    rounded(z) = complex(T(real(z)), T(imag(z)))
    value = [rounded(s[3]) for s in result.states]
    derivative = [rounded(s[4]) for s in result.states]
    metadata(v) = (; _coefficient_derivative_metadata(plan, v)...,
        method = result.propagated ? :differentiated_equation : :differentiated_expansion,
        radial_series_tail = _unwrap_swfloat(result.tail))
    mv, md = metadata(value), metadata(derivative)
    if c isa Real
        output = (; dvalue_dc = value, dderivative_dc = derivative)
        return diagnostics ?
               (; output..., metadata_value = mv, metadata_derivative = md) : output
    end
    output = (;
        dvalue_dcreal = value, dvalue_dcimag = complex.(-imag.(value), real.(value)),
        dderivative_dcreal = derivative, dderivative_dcimag = complex.(-imag.(derivative), real.(derivative)))
    return diagnostics ?
           (; output..., metadata_value_dcreal = mv, metadata_value_dcimag = mv,
        metadata_derivative_dcreal = md, metadata_derivative_dcimag = md) : output
end
