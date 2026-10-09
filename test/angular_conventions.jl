using SpheroidalWaves, Test

@testset "Complex angular coordinates" begin
    # Independent Ferrers polynomial recurrence fixes the square-root branch.
    function ferrers(m, degrees, z)
        previous = zero(z)
        current = (-1)^m*prod(oftype(z, k) for k in 1:2:(2m - 1); init = one(z))*(sqrt(1-z)*sqrt(1+z))^m
        values = typeof(z)[]
        for l in m:last(degrees)
            l in degrees && push!(values, current)
            previous, current = current, ((2l+1)*z*current-(l+m)*previous)/(l-m+1)
        end
        values
    end
    for (precision, spheroid) in ((:quad, :prolate), (:double, :oblate))
        T = precision===:quad ? BigFloat : Float64
        tolerance = precision===:quad ? big"2e-27" : 3e-12
        points = Complex{T}[0.3 + 0.2im, -0.4 - 0.3im, 1.4 + 0.2im]
        for c in (T(1.25), complex(T(1.25), T(0.125))), normalize in (false, true)

            coefficients = dmn(1, 2, c; spheroid, precision, normalize)
            basis = [ferrers(1, coefficients.degrees, z) for z in points]
            actual = smn(1, 2, c, points; spheroid, precision, normalize)
            @test actual.value ≈ [sum(coefficients.coefficients .* p) for p in basis] rtol=tolerance
            tangent = jacobian_smn(1, 2, c, points; spheroid, precision, normalize)
            @test (c isa Real ? tangent.dvalue_dc : tangent.dvalue_dcreal) ≈
                  [sum(coefficients.dcoefficients_dc .* p) for p in basis] rtol=tolerance
            @test actual.value isa Vector{Complex{T}}
        end
        z = first(points)
        q = smn(0, 0, 0, z; spheroid, precision, kind = 2,
            derivatives = 4, logderivative = true)
        @test only(q.value) ≈ atanh(z) rtol=tolerance
        @test only(q.derivative) ≈ inv(1-z^2) rtol=tolerance
        @test only(q.second_derivative) ≈ 2z/(1-z^2)^2 rtol=tolerance
        @test only(q.third_derivative) ≈ 2(1+3z^2)/(1-z^2)^3 rtol=tolerance
        @test only(q.fourth_derivative) ≈ 24z*(1+z^2)/(1-z^2)^4 rtol=tolerance
        @test only(q.logderivative) ≈ inv((1-z^2)*atanh(z)) rtol=tolerance
        p = smn(2, 3, 0, z; spheroid, precision, derivatives = 4)
        @test only(p.value) ≈ 15z*(1-z^2) rtol=tolerance
        @test only(p.derivative) ≈ 15-45z^2 rtol=tolerance
        @test only(p.second_derivative) ≈ -90z rtol=tolerance
        @test only(p.third_derivative) ≈ -90 rtol=tolerance
        @test abs(only(p.fourth_derivative)) < tolerance
        p1 = smn(1, 1, 0, z; spheroid, precision)
        @test only(p1.value) ≈ -sqrt(1-z)*sqrt(1+z) rtol=tolerance
        @test only(p1.derivative) ≈ z/(sqrt(1-z)*sqrt(1+z)) rtol=tolerance
        for kind in 1:2
            real_result = smn(
                1, 2, T(1.25), T[0.3]; spheroid, precision, kind, derivatives = 4)
            complex_result = smn(
                1, 2, T(1.25), Complex{T}[0.3]; spheroid, precision, kind, derivatives = 4)
            @test all(isapprox(a, b; rtol = tolerance)
            for (a, b) in zip(real_result, complex_result))
            @test smn(1, 2, T(1.25), conj.(points); spheroid, precision, kind).value ≈
                  conj.(smn(1, 2, T(1.25), points; spheroid, precision, kind).value) rtol=tolerance
        end
        ranged = smn(
            0, 0:2, T(1.25), z; spheroid, precision, scaled = true, derivatives = 3)
        @test size(ranged.value.mantissa) == (1, 3)
        for n in 0:2
            single = smn(
                0, n, T(1.25), z; spheroid, precision, scaled = true, derivatives = 3)
            @test ranged.third_derivative.mantissa[:, n + 1] ==
                  single.third_derivative.mantissa
            @test ranged.third_derivative.exponent[:, n + 1] ==
                  single.third_derivative.exponent
        end
        @test only(accuracy(0, 0, 1, [z]; spheroid, precision, target = :angular)) >=
              (precision===:quad ? 28 : 12)
    end
    for kind in 1:2
        z = ComplexF64[-1, 1, 0.25 + 0.125im]
        actual = smn(1, 2, 1.25+0.125im, z; kind, derivatives = 4,
            scaled = true, logderivative = true)
        reference = smn(1, 2, 1.25+0.125im, [-1, 1]; kind, derivatives = 4,
            scaled = true, logderivative = true)
        @test isequal(actual.derivative.mantissa[1:2], reference.derivative.mantissa)
        @test isequal(actual.fourth_derivative.mantissa[1:2], reference.fourth_derivative.mantissa)
        @test isequal(actual.logderivative[1:2], reference.logderivative)
    end
    @test smn(0, 1, 0, Number[0, 0.25im]).value == [0, 0.25im]
    @test isnan(only(smn(0, 1, 0, 0im; logderivative = true).logderivative))
    @test_throws DomainError smn(0, 0, 1, 2+0im)
    @test_throws DomainError smn(0, 0, 1, -2+0im)
    @test_throws ArgumentError smn(0, 0, 1, ComplexF64[])
    @test_throws ArgumentError smn(0, 0, 1, Inf+im)
    @test_throws ArgumentError smn(0, 0, 1, 0.3im; kind = 2, normalize = true)
    @test_throws ArgumentError smn(0, 0, 1, 0.3im; derivatives = 5)
end

@testset "Public expansion coefficients" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-27" : 2e-12
        x = T[-0.4, 0.25, 0.7]
        for c in (T(1.25), complex(T(1.25), T(0.125)), complex(zero(T), T(1.25))),
            normalize in (false, true)

            data = dmn(1, 2, c; spheroid, precision, normalize)
            @test data.converged && data.terms == length(data.degrees)
            @test all(iseven, data.degrees) && first(data.degrees) == 2
            basis = [smn(1, l, zero(T), x; precision).value for l in data.degrees]
            @test sum(d .* p for (d, p) in zip(data.coefficients, basis)) ≈
                  smn(1, 2, c, x; spheroid, precision, normalize).value rtol=tolerance
            tangent = jacobian_smn(1, 2, c, x; spheroid, precision, normalize)
            reference = c isa Real ? tangent.dvalue_dc : tangent.dvalue_dcreal
            @test sum(d .* p for (d, p) in zip(data.dcoefficients_dc, basis)) ≈ reference rtol=tolerance
            # The bilinear norm follows from Ferrers orthogonality.
            norm2(l) = T(2)*prod(T(k) for k in l:(l + 1))/(2l+1)
            @test sum(d^2*norm2(l) for (d, l) in zip(data.coefficients, data.degrees)) ≈
                  (normalize ? one(T) : norm2(2)) rtol=tolerance
            saved = copy(data.coefficients)
            fill!(data.coefficients, 0)
            fill!(data.dcoefficients_dc, 0)
            fill!(data.degrees, 0)
            @test dmn(1, 2, c; spheroid, precision, normalize).coefficients == saved
        end
        for (m, n) in ((0, 0), (0, 4), (1, 1), (2, 3)), c in (zero(T), complex(zero(T)))

            data = dmn(m, n, c; spheroid, precision)
            @test data.coefficients == [l == n ? 1 : 0 for l in data.degrees]
            @test all(iszero, data.dcoefficients_dc)
            @test data.eigenvalue == n*(n+1)
            @test amn(m, n, c; spheroid, precision) == 1
            @test amn(-m, n, c; spheroid, precision) == 1
            @test kmn(m, n, c; spheroid, precision) == (n == 0 ? 1 : 0)
        end
        c = precision === :quad ? big"1e-30" : 1e-8
        data = dmn(0, 0, c; spheroid, precision)
        sigma = spheroid === :prolate ? 1 : -1
        @test data.coefficients[2]/c^2 ≈ -T(sigma)/9 rtol=tolerance
        @test kmn(0, 1, c; spheroid, precision)/c ≈
              (spheroid === :prolate ? one(T) : complex(zero(T), one(T)))/3 rtol=tolerance
    end
    for f in (dmn, kmn, amn)
        @test_throws ArgumentError f(0, 0, 1; rtol = 0)
        @test_throws ArgumentError f(0, 4, 1; max_terms = 1)
        @test_throws r"did not converge" f(0, 2, 1; max_terms = 4)
        @test_throws r"spheroid must be" f(0, 0, 1; spheroid = :invalid)
        @test_throws r"precision must be" f(0, 0, 1; precision = :invalid)
        @test_throws r"0 <= m <= n" f(2, 1, 1)
    end
end

@testset "Joining and radial factor connection identities" begin
    SW = SpheroidalWaves
    plan = SW._coefficient_plan(0, 0, 0.0)
    @test_throws r"denominator is singular" SW._joining_factor_value(
        (; plan..., v = zero(plan.v)), SW._SWFloat(1))
    @test_throws r"origin normalization" SW._joining_factor_value(plan, SW._SWFloat(0))
    # Exterior Legendre polynomials, independently of the package's Ferrers
    # evaluator. z*sqrt(1-z^-2) specifies the exterior branch for odd order.
    function exterior(m, degrees, z)
        previous = zero(z)
        current = prod(oftype(z, k) for k in 1:2:(2m - 1); init = one(z))*(z*sqrt(1-inv(z^2)))^m
        values = typeof(z)[]
        for l in m:last(degrees)
            l in degrees && push!(values, current)
            previous, current = current, ((2l+1)*z*current-(l+m)*previous)/(l-m+1)
        end
        values
    end
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"2e-26" : 3e-11
        for (m, n) in ((0, 0), (0, 1), (1, 1), (1, 2), (2, 4)),
            c in (T(1.25), complex(T(1.25), T(0.125)))

            data = dmn(m, n, c; spheroid, precision)
            K = kmn(m, n, c; spheroid, precision)
            for x in T[1.25, 1.75]
                z = spheroid === :prolate ? x : -im*x
                continued = sum(data.coefficients .* exterior(m, data.degrees, z))
                @test K*continued ≈ only(rmn(m, n, c, x; spheroid, precision).value) rtol=tolerance
            end
            norm = sqrt(T(2)*prod(T(k) for k in (n - m + 1):(n + m); init = one(T))/(2n+1))
            @test kmn(m, n, c; spheroid, precision, normalize = true) ≈ K*norm rtol=tolerance
            Aplus = amn(m, n, c; spheroid, precision)
            Aminus = amn(-m, n, c; spheroid, precision)
            w(l) = prod(T(k) for k in (l - m + 1):(l + m); init = one(T))
            @test Aminus ≈
                  sum(d*w(l) for (d, l) in zip(data.coefficients, data.degrees))/w(n) rtol=tolerance
            p = smn(m, n, c, T(0.3); spheroid, precision)
            q = smn(m, n, c, T(0.3); spheroid, precision, kind = 2)
            @test (1-T(0.3)^2)*only(p.value .* q.derivative-p.derivative .* q.value) ≈
                  w(n)*Aplus*Aminus rtol=tolerance
            @test K isa (spheroid === :oblate || c isa Complex ? Complex{T} : T)
            @test Aplus isa (c isa Complex ? Complex{T} : T)
        end
    end
end

@testset "Native spherical and imaginary-parameter angular limits" begin
    SW = SpheroidalWaves
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        prefix = spheroid === :prolate ? :cprolate : :coblate
        opposite = spheroid === :prolate ? :oblate : :prolate
        for c in (0im, 1im)
            raw = SW._call_complex_smn_raw(prefix, 0, 1, c, [0.3]; precision)
            reference = smn(0, 1, abs(imag(c)), 0.3; spheroid = opposite, precision)
            @test raw.value ≈ reference.value
            @test raw.derivative ≈ reference.derivative
        end
        # Coefficient reconstruction must choose the same phase on the
        # imaginary axis as the real problem with opposite geometry.
        expansion = SW._angular_coefficients(1, 2, 1im; precision, spheroid)
        reconstructed = sum(d .* smn(1, l, 0, 0.3; precision).value
        for (d, l) in zip(expansion.coefficients, expansion.degrees))
        @test reconstructed ≈ smn(1, 2, 1im, 0.3; precision, spheroid).value rtol=(precision ===
                                                                                   :quad ?
                                                                                   big"1e-28" :
                                                                                   1e-12)
    end
end

@testset "Continuation endpoint anchors and unresolved paths" begin
    SW = SpheroidalWaves
    # New trackers reuse native data, but their profiles remain independent.
    first = SW._angular_phase_evaluator(:cprolate, 0, 1, :double)(1.25+0.1im)
    expected = copy(first.profile)
    fill!(first.profile, NaN)
    second = SW._angular_phase_evaluator(:cprolate, 0, 1, :double)(1.25+0.1im)
    @test second.profile == expected

    for precision in (:double, :quad)
        evaluate = SW._angular_phase_evaluator(
            :cprolate, 0, 1, precision; endpoint_anchor = true)
        state = evaluate(1.25+0.1im)
        @test state.endpoint_anchor && 0 < state.endpoint_point < 1
        @test state.anchor == state.profile[end - 1]
        @test length(state.profile) == 10
        @test_throws r"origin anchor" SW._require_angular_anchor(state)
    end
    state = (m = 0, n = 0, lambda = 0.0, anchor = 1.0, profile = [1.0, 0.0], sign = 1,
        endpoint_anchor = false, endpoint_point = 0.0)
    constant(c, n = 0) = (; state..., n, lambda = Float64(n*(n+1)))
    @test SW._transport_angular_phase(
        constant, 1.0, state, nextfloat(1.0); branch_window = 0).anchor == 1
    @test SW._transport_angular_phase(constant, 0.0, state, 0.1; branch_window = 1).n == 0
    @test_throws r"must be finite" SW._transport_angular_phase(constant, 0.0, state, Inf)
    # An unresolved jump must fail even when no representable midpoint remains.
    jump(c, n = 0) = (; state..., profile = [0.0, 1.0])
    @test_throws r"could not resolve" SW._transport_angular_phase(
        jump, 1.0, state, nextfloat(1.0); branch_window = 0)
    @test_throws ArgumentError SW._angular_coefficients(0, 4, 1; max_terms = 1)
    limited = SW._coefficient_plan(0, 2, 1; max_terms = 4)
    @test !limited.converged && limited.terms == 4
    plan = SW._coefficient_plan(0, 1, 1.25)
    for value in (big"0.0", BigFloat(NaN))
        unresolved = (; plan..., v = fill(value, length(plan.v)),
            phase = Ref{Union{Nothing, Int}}(nothing))
        @test_throws r"phase cannot be resolved" SW._coefficient_phase(unresolved, :double)
        @test unresolved.phase[] === nothing
    end
end

@testset "Angular phase across native backends and degree ranges" begin
    for precision in (:double, :quad)
        lib = SpheroidalWaves.backend_library(; precision)
        if lib === nothing || !isfile(lib)
            @info "Skipping native angular phase tests: backend unavailable." precision
            continue
        end
        # Independent reference for m=1, n=1, c=0.3, eta=0.
        value = only(smn(1, 1, 0.3, 0.0; precision).value)
        @test value ≈ -1.001795214382128 rtol=1e-14

        for spheroid in (:prolate, :oblate), complex_c in (false, true),
            normalize in (false, true)
            c = complex_c ? 0.3 + 0.0im : 0.3
            eta = [-0.25, 0.0, 0.25]
            for m in (1, 2)
                block = smn(m, m:(m + 3), c, eta; precision, spheroid, normalize)
                for (j, n) in enumerate(m:(m + 3))
                    single = smn(m, n, c, eta; precision, spheroid, normalize)
                    @test block.value[:, j]≈single.value rtol=1e-12 atol=1e-14
                    @test block.derivative[:, j]≈single.derivative rtol=1e-12 atol=1e-14
                    # DLMF 30.4(i), the phase specified after Eq. 30.4.1.
                    if iseven(n - m)
                        @test sign(real(block.value[2, j])) == (-1)^((n + m) ÷ 2)
                    else
                        @test sign(real(block.derivative[2, j])) == (-1)^((n + m - 1) ÷ 2)
                    end
                end
            end
            # Continuity to the closed-form spherical limit, independently of
            # whether the backend uses real or complex arithmetic.
            small_c = complex_c ? 1e-3 + 0.0im : 1e-3
            near = smn(1, 2, small_c, eta; precision, spheroid, normalize)
            exact = smn(1, 2, zero(small_c), eta; precision, spheroid, normalize)
            @test near.value≈exact.value rtol=1e-6 atol=1e-12
            @test near.derivative≈exact.derivative rtol=1e-6 atol=1e-12
        end

        # Approach zero through both sides of arg(c)=pi/4, where the old
        # complex oblate normalization selected opposite signs for n=m+2.
        for spheroid in (:prolate, :oblate), normalize in (false, true),
            c in (1e-3 + 1e-4im, 1e-4 + 1e-3im), n in (3, 4)
            eta = [-0.25, 0.25] # Phase probes are independent of this batch.
            near = smn(1, n, c, eta; precision, spheroid, normalize)
            exact = smn(1, n, zero(c), eta; precision, spheroid, normalize)
            @test near.value≈exact.value rtol=1e-6 atol=1e-12
            @test near.derivative≈exact.derivative rtol=1e-6 atol=1e-12
            with_anchor = smn(1, n, c, [0.0; eta]; precision, spheroid, normalize)
            @test near.value ≈ with_anchor.value[2:end] rtol=1e-13
            @test near.derivative ≈ with_anchor.derivative[2:end] rtol=1e-13
        end

        # The phase must propagate to the parameter Jacobian: the negative
        # S_11(c, 0) decreases locally around this reference point.
        jac = jacobian_smn(1, 1, 0.3, [0.0]; precision, h = 1e-4, adaptive = false)
        @test only(jac.dvalue_dc) < 0
        @test only(jac.dderivative_dc) == 0
    end
end
