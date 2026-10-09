using SpheroidalWaves, Test

@testset "One-sided complex endpoint derivatives" begin
    fields = (
        :value, :derivative, :second_derivative, :third_derivative, :fourth_derivative)
    for precision in (:double, :quad)
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"2e-27" : 2e-12
        for c in (zero(Complex{T}), complex(T(1.25), T(0.125)))
            for spheroid in (:prolate, :oblate), m in 0:9

                n = m+1
                r = smn(m, n, c, T[-1, 1]; spheroid, precision, derivatives = 4)
                close = smn(m, n, c, [-1+big"1e-12", 1-big"1e-12"];
                    spheroid, precision, derivatives = 4)
                for (k, field) in enumerate(fields)
                    values = getproperty(r, field)
                    @test !any(isnan, values)
                    if m>2(k-1)
                        @test all(iszero, values)
                    elseif isodd(m)
                        @test all(isinf, values)
                        @test signbit.(real.(values)) ==
                              signbit.(real.(getproperty(close, field)))
                    else
                        @test values≈getproperty(close, field) rtol=1e-8 atol=1e-8
                    end
                end
            end
        end
        # Exact endpoint derivatives of associated Legendre polynomials.
        for (m, n, third, fourth) in ((0, 4, [-105, 105], [105, 105]), (
            2, 3, [-90, -90], [0, 0]),
            (4, 4, [-2520, 2520], [2520, 2520]), (
            6, 6, [498960, -498960], [-2993760, -2993760]),
            (8, 8, [0, 0], fill(384*2027025, 2)))
            r = smn(m, n, zero(Complex{T}), T[-1, 1]; precision, derivatives = 4)
            @test r.third_derivative≈third rtol=tolerance atol=tolerance
            @test r.fourth_derivative≈fourth rtol=tolerance atol=tolerance
        end
        c = complex(T(1.25), T(0.125))
        for m in 0:4, kind in 1:4

            r = rmn(
                m, m+1, c, T[1, 2]; precision, kind, derivatives = 4, logderivative = true)
            nearby = rmn(m, m+1, c, [1+big"1e-10"]; precision, kind, derivatives = 4)
            for field in fields
                endpoint = first(getproperty(r, field))
                @test kind==1 ? !isnan(endpoint) : isnan(endpoint)
                if kind==1 && isfinite(endpoint)
                    @test endpoint≈only(getproperty(nearby, field)) rtol=1e-6 atol=1e-5
                end
            end
            @test r.fourth_derivative[2] ≈
                  only(rmn(m, m+1, c, T[2]; precision, kind, derivatives = 4).fourth_derivative) rtol=tolerance
            @test first(accuracy(m, m+1, c, T[1, 2]; precision, kind)) == -1
            kind==1 && m>0 && @test first(r.logderivative) == Inf
        end
        for kind in (1, 2), spheroid in (:prolate, :oblate)

            r = smn(1, 2, c, T[-1, 1]; precision, spheroid, kind,
                derivatives = 4, scaled = true, logderivative = true)
            @test all(isinf, r.third_derivative.mantissa)
            @test all(isinf, r.fourth_derivative.mantissa)
            @test r.logderivative == (kind==1 ? [Inf, -Inf] : [-Inf, Inf])
        end
        for m in (0, 2)
            imaginary = rmn(
                m, m+1, T(1.25)*im, [big"1", 1+big"1e-20"]; precision, derivatives = 4)
            @test all(field -> !any(isnan, getproperty(imaginary, field)), fields)
            @test imaginary.second_derivative[1]≈imaginary.second_derivative[2] rtol=1e-15 atol=1e-15
            @test imaginary.fourth_derivative[1]≈imaginary.fourth_derivative[2] rtol=1e-15 atol=1e-15
        end
        origin = rmn(0, 1, c, T[0]; spheroid = :oblate, precision, logderivative = true)
        @test only(origin.logderivative) == Inf
        batch = rmn(0, 0:1, c, T[1]; precision, derivatives = 4)
        @test all(isfinite, batch.fourth_derivative)
    end
end

@testset "Higher coordinate derivatives" begin
    offsets = collect(-4:4)
    matrix = Rational{BigInt}[big(k)^p for p in 0:8, k in offsets]
    weights(order) = matrix \
                     Rational{BigInt}[p==order ? factorial(big(order)) : 0 for p in 0:8]
    w3, w4 = weights(3), weights(4)
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T = precision === :quad ? BigFloat : Float64
        tol = precision === :quad ? big"2e-27" : 2e-12
        for x in (T(0.25), big"1"-big"1e-20", -big"1"+big"1e-20")
            polynomial = smn(0, 4, 0, [x]; precision, spheroid, derivatives = 4)
            @test only(polynomial.third_derivative) ≈ 105x rtol=tol
            @test only(polynomial.fourth_derivative) ≈ 105 rtol=tol
            associated = smn(2, 3, 0, [x]; precision, spheroid, derivatives = 4)
            @test only(associated.third_derivative) ≈ -90 rtol=tol
            @test abs(only(associated.fourth_derivative)) < tol
            root = smn(1, 1, 0, [x]; precision, spheroid, derivatives = 4)
            u = (1-x)*(1+x)
            @test only(root.third_derivative) ≈ 3x/u^(5//2) rtol=tol
            @test only(root.fourth_derivative) ≈ 3(1+4x^2)/u^(7//2) rtol=tol
        end
        for c in (T(1.25), complex(T(1.25), T(0.125)))
            for (f, x, kinds) in ((smn, T(0.3), (1, 2)), (rmn, T(2), (1, 2, 3, 4)))
                for kind in kinds
                    result = f(1, 2, c, [x]; spheroid, precision, kind, derivatives = 4)
                    # Independent nine-point differences of function values.
                    wide = c isa Real ? BigFloat(c) : Complex{BigFloat}(c)
                    h = big"0.0005"
                    grid = BigFloat(x) .+ h .* offsets
                    values = f(1, 2, wide, grid; spheroid, precision = :quad, kind).value
                    values .-= values[5]
                    @test only(result.third_derivative) ≈ sum(BigFloat.(w3) .* values)/h^3 rtol=2e-13
                    @test only(result.fourth_derivative) ≈ sum(BigFloat.(w4) .* values)/h^4 rtol=2e-13
                    @test eltype(result.fourth_derivative) ==
                          (f===rmn || c isa Complex ? Complex{T} : T)
                    scaled = f(1, 2, c, [x]; spheroid, precision, kind,
                        derivatives = 3, scaled = true, logderivative = true)
                    @test !hasproperty(scaled, :fourth_derivative)
                    @test scaled.third_derivative.mantissa .*
                          T(10) .^ scaled.third_derivative.exponent ≈
                          result.third_derivative rtol=tol
                    @test scaled.logderivative ≈ result.derivative ./ result.value rtol=tol
                end
            end
        end
        for f in (smn, rmn)
            x = f===smn ? T(0.3) : T(2)
            batch = f(0, 1:2, T(1.25), x; spheroid, precision, derivatives = 4)
            for (i, n) in enumerate(1:2)
                @test batch.fourth_derivative[:, i] ==
                      f(0, n, T(1.25), [x]; spheroid, precision, derivatives = 4).fourth_derivative
            end
            @test f(0, 1, T(1.25), [x]; spheroid, precision, derivatives = 2) ==
                  f(0, 1, T(1.25), [x]; spheroid, precision, second_derivative = true)
            @test f(0, 1:2, T(1.25), [x]; spheroid, precision, derivatives = 2) ==
                  f(0, 1:2, T(1.25), [x]; spheroid, precision, second_derivative = true)
        end
        for m in (0, 1, 2)
            x = big"1"+big"1e-20"
            result = rmn(
                m, m+1, T(1.25), [x]; spheroid = :prolate, precision, derivatives = 4)
            @test all(isfinite, result.third_derivative) &&
                  all(isfinite, result.fourth_derivative)
        end
        origin = rmn(
            0, 1, T(1.25), [zero(T)]; spheroid = :oblate, precision, derivatives = 4)
        @test only(origin.fourth_derivative) == 0
        # The regular exponent m/2 fixes these near-boundary derivative ratios.
        c, distance = complex(T(1.25), T(0.125)), big"1e-20"
        for (f, x, sign) in ((smn, 1-distance, 1), (smn, -1+distance, -1), (
            rmn, 1+distance, -1))
            r = f(1, 2, c, [x]; precision, derivatives = 4)
            @test only(r.fourth_derivative ./ r.third_derivative)*distance ≈ sign*5//2 rtol=1e-14
        end
    end
    for f in (smn, rmn), n in (1, 1:2)

        @test_throws ArgumentError f(0, n, 1, [1]; derivatives = 5)
        @test all(isfinite, f(0, n, 1, [1]; derivatives = 3).third_derivative)
    end
end

@testset "Independent coordinate residual stencils" begin
    SW = SpheroidalWaves
    for precision in (:double, :quad)
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? T(1e-20) : T(1e-6)
        # P_3 and its derivatives give exact polynomial stencil references.
        x = T[0.2, 0.6]
        diagnostic = SW._spheroidal_residual(0, 3, zero(T), x; precision)
        @test diagnostic.derivative ≈ (15x .^ 2 .- 3) ./ 2 rtol=tolerance
        @test diagnostic.second_derivative ≈ 15x rtol=tolerance
        @test maximum(diagnostic.relative_residual) < tolerance
        constant = SW._spheroidal_residual(0, 0, zero(T), T(0.3); precision)
        @test only(constant.residual) == only(constant.relative_residual) == 0
        for spheroid in (:prolate, :oblate)
            # Oblate x=0 requires a forward stencil, while x=2 is centered.
            points = spheroid === :oblate ? T[0, 2] : T[2]
            diagnostic = SW._spheroidal_residual(0, 0, T(1), points;
                target = :radial, spheroid, precision)
            direct = rmn(0, 0, T(1), points; spheroid, precision, second_derivative = true)
            @test diagnostic.derivative≈direct.derivative rtol=tolerance atol=tolerance
            @test diagnostic.second_derivative≈direct.second_derivative rtol=tolerance atol=tolerance
            @test maximum(diagnostic.relative_residual) < 10tolerance
            @test all(diagnostic.step .> 0)
        end
    end
    @test_throws ArgumentError SW._spheroidal_residual(0, 1, 1, 0.3; target = :bad)
    @test_throws DomainError SW._spheroidal_residual(0, 1, 1, 1)
    @test_throws DomainError SW._spheroidal_residual(0, 1, 1, 1; target = :radial)
    for h in (0, -1, Inf, NaN)
        @test_throws ArgumentError SW._spheroidal_residual(0, 1, 1, 0.3; h)
    end
    @test_throws r"below coordinate resolution" SW._spheroidal_residual(
        0, 1, 1, 0.3; h = eps(0.3)/8)
end

@testset "Taylor steps reject unconverged series" begin
    SW = SpheroidalWaves
    # Q_0(x)=atanh(x) has convergence radius one at x=0. A step to x=2
    # must fail its tail check rather than return a truncated polynomial.
    @test_throws r"did not converge" SW._qs_step(0, big"0.0", big"0.0", big"0.0",
        big"0.0", big"1.0", big"2.0")
    # At c=lambda=m=0 the radial equation admits log((x+1)/(x-1)).
    # Its Taylor series at x=2 cannot cross the singularity x=1.
    plan = (spheroid = :prolate, c = big"0", lambda = big"0", dlambda = big"0", m = 0)
    @test_throws r"did not converge" SW._radial_sensitivity_step(plan, big"2",
        (log(big"3"), -big"2"/3, big"0", big"0"), -big"2")
end

@testset "Regular endpoint second derivatives" begin
    SW = SpheroidalWaves
    setprecision(BigFloat, 256) do
        # P_3^2(x)=15x(1-x^2), so its endpoint slopes are -30 and its
        # second derivatives are +90 at -1 and -90 at +1.
        spherical = SW._evaluate_coefficient_vector((m = 2, degrees = [3]),
            [sqrt(big"240.0"/7)], [-1, 1]; second_derivative = true)
        @test spherical.value == [0, 0]
        @test spherical.derivative ≈ [-30, -30] rtol=big"1e-70"
        @test spherical.second_derivative ≈ [90, -90] rtol=big"1e-70"
        for precision in (:double, :quad), spheroid in (:prolate, :oblate)
            # At m=2, the endpoint limit of the angular equation gives
            # S''(x)=x*(lambda-sigma*c^2-3)*S'(x)/3 for x=+/-1.
            # c=5 selects coefficient evaluation through the public API.
            result = smn(2, 3, 5, [-1, 1]; precision, spheroid, second_derivative = true)
            lambda = eigenvalue(2, 3, 5; precision, spheroid)
            sigma = spheroid === :prolate ? 1 : -1
            expected = [-1, 1] .* (lambda-sigma*25-3) .* result.derivative ./ 3
            @test result.second_derivative ≈ expected rtol=(precision === :quad ?
                                                            big"1e-27" : 1e-12)
        end
        for precision in (:double, :quad), m in 0:5

            n, c = m+1, big"1.25"
            points = BigFloat[-1, -1 + big"1e-8", 1 - big"1e-8", 1]
            # A differentiated Legendre expansion provides an independent
            # reference for the native endpoint Taylor reconstruction.
            plan = SW._coefficient_plan(m, n, c; precision = :quad)
            scale = SW._coefficient_phase(plan, :quad)*sqrt(SW._ferrers_norm2(m, n, BigFloat))
            reference = SW._evaluate_coefficient_vector(plan, scale .* plan.v, points;
                second_derivative = true)
            result = smn(
                m, n, c, points; precision, second_derivative = true, scaled = true)
            for field in (:value, :derivative, :second_derivative)
                encoded = getproperty(result, field)
                actual = encoded.mantissa .* BigFloat(10) .^ encoded.exponent
                expected = getproperty(reference, field)
                @testset "$precision m=$m $field" begin
                    for i in eachindex(points)
                        @test actual[i] ≈ expected[i] rtol=(precision === :quad ?
                                                            big"1e-25" : 1e-11)
                    end
                end
            end
            # The regular radial solution has the same endpoint power.
            radial = rmn(
                m, n, c, BigFloat[1, 1 + big"1e-8"]; precision, second_derivative = true)
            if isodd(m) && m < 4
                @test isinf(real(radial.second_derivative[1]))
            elseif m > 4
                @test radial.second_derivative[1] == 0
            else
                @test isfinite(radial.second_derivative[1])
                @test radial.second_derivative[2] ≈ radial.second_derivative[1] rtol=1e-5
            end
        end
    end
end

@testset "Large bandwidth tolerance and second-kind residuals" begin
    # At c=375, exp(-2c) underflows in Float64 even though the requested
    # coefficient tolerance is meant to be formed at the guarded precision.
    # Wider-input reference values are inline to avoid rebuilding the large
    # quad expansions during every normal test run.
    for (spheroid, value, ve, derivative, de) in (
        (:prolate, big"1.5256488534902150635105594465576701032672", -7,
        big"-1.7954837073168552408315797203151455810388", -5),
        (:oblate, big"4.1896917085417494071462469474848952172303", -113,
        big"1.5679082330112565415126277769813336590425", -110))
        ordinary=smn(0, 0, 375.0, big".3"; spheroid, precision = :double, scaled = true)
        @test only(ordinary.value.mantissa)≈value rtol=1e-13
        @test only(ordinary.value.exponent)==ve
        @test only(ordinary.derivative.mantissa)≈derivative rtol=1e-13
        @test only(ordinary.derivative.exponent)==de
    end
    for spheroid in (:prolate, :oblate), kind in (1, 2)

        diagnostic=SpheroidalWaves._spheroidal_residual(0, 0, big"1", big".3";
            spheroid, kind, precision = :quad, h = big"1e-5")
        actual=smn(0, 0, big"1", big".3"; spheroid, kind,
            precision = :quad, second_derivative = true)
        @test diagnostic.derivative≈actual.derivative rtol=big"1e-20"
        @test diagnostic.second_derivative≈actual.second_derivative rtol=big"1e-20"
        @test maximum(diagnostic.relative_residual)<big"1e-20"
    end
end
@testset "Shared complex-coordinate continuation" begin
    SW = SpheroidalWaves
    F = SW._SWFloat
    # Rounded slopes can coincide for distinct rays. Their grouping keys must
    # retain enough precision to distinguish ratios of the input mantissas.
    SW._with_swprecision(80) do
        initial = SW._coordinate_initial_data(
            0, 0, F(0), :prolate, :quad, :angular, 1, false, :standard)
        u = eps(F)
        points = [complex(1+u, F(1)), complex(F(1), 1-u)]
        @test imag(points[1])/real(points[1]) == imag(points[2])/real(points[2])
        step_count = Ref(0)
        for z in points
            SW._coordinate_continuation(
                initial.equation, initial.state, initial.anchor, [z]; step_count)
        end
        @test SW._coordinate_batch(initial, points, :angular, :prolate).steps ==
              step_count[]
    end
    SW._with_swprecision(384) do
        for target in (:angular, :radial), spheroid in (:prolate, :oblate)

            initial = SW._coordinate_initial_data(
                1, 4, F(1.25), spheroid, :quad, target, 1, false, :standard)
            points = target===:angular ? [complex(F(k)/16, F(k)/32) for k in 1:12] :
                     [complex(F(2), F(k)/8) for k in 1:12]
            batch = SW._coordinate_batch(initial, points, target, spheroid)
            step_count = Ref(0)
            separate = [SW._coordinate_continuation(
                            initial.equation, initial.state, initial.anchor,
                            SW._coordinate_path(z, target, spheroid); step_count)
                        for z in points]
            @test batch.steps < step_count[]/2
            @test all(all(isapprox(a, b; rtol = F("1e-28"), atol = F("1e-50"))
                      for (a, b) in zip(x, y)) for (x, y) in zip(batch.states, separate))
            @test batch.states ==
                  reverse(SW._coordinate_batch(initial, reverse(points), target, spheroid).states)
            duplicates = SW._coordinate_batch(initial, [points; points], target, spheroid)
            @test duplicates.states == [batch.states; batch.states]
            @test duplicates.steps == batch.steps
        end
    end
    for spheroid in (:prolate, :oblate)
        points = Complex{BigFloat}[
            -2, -2 + 0.5im, 0.3 + 0.5im, 2 + 0.5im, 2 - 0.5im, 2 + 0.5im]
        batch = rmn(
            0, 0, 0, points; spheroid, precision = :quad, kind = 2, normalization = :static)
        expected = spheroid===:prolate ? -atanh.(inv.(points)) : atan.(points) .- big(pi)/2
        @test batch.value ≈ expected rtol=big"1e-28"
        @test batch.value ==
              reverse(rmn(0, 0, 0, reverse(points); spheroid, precision = :quad,
            kind = 2, normalization = :static).value)
    end
end
