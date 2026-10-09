using SpheroidalWaves, Test

@testset "Static radial normalization" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"2e-28" : 2e-12
        sigma = spheroid === :prolate ? 1 : -1
        points = spheroid === :prolate ? T[1.25, 2, 3] : T[0, 0.5, 2, 3]
        for m in 0:2, n in m:(m + 2), kind in 1:2
            limit = rmn(
                m, n, zero(T), points; spheroid, precision, kind, normalization = :static)
            small = rmn(
                m, n, T(1e-8), points; spheroid, precision, kind, normalization = :static)
            @test small.value≈limit.value rtol=1e-12 atol=1e-13
            @test small.derivative≈limit.derivative rtol=1e-12 atol=1e-13
            @test eltype(limit.value) == Complex{T}
            @test maximum(radial_wronskian(m, n, 0, points; spheroid, precision,
                normalization = :static, form = :error)) < tolerance
        end
        for x in points
            angle = spheroid === :prolate ? atanh(inv(x)) :
                    iszero(x) ? T(pi)/2 : atan(inv(x))
            u = x^2-sigma
            # Closed Laplace solutions, independent of the implemented series.
            for (m, n, v, d) in ((0, 0, -angle, 1/u), (
                0, 1, -3sigma*(x*angle-1), -3sigma*(angle-x/u)),
                (1, 1, 3sigma/2*sqrt(u)*(angle-x/u),
                3sigma/2*(x/sqrt(u)*(angle-x/u)+2sigma/u^(3//2))))
                r = rmn(
                    m, n, 0, [x]; spheroid, precision, kind = 2, normalization = :static)
                @test only(r.value)≈v rtol=tolerance atol=tolerance
                @test only(r.derivative)≈d rtol=tolerance atol=tolerance
            end
            @test only(rmn(0, 1, 0, [x]; spheroid, precision, normalization = :static).value) ≈
                  x/3 rtol=tolerance
        end
        for c in (zero(T), complex(T(1.25), T(0.125))), kind in 1:2

            result = rmn(
                1, 1:2, c, T[2]; spheroid, precision, kind, normalization = :static,
                derivatives = 4, scaled = true, logderivative = true)
            for (i, n) in enumerate(1:2)
                scalar = rmn(1, n, c, T[2]; spheroid, precision, kind,
                    normalization = :static, derivatives = 4, logderivative = true)
                for field in (:value, :derivative, :second_derivative, :third_derivative, :fourth_derivative)
                    scaled = getproperty(result, field)
                    @test scaled.mantissa[:, i] .* T(10) .^ scaled.exponent[:, i] ≈
                          getproperty(scalar, field) rtol=tolerance
                end
                @test result.logderivative[:, i] ≈ scalar.logderivative rtol=tolerance
            end
            @test maximum(radial_wronskian(1, 2, c, T[2]; spheroid, precision,
                normalization = :static, form = :error)) < tolerance
            @test only(radial_wronskian(
                1, 2, c, T[2]; spheroid, precision, normalization = :static)) ≈
                  T(1)/(4-sigma) rtol=tolerance
            @test only(radial_wronskian(1, 2, c, T[2]; spheroid, precision,
                normalization = :static, form = :normalized)) ≈ 1 rtol=tolerance
            @test accuracy(
                1, 2, c, T[2]; spheroid, precision, kind, normalization = :static) == [-1]
        end
        regular = rmn(1, 1, 0, T[1]; precision, normalization = :static,
            derivatives = 4, logderivative = true)
        @test only(regular.value) == 0
        @test real(only(regular.derivative)) == Inf
        @test only(regular.logderivative) == Inf
        singular = rmn(
            0, 0, 0, T[1, 2]; precision, kind = 2, normalization = :static, derivatives = 4)
        @test isnan(singular.value[1]) && isnan(singular.fourth_derivative[1])
        @test isfinite(singular.fourth_derivative[2])
        for c in (T(1.25), complex(T(1.25)), complex(zero(T), T(1.25)), complex(T(1.25), T(0.125)))
            endpoint = rmn(
                1, 2, c, T[1]; precision, normalization = :static, derivatives = 4)
            @test all(field -> !any(isnan, getproperty(endpoint, field)), propertynames(endpoint))
            @test isinf(only(endpoint.derivative))
        end
    end
    for n in (1, 1:2)
        @test_throws ArgumentError rmn(0, n, 0, [2]; normalization = :static, kind = 3)
        @test_throws ArgumentError rmn(0, n, 1, [2]; normalization = :bad)
    end
    @test_throws ArgumentError accuracy(
        0, 1, 0, [0.3]; target = :angular, normalization = :static)
    @test_throws r"did not converge" SpheroidalWaves._static_irregular_radial(
        0, 0, big"2.0", 1; max_terms = 1)
end

@testset "Mixed radial boundary and interior batches" begin
    SW = SpheroidalWaves
    for precision in (:double, :quad), kind in 2:4

        tolerance = precision === :quad ? big"1e-28" : 1e-12
        for scaled in (false, true), points in ([1.0], [1.0, 2.0])

            result = rmn(
                0, 1, 1.25, points; precision, kind, scaled, second_derivative = true)
            unpack(v) = scaled ? v.mantissa .* BigFloat(10) .^ v.exponent : v
            @test isnan(first(unpack(result.value)))
            @test isnan(first(unpack(result.derivative)))
            @test isnan(first(unpack(result.second_derivative)))
            if length(points) == 2
                interior = rmn(0, 1, 1.25, 2; precision, kind, second_derivative = true)
                for field in keys(interior)
                    @test unpack(getproperty(result, field))[2] ≈
                          only(getproperty(interior, field)) rtol=tolerance
                end
            end
        end
        if precision === :quad
            native = SW._call_real_rmn(
                :psms, 0, 1, 1.25, [1.0, 2.0]; precision, kind, with_accuracy = true)
            @test native.accuracy[1] == -1
            @test native.accuracy[2] == only(SW._call_real_rmn(
                :psms, 0, 1, 1.25, [2.0]; precision, kind, with_accuracy = true).accuracy)
        end
    end
end

@testset "Complex native and analytic radial transfers" begin
    SW = SpheroidalWaves
    for spheroid in (:prolate, :oblate), precision in (:double, :quad)

        prefix = spheroid === :prolate ? :cprolate : :coblate
        tolerance = precision === :quad ? big"1e-27" : 1e-11
        for c in (10im, 1.25+0.1im), kind in 3:4

            direct = rmn(0, 1, c, 2; spheroid, precision, kind)
            native = SW._call_complex_rmn_raw(prefix, 0, 1, c, [2]; precision, kind)
            @test native.value ≈ direct.value rtol=tolerance
            @test native.derivative ≈ direct.derivative rtol=tolerance
            if precision === :quad && !iszero(real(c))
                quad = SW._complex_quad_radial(prefix, 0, 1, c, [2], kind)
                @test quad.value ≈ direct.value rtol=tolerance
                @test quad.derivative ≈ direct.derivative rtol=tolerance
            end
            if iszero(real(c))
                scaled = SW._scaled_native_values(
                    0, 1, c, [2], spheroid, precision, :radial, kind)
                @test scaled.value ≈ direct.value rtol=tolerance
            end
        end
    end
    boundary = SW._scaled_native_values(
        0, 1, 1.25+0.1im, [1], :prolate, :double, :radial, 1)
    @test boundary.value ≈ rmn(0, 1, 1.25+0.1im, [1]).value
end

@testset "Radial functions exclude exact zero parameter" begin
    # Reject before calling a backend, for every radial kind and input family.
    # The former prolate m=0 fallback returned unrelated Legendre P/Q values.
    for spheroid in (:prolate, :oblate), precision in (:double, :quad),
        c in (0, 0.0, -0.0, 0.0 + 0.0im, big"0.0", complex(big"0.0", big"0.0"))
        for kind in 1:4
            @test_throws DomainError rmn(0, 1, c, [2.0]; spheroid, precision, kind)
            @test_throws DomainError accuracy(
                0, 1, c, [2.0]; spheroid, precision, kind, target = :radial)
        end
        @test_throws DomainError rmn(1, 1, c, [2.0]; spheroid, precision)
        for (m, n) in ((0, 0), (0, 2), (1, 2))
            @test_throws DomainError rmn(m, n, c, [2.0]; spheroid, precision)
        end
        @test_throws DomainError rmn(0, 0:2, c, [2.0]; spheroid, precision)
        @test_throws DomainError radial_wronskian(0, 1, c, [2.0]; spheroid, precision)
        # A finite-difference derivative must reject an undefined center too.
        @test_throws DomainError jacobian_rmn(0, 1, c, [2.0]; spheroid, precision)
    end

    # A singular-coordinate error must not conceal the parameter-domain error.
    err = try
        rmn(0, 0, 0.0, [1.0])
    catch e
        e
    end
    @test err isa DomainError
    @test occursin("c != 0", sprint(showerror, err))
end

@testset "Small nonzero radial parameters remain supported" begin
    # Independent signed reference values, with decimal inputs preserved
    # for the high-precision second-kind comparisons.
    for precision in (:double, :quad)
        lib = SpheroidalWaves.backend_library(; precision)
        if lib === nothing || !isfile(lib)
            @info "Skipping nonzero radial references: backend unavailable." precision
            continue
        end
        first = rmn(0, 1, big"0.1", [2]; precision, kind = 1)
        @test only(first.value) ≈ 0.06642698122211762 rtol=5e-15
        tolerance = precision === :quad ? big"1e-27" : big"5e-15"
        for (c, expected) in [
            (big"0.01", big"-2958.897611846449350102731791001293876590218151081837314"),
            (big"0.001", big"-295837.3950030077336871204221250744159129119120049634542")]
            second = rmn(0, 1, c, [2]; precision, kind = 2)
            @test only(second.value) ≈ expected rtol=tolerance
        end
    end
end
@testset "Complex radial coordinates" begin
    for (precision, spheroid) in ((:double, :prolate), (:quad, :oblate))
        T = precision===:quad ? BigFloat : Float64
        tolerance = precision===:quad ? big"2e-26" : 3e-12
        sigma = spheroid===:prolate ? 1 : -1
        points = Complex{T}[1.6 + 0.3im, 2 - 0.2im]
        for c in (T(1.25), complex(T(1.25), T(0.125)))
            regular = rmn(1, 2, c, points; spheroid, precision, derivatives = 4)
            second = rmn(1, 2, c, points; spheroid, precision, kind = 2)
            for kind in 3:4
                sign = kind==3 ? 1 : -1
                actual = rmn(1, 2, c, points; spheroid, precision,
                    kind, scaled = true, logderivative = true)
                value = actual.value.mantissa .* T(10) .^ actual.value.exponent
                derivative = actual.derivative.mantissa .*
                             T(10) .^ actual.derivative.exponent
                @test value ≈ regular.value+sign*im*second.value rtol=tolerance
                @test derivative ≈ regular.derivative+sign*im*second.derivative rtol=tolerance
                @test actual.logderivative ≈ derivative ./ value rtol=tolerance
            end
            @test radial_wronskian(1, 2, c, points; spheroid, precision) ≈
                  inv.(c .* (points .^ 2 .- sigma)) rtol=tolerance
            @test radial_wronskian(
                1, 2, c, points; spheroid, precision, form = :normalized) ≈ ones(2) rtol=tolerance
            @test maximum(radial_wronskian(
                1, 2, c, points; spheroid, precision, form = :error)) < tolerance
            @test only(radial_wronskian(1, 2, c, first(points); spheroid, precision)) ≈
                  inv(c*(first(points)^2-sigma)) rtol=tolerance
            # The differential equation gives an independent second-derivative check.
            lambda = eigenvalue(1, 2, c; spheroid, precision)
            expected = @. -(2points*regular.derivative+(c^2*points^2-lambda-sigma/(points^2-sigma))*regular.value)/(points^2-sigma)
            @test regular.second_derivative ≈ expected rtol=tolerance
            @test rmn(1, 2, c, reverse(points); spheroid, precision).value ==
                  reverse(regular.value)
        end
        for z in Complex{T}[-2, -2 + 0.3im, 0.3 + 0.5im, 2 - 0.3im]
            root = spheroid===:prolate ? sqrt(z-1)*sqrt(z+1) : sqrt(1-im*z)*sqrt(1+im*z)
            regular = rmn(1, 1, 0, z; spheroid, precision, normalization = :static)
            @test only(regular.value) ≈ root/3 rtol=tolerance
            @test only(regular.derivative) ≈ z/(3root) rtol=tolerance
            result = rmn(0, 0, 0, z; spheroid, precision, kind = 2,
                normalization = :static, derivatives = 4)
            expected = spheroid===:prolate ? -atanh(inv(z)) : atan(z)-T(pi)/2
            @test only(result.value) ≈ expected rtol=tolerance
            @test only(result.derivative) ≈ inv(z^2-sigma) rtol=tolerance
            @test only(result.second_derivative) ≈ -2z/(z^2-sigma)^2 rtol=tolerance
            @test only(result.third_derivative) ≈ 2(3z^2+sigma)/(z^2-sigma)^3 rtol=tolerance
            @test only(result.fourth_derivative) ≈ -24z*(z^2+sigma)/(z^2-sigma)^4 rtol=tolerance
            @test only(radial_wronskian(
                0, 0, 0, z; spheroid, precision, normalization = :static, form = :error)) <
                  tolerance
        end
        for kind in 1:4
            real_result = rmn(
                1, 2, T(1.25), T[2]; spheroid, precision, kind, derivatives = 4)
            complex_result = rmn(
                1, 2, T(1.25), Complex{T}[2]; spheroid, precision, kind, derivatives = 4)
            @test all(isapprox(a, b; rtol = tolerance)
            for (a, b) in zip(real_result, complex_result))
        end
        ranged = rmn(0, 0:2, T(1.25), first(points); spheroid, precision,
            normalization = :static, second_derivative = true)
        @test size(ranged.value) == (1, 3)
        for n in 0:2
            @test ranged.value[:, n + 1] == rmn(
                0, n, T(1.25), first(points); spheroid, precision, normalization = :static).value
        end
        @test all(d -> d>=(precision===:quad ? 28 : 12), accuracy(
            0, 0, 1, points; spheroid, precision))
    end
    for kind in 1:4
        result = rmn(1, 2, 1.25+0.125im, 1+0im; kind, derivatives = 4, logderivative = true)
        reference = rmn(1, 2, 1.25+0.125im, 1; kind, derivatives = 4, logderivative = true)
        @test all(isequal(getproperty(result, key), getproperty(reference, key))
        for key in keys(reference))
    end
    @test only(rmn(0, 1, 1.25, 0im; spheroid = :oblate, logderivative = true).logderivative) ==
          Inf
    @test_throws DomainError rmn(0, 0, 1, 0im)
    @test_throws DomainError rmn(0, 0, 1, -1+0im)
    @test_throws DomainError rmn(0, 0, 1, 2im; spheroid = :oblate)
    @test_throws DomainError rmn(0, 0, 1, -2im; spheroid = :oblate)
    @test_throws ArgumentError rmn(0, 0, 1, 2im; normalization = :static, kind = 3)
    @test_throws ArgumentError radial_wronskian(0, 0, 1, 2im; form = :invalid)
    @test_throws ArgumentError accuracy(0, 0, 1, [2im]; target = :invalid)
end
