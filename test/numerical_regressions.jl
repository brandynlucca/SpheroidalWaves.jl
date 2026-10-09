using SpheroidalWaves, Test

@testset "Coordinate zeros and stationary points" begin
    SW = SpheroidalWaves
    for precision in (:double, :quad)
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? T(1e-27) : T(1e-11)
        # P_3'(x) = (15x^2-3)/2 and P_2(x) = (3x^2-1)/2.
        @test SW.angular_zeros(0, 3, 0; stationary = true, precision) ≈
              [-inv(sqrt(T(5))), inv(sqrt(T(5)))] atol=tolerance
        @test isempty(SW.angular_zeros(2, 2, 0; precision))
        for spheroid in (:prolate, :oblate), stationary in (false, true), kind in (1, 2)
            roots = SW.radial_zeros(0, 0, T(2), (T(2), T(4));
                precision, spheroid, kind, stationary)
            @test !isempty(roots)
            @test issorted(roots) && all(2 .<= roots .<= 4)
            values = rmn(0, 0, T(2), roots; precision, spheroid, kind)
            @test all(abs.(getproperty(values, stationary ? :derivative : :value)) .<
                      20tolerance)
        end
    end
    @test_throws DomainError SW.angular_zeros(0, 0, 0; stationary = true)
    for rtol in (0, -1, NaN, Inf)
        @test_throws ArgumentError SW.angular_zeros(0, 1, 0; rtol)
        @test_throws ArgumentError SW.radial_zeros(0, 0, 1, (2, 3); rtol)
    end
    @test_throws ArgumentError SW.angular_zeros(0, 1, 0; max_points = 32)
    @test_throws ArgumentError SW.radial_zeros(0, 0, 1, (2, 3); max_points = 32)
    @test_throws ArgumentError SW.radial_zeros(0, 0, 1, (2, 3); kind = 3)
    @test_throws ArgumentError SW.radial_zeros(0, 0, 1, (1, 2))
    @test_throws ArgumentError SW.radial_zeros(0, 0, 1, (3, 2))
    # Synthetic polynomials exercise exact grid roots and failed searches.
    polynomial(x) = x .* (x .- 1)
    @test SW._coordinate_zeros(polynomial, 0.0, 1.0, 1e-12, 4, 8) == [0.0, 1.0]
    @test SW._bisect_wave(identity, 0.0, 1.0, 0.0, 1.0, 1e-12) == 0
    @test SW._bisect_wave(identity, -1.0, 0.0, -1.0, 0.0, 1e-12) == 0
    @test_throws ErrorException SW._bisect_wave(
        x -> fill(NaN, length(x)), -1.0, 1.0, -1.0, 1.0, 1e-12)
    @test_throws ErrorException SW._coordinate_zeros(
        x -> fill(NaN, length(x)), 0.0, 1.0, 1e-12, 4, 8)
    @test_throws DomainError SW._coordinate_zeros(x -> zero(x), 0.0, 1.0, 1e-12, 4, 8)
    @test_throws ErrorException SW._coordinate_zeros(polynomial, 0.0, 1.0, 1e-12, 4, 4)
    # An unattainable tolerance must exhaust the bounded bisection search.
    setprecision(BigFloat, 2048) do
        @test_throws r"did not converge" SW._bisect_wave(x -> 3x .- 1,
            big"0", big"1", big"-1", big"2", big"1e-1000")
    end
end

@testset "Focused numerical regressions" begin
    for precision in (:double, :quad)
        lib = SpheroidalWaves.backend_library(; precision)
        (lib === nothing || !isfile(lib)) && continue
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-27" : big"2e-11"

        # Independent spectral value, embedded here so CI need not reconstruct
        # a high-precision reference or execute the extensive path matrices.
        c = 2.0+3.0im
        expected_lambda = complex(big"1.8526077784286968122645450083859036836748134683954",
            big"3.1017112824466201071716399658076138048679053435607")
        expected_s = complex(big"1.4672095477135660718169187632852425554599094354553",
            big"0.02551066455758449011154602903404097070828353999752")
        @test eigenvalue(0, 0, c; precision) ≈ expected_lambda rtol=tolerance
        @test only(smn(0, 0, c, T(3)/10; precision).value) ≈ expected_s rtol=tolerance
        @test only(radial_wronskian(0, 0, c, T(2); precision, form = :normalized)) ≈ 1 rtol=100tolerance

        for spheroid in (:prolate, :oblate)
            parameter = T(5)/4
            r = rmn(1, 2, parameter, T(2); precision, spheroid,
                scaled = true, logderivative = true)
            ordinary = rmn(1, 2, parameter, T(2); precision, spheroid)
            @test r.value.mantissa .* BigFloat(10) .^ r.value.exponent ≈ ordinary.value rtol=tolerance
            @test r.logderivative ≈ ordinary.derivative ./ ordinary.value rtol=tolerance
            @test only(radial_wronskian(
                1, 2, parameter, T(2); precision, spheroid, form = :error)) < tolerance
        end

        # Exact spherical polynomials independently check coordinate derivatives
        # and zeros, without solving a numerical reference eigenproblem.
        s = smn(0, 3, zero(T), T(1)/3; precision, second_derivative = true)
        @test only(s.second_derivative) ≈ T(5) rtol=tolerance
        roots = SpheroidalWaves.angular_zeros(0, 3, zero(T); precision)
        a = sqrt(T(3)/5)
        @test roots≈T[-a, 0, a] rtol=tolerance atol=tolerance
    end

    # Original exact-decimal radial boundary regressions, kept in the normal
    # suite even though the extensive boundary integration tests are local.
    lib = SpheroidalWaves.backend_library(; precision = :quad)
    if lib !== nothing && isfile(lib)
        for (m, n, value) in ((
            0, 1, big"-28.19279152071340627099348388765367655558239286251609457"),
            (1, 2, big"-2846358.502301398239201169080317823638197684953289129961"))
            @test only(rmn(
                m, n, big"1.25", big"1.000000000001"; precision = :quad, kind = 2).value) ≈
                  value rtol=big"1e-27"
        end
    end

    s = smn(200, 200, 0.0, 0.3; scaled = true, logderivative = true)
    @test only(s.value.exponent) > 308
    @test only(s.logderivative) ≈ -200*0.3/(1-0.3^2) rtol=2e-13
end
