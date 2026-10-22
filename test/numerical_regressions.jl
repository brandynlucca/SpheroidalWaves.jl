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

    # Rodrigues' formula gives P_200^200(x) = 399!!*(1-x^2)^100 exactly.
    # Preserve the represented Float64 input as a rational before evaluating it.
    x = Rational{BigInt}(0.3)
    exact = prod(big(k) for k in 1:2:399)*(1-x^2)^100
    slope = -200x/(1-x^2)*exact
    expected = BigFloat(numerator(exact))/BigFloat(denominator(exact))
    dexpected = BigFloat(numerator(slope))/BigFloat(denominator(slope))
    for precision in (:double, :quad)
        s = smn(200, 200, 0.0, 0.3; precision, scaled = true, logderivative = true)
        tolerance = precision === :quad ? big"1e-28" : 2e-13
        @test only(s.value.exponent) > 308
        @test only(s.logderivative) ≈ -200*0.3/(1-0.3^2) rtol=2e-13
        @test only(s.value.mantissa)*big(10)^only(s.value.exponent) ≈ expected rtol=tolerance
        @test only(s.derivative.mantissa)*big(10)^only(s.derivative.exponent) ≈ dexpected rtol=tolerance
    end
end

@testset "Nonzero-parameter angular zeros and stationary points" begin
    # Independent mpmath 1.3.0 dense Legendre eigensolves, m=0 and c=1.
    # scripts/generate_accuracy_references.py reproduces these fixtures without
    # loading this package. Precision 70 -> 100 digits and truncation 24 -> 32
    # terms changed all entries by less than 8e-71 and 6e-75, respectively.
    # DLMF 30.4 normalization: integral S^2 dx = 2/(2n+1), with P_n phase.
    # Each tuple stores geometry, n, k, lambda, S(x1), S'(x1), S(x2), S'(x2).
    # x1=k/2^32 and x2=(k+1)/2^32 are exact in both tested precisions.
    references = (
        (:prolate, 2, 2402525102,
            big"6.5334718005237964814925379360799328556369815722767",
            big"-0.0000000003646498962489905788454236918058790314910610385883",
            big"1.703323695364498173128344421120509893166775201302",
            big"0.000000000031936056168787758739232927394918200909331216420251",
            big"1.7033236960102420046654208441424825333890648774496"),
        (:prolate, 3, 1889698684,
            big"12.514462145094064809261980883628730057332910207698",
            big"-0.44292097656532255188482591614237266717618479737531",
            big"-0.00000000038207851043720495748415690501344184277370260476039",
            big"-0.44292097656532255179036066283900872297066311925656",
            big"0.0000000011935288576588813378344880644357212632729745748638"),
        (:oblate, 2, 2559710628,
            big"5.4868000538186858525400642095821339086324379566395",
            big"-0.00000000022586337863787453748872263184262640129052089928229",
            big"1.7661129448340532706054786121628771889090325546217",
            big"0.00000000018534183516176014934265906632273885987885759586284",
            big"1.7661129455941844827844704379333764287770031145608"),
        (:oblate, 3, 1952290389,
            big"11.492120902297953778872839994203848799719935527161",
            big"-0.45130974952979024172775506093481452827491145456698",
            big"-0.00000000088576585426959176119068231968610374136907004850296",
            big"-0.45130974952979024175361115885493258229135506392983",
            big"0.00000000066366366447389452958718345957023128915676149551127")
    )
    for precision in (:double, :quad),
        (spheroid, n, k, lambda, v1, d1, v2, d2) in references

        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-28" : 2e-14
        points = T[k, k + 1] ./ T(2)^32
        @test eigenvalue(0, n, 1; precision, spheroid) ≈ lambda rtol=tolerance
        for normalize in (false, true)
            result = smn(0, n, 1, points; precision, spheroid, normalize)
            scale = normalize ? sqrt(BigFloat(2n+1)/2) : BigFloat(1)
            for (i, (value, derivative)) in enumerate(((v1, d1), (v2, d2)))
                # Absolute error for the vanishing component, relative for the other.
                if n == 2
                    @test abs(result.value[i]-scale*value) <= tolerance*scale
                    @test result.derivative[i]≈scale*derivative rtol=tolerance atol=0
                else
                    @test result.value[i]≈scale*value rtol=tolerance atol=0
                    @test abs(result.derivative[i]-scale*derivative) <= tolerance*scale
                end
            end
        end
    end
end

@testset "Absolute errors near angular zeros" begin
    # Exact P_2 and P_3 polynomials, independent of the solver's recurrence.
    # Rational arithmetic captures the actual input bits. Only the final ratio
    # and unit-norm factor use BigFloat arithmetic (256 bits in the test suite).
    # Condon-Shortley phase is +1 for m=0. DLMF 14.7.1 and 14.7.13.
    for precision in (:double, :quad), spheroid in (:prolate, :oblate),
        normalize in (false, true), n in (2, 3)
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-28" : 2e-14
        delta = T(2)^(precision === :quad ? -80 : -20)
        root = inv(sqrt(T(n == 2 ? 3 : 5)))
        points = T[-root - delta, -root + delta, 0, root - delta, root + delta]
        result = smn(0, n, 0, points; precision, spheroid, normalize)
        scale = normalize ? sqrt(BigFloat(2n+1)/2) : BigFloat(1)
        for (i, point) in enumerate(points)
            x = Rational{BigInt}(point)
            value = n == 2 ? (3x^2-1)/2 : (5x^3-3x)/2
            derivative = n == 2 ? 3x : (15x^2-3)/2
            expected = scale*BigFloat(numerator(value))/BigFloat(denominator(value))
            dexpected = scale*BigFloat(numerator(derivative))/BigFloat(denominator(derivative))
            # Absolute error at value zeros or stationary points, relative error
            # for the other component. Scaling accounts for unit normalization.
            if n == 2 || iszero(point)
                @test abs(result.value[i]-expected) <= tolerance*scale
            else
                @test result.value[i]≈expected rtol=tolerance atol=0
            end
            if n == 3 || iszero(point)
                @test abs(result.derivative[i]-dexpected) <= tolerance*scale
            else
                @test result.derivative[i]≈dexpected rtol=tolerance atol=0
            end
        end
    end
end
