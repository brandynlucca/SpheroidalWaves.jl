using SpheroidalWaves,Test

@testset "Independent coordinate residual stencils" begin
    SW = SpheroidalWaves
    for precision in (:double, :quad)
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? T(1e-20) : T(1e-6)
        # P_3 and its derivatives give exact polynomial stencil references.
        x = T[0.2, 0.6]
        diagnostic = SW._spheroidal_residual(0, 3, zero(T), x; precision)
        @test diagnostic.derivative ≈ (15x.^2 .- 3)./2 rtol=tolerance
        @test diagnostic.second_derivative ≈ 15x rtol=tolerance
        @test maximum(diagnostic.relative_residual) < tolerance
        constant = SW._spheroidal_residual(0, 0, zero(T), T(0.3); precision)
        @test only(constant.residual) == only(constant.relative_residual) == 0
        for spheroid in (:prolate, :oblate)
            # Oblate x=0 requires a forward stencil, while x=2 is centered.
            points = spheroid === :oblate ? T[0, 2] : T[2]
            diagnostic = SW._spheroidal_residual(0, 0, T(1), points;
                                                target=:radial, spheroid, precision)
            direct = rmn(0, 0, T(1), points; spheroid, precision, second_derivative=true)
            @test diagnostic.derivative ≈ direct.derivative rtol=tolerance atol=tolerance
            @test diagnostic.second_derivative ≈ direct.second_derivative rtol=tolerance atol=tolerance
            @test maximum(diagnostic.relative_residual) < 10tolerance
            @test all(diagnostic.step .> 0)
        end
    end
    @test_throws ArgumentError SW._spheroidal_residual(0, 1, 1, 0.3; target=:bad)
    @test_throws DomainError SW._spheroidal_residual(0, 1, 1, 1)
    @test_throws DomainError SW._spheroidal_residual(0, 1, 1, 1; target=:radial)
    for h in (0, -1, Inf, NaN)
        @test_throws ArgumentError SW._spheroidal_residual(0, 1, 1, 0.3; h)
    end
    @test_throws r"below coordinate resolution" SW._spheroidal_residual(0, 1, 1, 0.3; h=eps(0.3)/8)
end

@testset "Taylor steps reject unconverged series" begin
    SW = SpheroidalWaves
    # Q_0(x)=atanh(x) has convergence radius one at x=0. A step to x=2
    # must fail its tail check rather than return a truncated polynomial.
    @test_throws r"did not converge" SW._qs_step(0, big"0.0", big"0.0", big"0.0",
                                                big"0.0", big"1.0", big"2.0")
    # At c=lambda=m=0 the radial equation admits log((x+1)/(x-1)).
    # Its Taylor series at x=2 cannot cross the singularity x=1.
    plan = (spheroid=:prolate, c=big"0", lambda=big"0", dlambda=big"0", m=0)
    @test_throws r"did not converge" SW._radial_sensitivity_step(plan, big"2",
        (log(big"3"), -big"2"/3, big"0", big"0"), -big"2")
end

@testset "Regular endpoint second derivatives" begin
    SW = SpheroidalWaves
    setprecision(BigFloat, 256) do
        # P_3^2(x)=15x(1-x^2), so its endpoint slopes are -30 and its
        # second derivatives are +90 at -1 and -90 at +1.
        spherical = SW._evaluate_coefficient_vector((m=2, degrees=[3]),
            [sqrt(big"240.0"/7)], [-1, 1]; second_derivative=true)
        @test spherical.value == [0, 0]
        @test spherical.derivative ≈ [-30, -30] rtol=big"1e-70"
        @test spherical.second_derivative ≈ [90, -90] rtol=big"1e-70"
        for precision in (:double, :quad), spheroid in (:prolate, :oblate)
            # At m=2, the endpoint limit of the angular equation gives
            # S''(x)=x*(lambda-sigma*c^2-3)*S'(x)/3 for x=+/-1.
            # c=5 selects coefficient evaluation through the public API.
            result = smn(2, 3, 5, [-1, 1]; precision, spheroid, second_derivative=true)
            lambda = eigenvalue(2, 3, 5; precision, spheroid)
            sigma = spheroid === :prolate ? 1 : -1
            expected = [-1, 1].*(lambda-sigma*25-3).*result.derivative./3
            @test result.second_derivative ≈ expected rtol=(precision === :quad ? big"1e-27" : 1e-12)
        end
        for precision in (:double, :quad), m in 0:5
            n, c = m+1, big"1.25"
            points = BigFloat[-1, -1+big"1e-8", 1-big"1e-8", 1]
            # A differentiated Legendre expansion provides an independent
            # reference for the native endpoint Taylor reconstruction.
            plan = SW._coefficient_plan(m, n, c; precision=:quad)
            scale = SW._coefficient_phase(plan, :quad)*sqrt(SW._ferrers_norm2(m,n,BigFloat))
            reference = SW._evaluate_coefficient_vector(plan, scale.*plan.v, points;
                                                        second_derivative=true)
            result = smn(m, n, c, points; precision, second_derivative=true, scaled=true)
            for field in (:value, :derivative, :second_derivative)
                encoded = getproperty(result, field)
                actual = encoded.mantissa .* BigFloat(10).^encoded.exponent
                expected = getproperty(reference, field)
                @testset "$precision m=$m $field" begin
                    for i in eachindex(points)
                        @test actual[i] ≈ expected[i] rtol=(precision === :quad ? big"1e-25" : 1e-11)
                    end
                end
            end
            # The regular radial solution has the same endpoint power.
            radial = rmn(m, n, c, BigFloat[1, 1+big"1e-8"]; precision, second_derivative=true)
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
    for (spheroid,value,ve,derivative,de) in (
        (:prolate,big"1.5256488534902150635105594465576701032672",-7,
                  big"-1.7954837073168552408315797203151455810388",-5),
        (:oblate,big"4.1896917085417494071462469474848952172303",-113,
                 big"1.5679082330112565415126277769813336590425",-110))
        ordinary=smn(0,0,375.,big".3";spheroid,precision=:double,scaled=true)
        @test only(ordinary.value.mantissa)≈value rtol=1e-13
        @test only(ordinary.value.exponent)==ve
        @test only(ordinary.derivative.mantissa)≈derivative rtol=1e-13
        @test only(ordinary.derivative.exponent)==de
    end
    for spheroid in (:prolate,:oblate),kind in (1,2)
        diagnostic=SpheroidalWaves._spheroidal_residual(0,0,big"1",big".3";
            spheroid,kind,precision=:quad,h=big"1e-5")
        actual=smn(0,0,big"1",big".3";spheroid,kind,precision=:quad,second_derivative=true)
        @test diagnostic.derivative≈actual.derivative rtol=big"1e-20"
        @test diagnostic.second_derivative≈actual.second_derivative rtol=big"1e-20"
        @test maximum(diagnostic.relative_residual)<big"1e-20"
    end
end
