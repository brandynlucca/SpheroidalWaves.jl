using SpheroidalWaves, Test

precision_api_cases = split(get(ENV, "SPHEROIDALWAVES_PRECISION_API_CASE", "accuracy,roots,curvature,complex-curvature,differences"), ',')
precision_api_settings = (precision(BigFloat), rounding(BigFloat))
if "accuracy" in precision_api_cases
    @testset "Arbitrary-precision accuracy estimates" begin
        for (target, c, points, kwargs) in (
                (:angular, 1, [1//4], (;)),
                (:angular, 1+1//8*im, [complex(1//4, 1//8)], (; kind = 2)),
                (:radial, 1, [2], (;)),
                (:radial, 1, [complex(2, 1//8)], (; spheroid = :oblate, kind = 2)),
                (:radial, 0, [2], (; normalization = :static)))
            result = accuracy(0, 0, c, points; target, kwargs..., precision_bits = 128, diagnostics = true)
            record = only(result.diagnostics)
            @test record.value.status===:estimated
            @test record.derivative.status===:estimated
            @test only(result.digits)>=36
            @test record.working_precision==(160, 192)
            @test precision(record.value.absolute_error)==128
            evaluate = target===:angular ? smn : rmn
            actual = evaluate(0, 0, c, points; kwargs..., precision_bits = 128)
            reference = evaluate(0, 0, c, points; kwargs..., precision_bits = 224)
            setprecision(BigFloat, 256) do
                for name in (:value, :derivative)
                    @test abs(only(getproperty(actual, name))-only(getproperty(reference, name)))<=getproperty(record, name).absolute_error
                end
            end
        end
        @test only(accuracy(0, 0, 1, [complex(1//4)]; target = :angular, precision_bits = 160))>33
        @test only(accuracy(0, 0, 1, [1//4]; target = :angular, precision_bits = 160))>33
        zero_value = accuracy(0, 1, 0, [0]; target = :angular, precision_bits = 80, diagnostics = true)
        @test only(zero_value.digits)==-1
        @test only(zero_value.diagnostics).value.absolute_error==0
        singular = accuracy(1, 1, 1, [1]; target = :angular, precision_bits = 80, diagnostics = true)
        @test only(singular.diagnostics).derivative.status===:singularity
        @test_throws ArgumentError accuracy(0, 0, 1, [0]; target = :bad, precision_bits = 80)
        @test_throws ArgumentError accuracy(0, 0, 1, [complex(0)]; target = :bad, precision_bits = 80)
        @test_throws ArgumentError accuracy(0, 0, 1, [0]; target = :angular, normalization = :static, precision_bits = 80)
        @test_throws ArgumentError accuracy(0, 0, 1, [0]; target = :angular, precision_bits = 1)
        failed = SpheroidalWaves._accuracy_at_precision(0, 0, 1, [0], :prolate, :double, :angular, 1,
            false, :standard, true, 80, 512; compare = _ -> nothing)
        @test only(failed.diagnostics).value.status===:refinement_failed
        SpheroidalWaves._with_swprecision(128) do
            estimate(a, b, r) = SpheroidalWaves._accuracy_component_at_precision(a, b, r, 80, false)
            @test estimate(Inf, 1, 1).status===:nonfinite
            @test estimate(1, Inf, 1).status===:refinement_failed
            @test estimate(1, 1, 2).status===:refinement_failed
        end
    end
end
if "roots" in precision_api_cases
    @testset "Arbitrary-precision coordinate zeros" begin
        for bits in (80, 160, 544)
            roots = setprecision(BigFloat, 32) do
                SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = bits)
            end
            @test length(roots)==2
            @test all(precision(x)==bits for x in roots)
            setprecision(BigFloat, bits+64) do
                reference = inv(sqrt(BigFloat(3)))
                @test abs(roots[1]+reference)<ldexp(BigFloat(1), -bits+2)
                @test abs(roots[2]-reference)<ldexp(BigFloat(1), -bits+2)
            end
        end
        stationary = SpheroidalWaves.angular_zeros(0, 3, 0; stationary = true, precision_bits = 96)
        setprecision(BigFloat, 160) do
            @test maximum(abs.(stationary .- [-1, 1] ./ sqrt(BigFloat(5))))<big"1e-28"
        end
        for spheroid in (:prolate, :oblate)
            angular = SpheroidalWaves.angular_zeros(0, 2, 1; spheroid, precision_bits = 80)
            @test maximum(abs, smn(0, 2, 1, angular; spheroid, precision_bits = 128).value)<1e-22
            for stationary in (false, true)
                roots = SpheroidalWaves.radial_zeros(0, 0, 1, (2, 5); spheroid, stationary, precision_bits = 80)
                @test !isempty(roots)
                values = rmn(0, 0, 1, roots; spheroid, precision_bits = 128)
                @test maximum(abs, getproperty(values, stationary ? :derivative : :value))<1e-22
            end
        end
        @test isempty(SpheroidalWaves.angular_zeros(0, 0, 1; precision_bits = 80))
        @test isempty(SpheroidalWaves.radial_zeros(0, 0, 1, (2, 21//10); precision_bits = 80))
        loose = SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = 80, rtol = 1//10^8)
        @test length(loose)==2
        @test_throws DomainError SpheroidalWaves.angular_zeros(0, 0, 0; stationary = true, precision_bits = 80)
        @test_throws ArgumentError SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = 1)
        @test_throws ArgumentError SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = 80, rtol = 0)
        @test_throws ArgumentError SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = 80, rtol = Inf)
        @test_throws ArgumentError SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = 80, max_points = 16)
        @test_throws ArgumentError SpheroidalWaves.radial_zeros(0, 0, 1, (2, 1); precision_bits = 80)
        @test_throws r"precision refinement" SpheroidalWaves._coordinate_zeros_at_precision(
            0, 2, 0, (-1, 1), :angular, :prolate, :double, 1, false,
            nothing, 16384, 80, 512; max_refinements = 1)
        tight = SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = 80, rtol = big"1e-100")
        @test tight==SpheroidalWaves.angular_zeros(0, 2, 0; precision_bits = 80)
        irregular = SpheroidalWaves.radial_zeros(0, 0, 1, (2, 5); kind = 2, precision_bits = 80)
        @test !isempty(irregular)
        @test maximum(abs, rmn(0, 0, 1, irregular; kind = 2, precision_bits = 128).value)<1e-22
    end
end
if "curvature" in precision_api_cases
    @testset "Arbitrary-precision wave curvature" begin
        for spheroid in (:prolate, :oblate), normalize in (false, true)
            value = jacobian_smn(0, 0, 0, [1//4]; spheroid, normalize, order = 2,
                precision_bits = 96, diagnostics = true)
            setprecision(BigFloat, 192) do
                factor = (spheroid===:prolate ? 1 : -1)/(normalize ? sqrt(BigFloat(2)) : BigFloat(1))
                @test abs(only(value.d2value_dc2)-factor*13/144)<big"1e-28"
                @test abs(only(value.d2derivative_dc2)+factor/6)<big"1e-28"
            end
            @test precision(only(value.d2value_dc2))==96
            @test value.metadata.method===:differenced_sensitivity
        end
        for evaluate in (jacobian_smn, jacobian_rmn)
            endpoint = evaluate(1, 1, 1, [1]; order = 2, precision_bits = 80, diagnostics = true)
            @test only(endpoint.d2value_dc2)==0
            @test isinf(real(only(endpoint.d2derivative_dc2)))
            @test endpoint.metadata.conditioning_flag===:singular
        end
        static = jacobian_rmn(1, 1, 0, [1]; normalization = :static, order = 2,
            precision_bits = 64, diagnostics = true)
        @test only(static.d2value_dc2)==0
        @test isinf(real(only(static.d2derivative_dc2)))
        @test_throws r"step refinement" SpheroidalWaves._precision_parameter_stencil(
            (z, _) -> reshape([exp(100z)], 1, 1), 0, 80, nothing, 2, true, 1e-6, 1e-10;
            max_refinements = 1)
        kinds = [jacobian_rmn(0, 0, 1, [2]; kind, order = 2, precision_bits = 64) for kind in 1:4]
        @test kinds[3].d2value_dc2≈kinds[1].d2value_dc2+im*kinds[2].d2value_dc2
        @test kinds[4].d2derivative_dc2≈kinds[1].d2derivative_dc2-im*kinds[2].d2derivative_dc2
        for mode in (RoundDown, RoundUp)
            result = setprecision(BigFloat, 32) do
                setrounding(BigFloat, mode) do
                    jacobian_smn(0, 0, 0, [0]; order = 2, precision_bits = 80)
                end
            end
            expected = setprecision(BigFloat, 256) do
                BigFloat(BigFloat(1)/9, mode; precision = 80)
            end
            @test only(result.d2value_dc2)==expected
        end
        tasks = [Threads.@spawn begin
            value = jacobian_smn(0, 0, 0, [0]; order = 2, precision_bits = p)
            (precision(only(value.d2value_dc2)), precision(BigFloat), rounding(BigFloat))
        end for p in (64, 96)]
        @test fetch.(tasks)==[(p, precision(BigFloat), rounding(BigFloat)) for p in (64, 96)]
    end
end
if "complex-curvature" in precision_api_cases
    @testset "Arbitrary-precision complex wave curvature" begin
        for (evaluate, point, kwargs) in ((jacobian_smn, 1//4, (;)),
                (jacobian_smn, complex(1//4, 1//8), (; kind = 2)),
                (jacobian_rmn, 2, (;)),
                (jacobian_rmn, complex(2, 1//8), (; spheroid = :oblate, normalization = :static)))
            c = complex(1//1, 1//8)
            result = evaluate(0, 0, c, [point]; kwargs..., order = 2, precision_bits = 80, diagnostics = true)
            @test result.d2value_dcreal_dcimag≈im .* result.d2value_dcreal2
            @test result.d2value_dcimag2≈-result.d2value_dcreal2
            @test result.d2derivative_dcreal_dcimag≈im .* result.d2derivative_dcreal2
            @test result.d2derivative_dcimag2≈-result.d2derivative_dcreal2
            setprecision(BigFloat, 256) do
                h = BigFloat(2)^(-60)
                high = evaluate(0, 0, c+h, [point]; kwargs..., precision_bits = 192)
                low = evaluate(0, 0, c-h, [point]; kwargs..., precision_bits = 192)
                @test abs(only(result.d2value_dcreal2)-only(high.dvalue_dcreal-low.dvalue_dcreal)/(2h))<big"1e-23"
                @test abs(only(result.d2derivative_dcreal2)-only(high.dderivative_dcreal-low.dderivative_dcreal)/(2h))<big"1e-23"
            end
            @test precision(real(only(result.d2value_dcreal2)))==80
        end
    end
end
if "differences" in precision_api_cases
    @testset "Arbitrary-precision finite differences" begin
        for c in (1//1, complex(1//1, 1//8)), order in (1, 2), adaptive in (false, true)
            result = jacobian_eigen(0, 0, c; precision_bits = 96, h = 1//10^10,
                order, adaptive, diagnostics = true)
            analytic = jacobian_eigen(0, 0, c; precision_bits = 160, order)
            name = order==1 ? :d_dcreal : :d2_dcreal2
            actual = c isa Real ? result.derivative : getproperty(result, name)
            reference = c isa Real ? analytic : getproperty(analytic, name)
            @test isapprox(actual, reference; rtol = 1e-18)
            @test precision(real(actual))==96
            @test c isa Real || order==2 || isapprox(result.d_dcimag, analytic.d_dcimag; rtol = 1e-18)
        end
        for (evaluate, point, kwargs) in ((jacobian_smn, 1//4, (;)),
                (jacobian_smn, complex(1//4, 1//8), (; kind = 2)),
                (jacobian_rmn, 2, (;)),
                (jacobian_rmn, complex(2, 1//8), (; normalization = :static)))
            for c in (1, complex(1, 1//8))
                result = evaluate(0, 0, c, [point]; kwargs..., precision_bits = 80,
                    h = 1//10^8, diagnostics = true)
                analytic = evaluate(0, 0, c, [point]; kwargs..., precision_bits = 128)
                for name in keys(analytic)
                    @test isapprox(getproperty(result, name), getproperty(analytic, name); rtol = 1e-14)
                    @test precision(real(only(getproperty(result, name))))==80
                end
            end
        end
        for operator in (:concentration, :fourier), order in (1, 2), c in (0, 1)
            result = jacobian_eigen(0, 0, c; operator, order, precision_bits = 80, h = 1//10^8)
            analytic = jacobian_eigen(0, 0, c; operator, order, precision_bits = 128)
            @test isapprox(result, analytic; rtol = 1e-14, atol = 1e-14)
        end
        for form in (:complement, :log)
            result = jacobian_eigen(0, 0, 1; operator = :concentration, form, precision_bits = 80, h = 1//10^8)
            analytic = jacobian_eigen(0, 0, 1; operator = :concentration, form, precision_bits = 128)
            @test isapprox(result, analytic; rtol = 1e-14)
        end
        for h in (0, -1, Inf, NaN)
            @test_throws ArgumentError jacobian_eigen(0, 0, 1; precision_bits = 80, h)
        end
        @test_throws ArgumentError jacobian_eigen(0, 0, 1; operator = :concentration, precision_bits = 80, h = 1)
        @test_throws DomainError jacobian_eigen(0, 0, 0; operator = :concentration, form = :log, precision_bits = 80, h = 1//100)
        @test_throws ArgumentError jacobian_rmn(0, 0, 1, [2]; precision_bits = 80, h = 1, order = 2)
        @test_throws ArgumentError jacobian_smn(0, 0, 1, [complex(0)]; precision_bits = 80, order = 3)
        @test_throws ArgumentError jacobian_rmn(0, 0, 1, [complex(2)]; precision_bits = 80, order = 3)
        poor = jacobian_eigen(0, 0, 1; precision_bits = 80, h = 1//4,
            rtol = 1e-30, diagnostics = true)
        @test poor.metadata.conditioning_flag===:poor
        @test poor.metadata.suggested_action===:retry_smaller_h
        @test poor.metadata.step_used==1//16
        for (evaluate, kwargs) in ((jacobian_smn, (;)), (jacobian_rmn, (; normalization = :static)))
            endpoint = evaluate(1, 1, 0, [1]; kwargs..., precision_bits = 80,
                h = 1//1000, diagnostics = true)
            @test only(endpoint.dvalue_dc)==0
            @test only(endpoint.dderivative_dc)==0
            result = evaluate(0, 0, 1, [1]; kwargs..., precision_bits = 80,
                h = 1//1000, adaptive = false)
            @test precision(real(only(result.dvalue_dc)))==80
        end
        endpoint = jacobian_smn(1, 1, complex(1, 1//8), [complex(1)]; precision_bits = 80,
            h = 1//1000, diagnostics = true)
        @test !endpoint.metadata_derivative_dcimag.finite_flag
        @test endpoint.metadata_derivative_dcimag.conditioning_flag===:singular
        @test keys(endpoint)[1:4]==(:dvalue_dcreal, :dvalue_dcimag, :dderivative_dcreal, :dderivative_dcimag)
        low_context = setprecision(BigFloat, 32) do
            jacobian_rmn(1, 1, 1, [1]; precision_bits = 96, h = 1//10^12)
        end
        @test isinf(real(only(low_context.dderivative_dc)))
        SpheroidalWaves._with_swprecision(128) do
            metadata = SpheroidalWaves._precision_stencil_metadata([SpheroidalWaves._SWFloat(Inf)], [SpheroidalWaves._SWFloat(Inf)],
                SpheroidalWaves._SWFloat(1//10), 80, 1e-6, 1e-10, 1)
            @test !metadata.finite_flag && metadata.conditioning_flag===:poor
        end
    end
end
@test (precision(BigFloat), rounding(BigFloat))==precision_api_settings
