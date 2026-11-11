using SpheroidalWaves, Test

let cases = split(get(ENV, "SPHEROIDALWAVES_PRECISION_API_CASE", "accuracy,roots,curvature,differences"), ',')
    SW = SpheroidalWaves
    settings = (precision(BigFloat), rounding(BigFloat))
    if "accuracy" in cases
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
            failed = SW._accuracy_at_precision(0, 0, 1, [0], :prolate, :double, :angular, 1,
                false, :standard, true, 80, 512; compare = _ -> nothing)
            @test only(failed.diagnostics).value.status===:refinement_failed
        end
    end
    if "roots" in cases
        @testset "Arbitrary-precision coordinate zeros" begin
            for bits in (80, 160, 544)
                roots = setprecision(BigFloat, 32) do
                    SW.angular_zeros(0, 2, 0; precision_bits = bits)
                end
                @test length(roots)==2
                @test all(precision(x)==bits for x in roots)
                setprecision(BigFloat, bits+64) do
                    reference = inv(sqrt(BigFloat(3)))
                    @test abs(roots[1]+reference)<ldexp(BigFloat(1), -bits+2)
                    @test abs(roots[2]-reference)<ldexp(BigFloat(1), -bits+2)
                end
            end
            stationary = SW.angular_zeros(0, 3, 0; stationary = true, precision_bits = 96)
            setprecision(BigFloat, 160) do
                @test maximum(abs.(stationary .- [-1, 1] ./ sqrt(BigFloat(5))))<big"1e-28"
            end
            for spheroid in (:prolate, :oblate)
                angular = SW.angular_zeros(0, 2, 1; spheroid, precision_bits = 80)
                @test maximum(abs, smn(0, 2, 1, angular; spheroid, precision_bits = 128).value)<1e-22
                for stationary in (false, true)
                    roots = SW.radial_zeros(0, 0, 1, (2, 5); spheroid, stationary, precision_bits = 80)
                    @test !isempty(roots)
                    values = rmn(0, 0, 1, roots; spheroid, precision_bits = 128)
                    @test maximum(abs, getproperty(values, stationary ? :derivative : :value))<1e-22
                end
            end
            @test isempty(SW.angular_zeros(0, 0, 1; precision_bits = 80))
            @test isempty(SW.radial_zeros(0, 0, 1, (2, 21//10); precision_bits = 80))
            loose = SW.angular_zeros(0, 2, 0; precision_bits = 80, rtol = 1//10^8)
            @test length(loose)==2
            @test_throws DomainError SW.angular_zeros(0, 0, 0; stationary = true, precision_bits = 80)
            @test_throws ArgumentError SW.angular_zeros(0, 2, 0; precision_bits = 1)
            @test_throws ArgumentError SW.angular_zeros(0, 2, 0; precision_bits = 80, rtol = 0)
            @test_throws ArgumentError SW.angular_zeros(0, 2, 0; precision_bits = 80, rtol = Inf)
            @test_throws ArgumentError SW.angular_zeros(0, 2, 0; precision_bits = 80, max_points = 16)
            @test_throws ArgumentError SW.radial_zeros(0, 0, 1, (2, 1); precision_bits = 80)
        end
    end
    if "curvature" in cases
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
            for evaluate in (jacobian_smn, jacobian_rmn)
                endpoint = evaluate(1, 1, 1, [1]; order = 2, precision_bits = 80, diagnostics = true)
                @test only(endpoint.d2value_dc2)==0
                @test isinf(real(only(endpoint.d2derivative_dc2)))
                @test endpoint.metadata.conditioning_flag===:singular
            end
        end
    end
    if "differences" in cases
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
            for h in (0, -1, Inf, NaN)
                @test_throws ArgumentError jacobian_eigen(0, 0, 1; precision_bits = 80, h)
            end
            @test_throws ArgumentError jacobian_eigen(0, 0, 1; operator = :concentration, precision_bits = 80, h = 1)
            @test_throws DomainError jacobian_eigen(0, 0, 0; operator = :concentration, form = :log, precision_bits = 80, h = 1//100)
            @test_throws ArgumentError jacobian_rmn(0, 0, 1, [2]; precision_bits = 80, h = 1, order = 2)
            @test_throws ArgumentError jacobian_smn(0, 0, 1, [complex(0)]; precision_bits = 80, order = 3)
            @test_throws ArgumentError jacobian_rmn(0, 0, 1, [complex(2)]; precision_bits = 80, order = 3)
        end
    end
    @test (precision(BigFloat), rounding(BigFloat))==settings
end
