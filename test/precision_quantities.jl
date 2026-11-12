using SpheroidalWaves, Test

@testset "Arbitrary-precision quantities" begin
    settings = (precision(BigFloat), rounding(BigFloat))
    setprecision(BigFloat, 320) do
        bits = 160
        tolerance = BigFloat(2)^(-bits+6)
        @testset "Integral reference and forms" begin
            # Independent sinc-kernel quadrature reference from integral_eigenvalues.jl, confirmed beyond 38 digits.
            reference = parse(BigFloat, "0.5725817806378951222239686534549344836027845059355985016921")
            concentration, logarithm, complement, fourier = setprecision(BigFloat, 80) do
                (eigenvalue(0, 0, 1; operator = :concentration, precision_bits = bits),
                    eigenvalue(0, 0, 1; operator = :concentration, form = :log, precision_bits = bits),
                    eigenvalue(0, 0, 1; operator = :concentration, form = :complement, precision_bits = bits),
                    eigenvalue(0, 0, 1; operator = :fourier, precision_bits = bits))
            end
            @test isapprox(concentration, reference; rtol = big"1e-37")
            @test isapprox(exp(logarithm), concentration; rtol = tolerance)
            @test isapprox(concentration+complement, 1; rtol = tolerance)
            @test isapprox(abs2(fourier)/(2BigFloat(pi)), concentration; rtol = tolerance)
            @test all(precision(real(x))==bits for x in (concentration, logarithm, complement, fourier))
            @test precision(imag(fourier))==bits
            wider = eigenvalue(0, 0, 1; operator = :concentration, precision_bits = 224)
            @test isapprox(wider, concentration; rtol = tolerance)
        end

        @testset "Exact limits and analytic derivatives" begin
            for spheroid in (:prolate, :oblate)
                @test amn(-1, 2, 0; spheroid, precision_bits = bits)==1
                @test amn(1, 2, 0; spheroid, precision_bits = bits)==1
                @test kmn(0, 0, 0; spheroid, precision_bits = bits)==1
                @test kmn(1, 2, 0; spheroid, precision_bits = bits)==0
                @test jacobian_eigen(0, 0, 0; spheroid, precision_bits = bits)==0
                curvature = jacobian_eigen(0, 0, 0; spheroid, precision_bits = bits, order = 2)
                @test precision(curvature)==bits
                @test isapprox(curvature, (spheroid===:prolate ? 1 : -1)*BigFloat(2)/3; rtol = tolerance)
            end
            for order in (1, 2)
                derivative = jacobian_eigen(0, 0, 0; operator = :concentration,
                    precision_bits = bits, order, diagnostics = true)
                @test derivative.metadata.method===:right_limit
                @test isapprox(derivative.derivative, order==1 ? 2/BigFloat(pi) : 0; rtol = tolerance)
            end
            @test jacobian_eigen(0, 0, 0; operator = :concentration,
                form = :log, precision_bits = bits)==Inf
            @test isfinite(jacobian_eigen(0, 0, 1; precision_bits = bits, h = 1//100))
        end

        @testset "Directed rounding of complex Jacobians" begin
            parameter = complex(1//1, 1//4)
            evaluators = (
                p -> jacobian_eigen(0, 0, parameter; precision_bits = p),
                p -> jacobian_eigen(0, 0, parameter; precision_bits = p, order = 2),
                p -> jacobian_rmn(0, 0, parameter, [2]; precision_bits = p))
            for evaluate in evaluators
                reference = setrounding(BigFloat, RoundNearest) do
                    evaluate(256)
                end
                for mode in (RoundDown, RoundUp)
                    rounded(z) = complex(BigFloat(real(z), mode; precision = bits),
                        BigFloat(imag(z), mode; precision = bits))
                    actual = setprecision(BigFloat, 80) do
                        setrounding(BigFloat, mode) do
                            evaluate(bits)
                        end
                    end
                    for name in keys(reference)
                        high, value = getproperty(reference, name), getproperty(actual, name)
                        if high isa AbstractArray
                            @test value==rounded.(high)
                            @test all(precision(real(z))==bits && precision(imag(z))==bits for z in value)
                        else
                            @test value==rounded(high)
                            @test precision(real(value))==precision(imag(value))==bits
                        end
                    end
                end
            end
        end

        @testset "Inverse and continuation" begin
            target = eigenvalue(0, 0, 1; precision_bits = 224)
            solved = find_c_for_eigenvalue(0, 0, target; bracket = (1//2, 3//2), precision_bits = bits)
            @test solved.converged
            @test precision(solved.c)==bits
            @test isapprox(solved.c, 1; rtol = tolerance)
            concentration = eigenvalue(0, 0, 1; operator = :concentration, precision_bits = 224)
            bandwidth = find_c_for_eigenvalue(0, 0, concentration; bracket = (1//2, 3//2),
                operator = :concentration, precision_bits = bits)
            @test bandwidth.converged
            @test isapprox(bandwidth.c, 1; rtol = tolerance)
            unfinished = find_c_for_eigenvalue(0, 0, target; bracket = (1//2, 2),
                precision_bits = bits, maxiter = 1)
            @test !unfinished.converged
            @test unfinished.method===:maxiter
            real_path = eigenvalue_sweep(0, 0, [1//1, 5//4]; precision_bits = bits, branch_lock = false)
            @test all(==(0), real_path.selected_n)
            @test all(precision(x)==bits for x in [real_path.c; real_path.lambda])
            path = complex.([1//1, 1//1, 1//1], [0//1, 1//8, 0//1])
            complex_path = eigenvalue_sweep(0, 0, path; precision_bits = bits)
            @test isapprox(first(complex_path.lambda), last(complex_path.lambda); rtol = tolerance)
            @test isapprox(real(first(complex_path.lambda)), target; rtol = tolerance)
            @test all(precision(real(x))==bits && precision(imag(x))==bits
                for x in [complex_path.c; complex_path.lambda])
        end

        @testset "Wave dispatch and combined results" begin
            ranged = smn(0, 0:1, 0, [1//2]; precision_bits = bits)
            @test size(ranged.value)==(1, 2)
            @test ranged.value==[1 1//2]
            @test all(precision(x)==bits for x in ranged.value)
            combined = spheroidal(0, 0:1, 0; eta = 1//2, coefficients = true,
                derivatives = false, precision_bits = bits)
            @test combined.angular.value==ranged.value
            @test !hasproperty(combined.angular, :derivative)
            @test length(combined.coefficients)==2
            @test all(x -> x.dcoefficients_dc===nothing, combined.coefficients)
            z = complex(1//2, 1//8)
            angular = smn(0, 1, 0, [z]; precision_bits = bits, derivatives = 4)
            @test only(angular.value)==z
            @test only(angular.derivative)==1
            @test only(angular.fourth_derivative)==0
            @test all(isfinite, jacobian_smn(0, 0, 1, [0]; precision_bits = bits, order = 2).d2value_dc2)
            @test all(isfinite, jacobian_rmn(0, 0, 1, [2]; precision_bits = bits, h = 1//100).dvalue_dc)
            wronskian = radial_wronskian(0, 0, 0, [2]; normalization = :static,
                form = :normalized, precision_bits = bits)
            @test isapprox(only(wronskian), 1; rtol = tolerance)
            @test precision(real(only(wronskian)))==bits
        end
    end
    @test (precision(BigFloat), rounding(BigFloat))==settings
end
