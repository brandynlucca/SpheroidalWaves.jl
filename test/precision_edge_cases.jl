using SpheroidalWaves, Test

@testset "Arbitrary-precision endpoints and continuation" begin
    SW = SpheroidalWaves
    bits = 96
    setprecision(BigFloat, 256) do
        tolerance = BigFloat(2)^(-bits+8)
        @testset "Integral limits and curvature" begin
            for n in 0:2
                @test eigenvalue(0, n, 0; operator = :fourier, precision_bits = bits)==(n==0 ? 2 : 0)
                for (form, expected) in ((:value, 0), (:complement, 1), (:log, -Inf))
                    @test eigenvalue(0, n, 0; operator = :concentration, form, precision_bits = bits)==expected
                end
                @test isapprox(jacobian_eigen(0, n, 0; operator = :fourier, precision_bits = bits),
                    n==1 ? -2im/BigFloat(3) : 0; rtol = tolerance)
            end
            c, h = BigFloat(1)/4, BigFloat(2)^(-70)
            for operator in (:fourier, :concentration), order in (1, 2)
                value = jacobian_eigen(0, 0, c; operator, order, precision_bits = bits)
                reference = order==1 ?
                    (eigenvalue(0, 0, c+h; operator, precision_bits = 224)-
                     eigenvalue(0, 0, c-h; operator, precision_bits = 224))/(2h) :
                    (jacobian_eigen(0, 0, c+h; operator, precision_bits = 224)-
                     jacobian_eigen(0, 0, c-h; operator, precision_bits = 224))/(2h)
                @test isapprox(value, reference; rtol = tolerance)
            end
        end

        @testset "Inverse iterations and zero endpoints" begin
            for (form, target) in ((:value, 0), (:complement, 1), (:log, -Inf))
                result = find_c_for_eigenvalue(0, 0, target; operator = :concentration,
                    form, bracket = (0, 2), precision_bits = bits)
                @test result.converged && result.c==0 && result.iterations==0
                @test_throws ArgumentError find_c_for_eigenvalue(0, 0, target;
                    operator = :concentration, form, bracket = (1, 2), precision_bits = bits)
            end
            for operator in (:separation, :concentration), use_jacobian in (false, true)
                target = eigenvalue(0, 0, 1; operator, precision_bits = 160)
                result = find_c_for_eigenvalue(0, 0, target; operator, use_jacobian,
                    bracket = (0, 3), precision_bits = bits)
                @test result.converged && result.iterations>1
                @test isapprox(result.c, 1; rtol = tolerance)
            end
            path = [1//1, 9//8, 5//4]
            result = eigenvalue_sweep(0, 0, path; precision_bits = bits, branch_lock = false)
            @test result.lambda==[eigenvalue(0, 0, c; precision_bits = bits) for c in path]
            complex_path = complex.(path, 1//8)
            unlocked = eigenvalue_sweep(0, 0, complex_path; precision_bits = bits, branch_lock = false)
            @test unlocked.selected_n==[0, 0, 0]
            @test all(isfinite, unlocked.lambda)
        end

        @testset "Radial endpoint derivatives and normalization" begin
            static = rmn(0, 0, 0, [1]; normalization = :static, derivatives = 4, precision_bits = bits)
            @test only(static.value)==1
            @test all(iszero(only(channel)) for channel in values(static)[2:end])
            for normalization in (:standard, :static)
                regular = rmn(0, 0, 1, [1]; normalization, derivatives = 4, precision_bits = bits)
                lambda = eigenvalue(0, 0, 1; precision_bits = 160)
                y, dy, ddy, third, fourth = only.(values(regular))
                @test isapprox(2dy, (lambda-1)*y; rtol = tolerance)
                @test isapprox(4ddy, (lambda-3)*dy-2y; rtol = tolerance)
                @test isapprox(6third, (lambda-7)*ddy-4dy-2y; rtol = tolerance)
                @test isapprox(8fourth, (lambda-13)*third-6ddy-6dy; rtol = tolerance)
            end
            irregular = rmn(0, 0, 0, [1]; kind = 2, normalization = :static,
                derivatives = 4, precision_bits = bits)
            @test all(isnan(real(only(channel))) for channel in values(irregular)[3:end])
            for kind in (1, 2)
                standard = rmn(0, 1, 5//4, [2]; kind, logderivative = true, precision_bits = bits)
                continued = rmn(0, 1, 5//4, [complex(2)]; kind, normalization = :static,
                    logderivative = true, precision_bits = bits)
                factor = BigFloat(5//4)^(kind==1 ? -1 : 2)
                @test isapprox(only(continued.value), factor*only(standard.value); rtol = tolerance)
                @test isapprox(only(continued.logderivative), only(standard.logderivative); rtol = tolerance)
            end
            @test isapprox(only(radial_wronskian(0, 0, 5//4, [2]; form = :normalized,
                precision_bits = bits)), 1; rtol = tolerance)
        end

        @testset "Angular origins, endpoint ratios, and Jacobians" begin
            origin = smn(0, 0, 0, [0]; kind = 2, derivatives = 4, precision_bits = bits)
            @test only.(values(origin))==(0, 1, 0, 2, 0)
            for kind in (1, 2)
                endpoints = smn(1, 1, 0, [-1, 1]; kind, logderivative = true, precision_bits = bits)
                @test endpoints.logderivative==(kind==1 ? [Inf, -Inf] : [-Inf, Inf])
            end
            derivative = jacobian_smn(0, 0, 0, [1//4]; precision_bits = bits)
            @test iszero(only(derivative.dvalue_dc)) && iszero(only(derivative.dderivative_dc))
            SW._with_swprecision(256) do
                x = one(SW._SWFloat)-ldexp(one(SW._SWFloat), -180)
                @test SW._precision_input_bits(x)==180
                expected = sqrt((1-x)*(1+x))
                near_endpoint = smn(1, 1, 0, [x]; precision_bits = bits)
                @test isapprox(abs(only(near_endpoint.value)), BigFloat(expected); rtol = tolerance)
                @test isapprox(abs(only(near_endpoint.derivative)), BigFloat(x/expected); rtol = tolerance)
            end
        end

        @testset "Scaled rounding and bounded coefficient refinement" begin
            SW._with_swprecision(192) do
                scaled = SW._precision_wave_output((; value = [SW._SWFloat("9.999")]), 8; scaled = true)
                @test only(scaled.value.mantissa)==1 && only(scaled.value.exponent)==1
                pair = SW._refine_coefficient_truncation(0, 0, 1, :prolate,
                    SW._mpfr_eigenvalue(0, 0, 1, :prolate, 192), SW._SWFloat(2)^(-80), 16, false)
                @test pair.converged && pair.terms==16
            end
            @test_throws r"quantity did not converge" SW._refine_quantity_at_precision(_ -> (; value = NaN), bits)
        end

        @testset "Imaginary parameters and endpoint phase anchors" begin
            SW._with_swprecision(192) do
                c = complex(SW._SWFloat(0), SW._SWFloat(1))
                pair = SW._solve_coefficient_eigenpair_at_precision(0, 1, c,
                    :prolate, :double, bits, nothing, 128, false; native_label = true)
                reference = SW._call_mpfr_batch(0, 1, 1, [1//4], :oblate, 1, true, 192)
                phase = SW._coefficient_phase_at_precision(pair)
                value = only(SW._evaluate_coefficient_vector(pair, pair.v, [1//4]).value)
                @test isapprox(phase*value, only(reference.value); rtol = tolerance)
                @test isapprox(pair.lambda, only(reference.eigenvalue); rtol = tolerance)
                evaluator = SW._mpfr_angular_phase_evaluator(0, 0, :prolate, :double, bits, 128)
                state = evaluator(complex(SW._SWFloat(1)), 1)
                @test state.endpoint_anchor && 0<state.endpoint_point<1
                @test length(state.profile)==10 && state.anchor==state.profile[end-1]
                endpoint_pair = SW._solve_coefficient_eigenpair_at_precision(0, 1, complex(1),
                    :prolate, :double, bits, nothing, 128, false; native_label = true)
                anchored = (; endpoint_pair..., mode_state = state, phase = typeof(pair.phase)(nothing))
                endpoint_phase = SW._coefficient_phase_at_precision(anchored)
                endpoint = only(SW._evaluate_coefficient_vector(anchored, anchored.v, [state.endpoint_point]).value)
                @test isapprox(endpoint_phase*endpoint, state.anchor; rtol = tolerance)
                unresolved = (; pair..., v = zero.(pair.v), phase = typeof(pair.phase)(nothing))
                @test_throws r"phase cannot be resolved" SW._coefficient_phase_at_precision(unresolved)
                @test unresolved.phase[]===nothing
            end
        end

        @testset "Refinement requires an independent confirmation" begin
            @test_throws r"working-precision refinement" SW._solve_coefficient_eigenpair_at_precision(
                0, 0, 0, :prolate, :double, bits, nothing, 16, false; max_refinements = 1)
            pair = SW._solve_coefficient_eigenpair_at_precision(
                0, 0, 0, :prolate, :double, bits, nothing, 16, false; max_refinements = 2)
            @test pair.converged && pair.lambda==0
            @test_throws r"radial precision did not converge" SW._refine_radial_values_at_precision(
                0, 0, 0, [2], :prolate, :double, 1, :static, bits, 16, 1, false; max_refinements = 1)
            radial = SW._refine_radial_values_at_precision(
                0, 0, 0, [2], :prolate, :double, 1, :static, bits, 16, 1, false; max_refinements = 2)
            @test only(radial.value)==1 && only(radial.derivative)==0
            @test_throws r"Angular values did not converge" SW._refine_angular_values_at_precision(
                0, 0, 0, [1//4], :prolate, :double, bits, 16, 1, false, 1; max_refinements = 1)
            angular = SW._refine_angular_values_at_precision(
                0, 0, 0, [1//4], :prolate, :double, bits, 16, 1, false, 1; max_refinements = 2)
            @test isapprox(only(angular.result.value), 1; rtol = tolerance)
            @test only(angular.result.derivative)==0
            @test_throws ArgumentError SW._refine_angular_values_at_precision(
                0, 0, 0, [1//4], :prolate, :double, bits, 16, 2, true, 1)
        end
    end
end

@testset "Optional MPFR JLL product discovery" begin
    code = raw"""
    using SpheroidalWaves
    SW = SpheroidalWaves
    JLL = SW.SpheroidalWaves_jll
    if !isdefined(JLL, :libspheroidal_batch_mpfr_path)
        path = SW.backend_library(; precision=:mpfr)
        @assert isfile(path)
        Core.eval(JLL, :(const libspheroidal_batch_mpfr_path = $path))
    end
    @assert SW._configure_backends_from_jll!()
    @assert SW.backend_library(; precision=:mpfr)==JLL.libspheroidal_batch_mpfr_path
    @assert isfinite(eigenvalue(0, 0, 1; precision_bits=64))
    """
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $code`
    @test success(addenv(command, "JULIA_NUM_THREADS" => "1"))
end
