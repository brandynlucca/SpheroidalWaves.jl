using SpheroidalWaves
using Test

@testset "SpheroidalWaves.jl" begin
    assert_allclose(actual, expected; atol=1e-10, rtol=0.0) = begin
        @test length(actual) == length(expected)
        @test all(isapprox.(actual, expected; atol=atol, rtol=rtol))
    end

    has_backend(precision::Symbol) = begin
        lib = SpheroidalWaves.backend_library(precision=precision)
        lib !== nothing && isfile(lib)
    end

    @testset "Public API Surface" begin
        exported = names(SpheroidalWaves)
        @test :smn in exported
        @test :rmn in exported
        @test :radial_wronskian in exported
        @test :accuracy in exported
        @test :eigenvalue in exported
        @test :eigenvalue_sweep in exported
        @test :jacobian_eigen in exported
        @test :jacobian_smn in exported
        @test :jacobian_rmn in exported
        @test :find_c_for_eigenvalue in exported
        @test !(:set_backend_library! in exported)
        @test !(:backend_library in exported)
    end

    @testset "Local Backend Configuration Does Not Evaluate Generated Code" begin
        @test !isdefined(SpheroidalWaves, :SPHEROIDAL_BATCH_LIBRARY_DOUBLE)
        @test !isdefined(SpheroidalWaves, :SPHEROIDAL_BATCH_LIBRARY_QUAD)
    end

    @testset "Backend configuration and ABI failures" begin
        SW = SpheroidalWaves
        saved = copy(SW._backend_libraries)
        try
            @test SW._configure_backends_from_jll!()
            @test all(has_backend, (:double, :quad))
            @test !SW._set_backend_from_candidate(nothing, :double, "test")
            @test_logs (:warn, r"not a string") !SW._set_backend_from_candidate(1, :double, "test")
            missing = joinpath(tempdir(), "spheroidalwaves-nonexistent-library")
            @test_logs (:warn, r"does not exist") !SW._set_backend_from_candidate(missing, :double, "test")
            withenv(SW._ENV_BACKEND_DOUBLE => saved[:double], SW._ENV_BACKEND_QUAD => saved[:quad]) do
                @test SW._configure_backends_from_env!()
                @test SW.backend_library() == saved[:double]
                @test SW.backend_library(precision=:quad) == saved[:quad]
            end
            withenv(SW._ENV_BACKEND_DOUBLE => nothing, SW._ENV_BACKEND_QUAD => nothing) do
                @test !SW._configure_backends_from_env!()
            end
            SW._backend_libraries[:double] = nothing
            @test_throws r"No backend library configured" SW._require_backend_library(:double)
            SW._configure_backends_from_local_build!()
            candidate = SW.backend_library()
            @test candidate === nothing || isfile(candidate)
            @test SW._check_scalar_status(Cint(0)) === nothing
            @test SW._check_vector_status(Cint[0, 0]) === nothing
            @test_throws r"status 3" SW._check_scalar_status(Cint(3))
            @test_throws r"non-zero status" SW._check_vector_status(Cint[0, 1])
            @test SW._has_required_quad_abi(saved[:quad])
            @test !SW._has_required_quad_abi(saved[:double])
            @test_throws r"Rebuild" SW._quad_symbol_pointer(saved[:quad], :nonexistent_spheroidal_symbol)
            for (kernel, expected) in ((:NT, "wave.dll"), (:Darwin, "libwave.dylib"),
                                        (:Linux, "libwave.so"), (:FreeBSD, "libwave.so"))
                @test SW._backend_filename("wave", kernel) == expected
            end
            mktempdir() do directory
                SW._backend_libraries[:double] = SW._backend_libraries[:quad] = nothing
                @test !SW._configure_backends_from_local_build!(directory)
                @test SW.backend_library() === nothing
                @test SW.backend_library(precision=:quad) === nothing
                # Discovery only checks filenames. Never load this placeholder.
                candidate = joinpath(directory, SW._backend_filename("spheroidal_batch_double"))
                touch(candidate)
                @test SW._configure_backends_from_local_build!(directory)
                @test SW.backend_library() == candidate
                @test SW.backend_library(precision=:quad) === nothing
            end
        finally
            merge!(SW._backend_libraries, saved)
        end
    end

    @testset "Backend initialization order and failure reporting" begin
        SW = SpheroidalWaves
        calls = Symbol[]
        SW._initialize_backends!(() -> push!(calls, :jll),
                                 () -> push!(calls, :local), () -> push!(calls, :env))
        @test calls == [:jll, :local, :env]
        empty!(calls)
        @test_logs (:warn, r"Failed to configure backend libraries.*injected initialization failure") SW._initialize_backends!(
            () -> error("injected initialization failure"),
            () -> push!(calls, :local), () -> push!(calls, :env))
        @test isempty(calls)
    end

    @testset "Backend Does Not Create fort.60" begin
        for precision in (:double, :quad)
            if !has_backend(precision)
                @info "Skipping fort.60 regression test: backend library not available." precision
                continue
            end

            mktempdir() do dir
                cd(dir) do
                    rmn(0, 0, 13.995744383559012, [1.1547005383792515];
                        spheroid=:prolate, precision=precision, kind=2)
                end
                @test !isfile(joinpath(dir, "fort.60"))
            end
        end
    end

    @testset "Concurrent Backend Calls" begin
        if Threads.nthreads() == 1
            @info "Skipping concurrent backend regression test: Julia has one thread."
        else
            for precision in (:double, :quad)
                if !has_backend(precision)
                    @info "Skipping concurrent backend regression test: backend library not available." precision
                    continue
                end

                cases = [(
                    m=mod(i, 2),
                    n=mod(i, 2) + 2,
                    c=is_complex ? 2.0 + 0.05im + 0.02i : 2.0 + 0.02i,
                    eta=[-0.7 + 0.01i, 0.1 + 0.002i, 0.65 - 0.003i],
                    x=spheroid === :prolate ? [1.1 + 0.001i, 1.3 + 0.002i] : [0.1 + 0.001i, 0.3 + 0.002i],
                    spheroid=spheroid,
                    kind=is_complex ? 1 : 2,
                ) for spheroid in (:prolate, :oblate), is_complex in (false, true), i in 1:2]
                serial = [
                    (
                        angular=smn(case.m, case.n, case.c, case.eta;
                                    spheroid=case.spheroid, precision=precision),
                        radial=rmn(case.m, case.n, case.c, case.x;
                                   spheroid=case.spheroid, precision=precision, kind=case.kind),
                    )
                    for case in cases
                ]

                # One overlapping batch covers both orders, geometries and
                # parameter types; repeating it starts the same native work again.
                tasks = [
                    Threads.@spawn begin
                        angular = smn(case.m, case.n, case.c, case.eta;
                                      spheroid=case.spheroid, precision=precision)
                        radial = rmn(case.m, case.n, case.c, case.x;
                                     spheroid=case.spheroid, precision=precision, kind=case.kind)
                        (; angular, radial)
                    end
                    for case in cases
                ]
                concurrent = fetch.(tasks)

                for (actual, expected) in zip(concurrent, serial)
                    @test actual.angular.value == expected.angular.value
                    @test actual.angular.derivative == expected.angular.derivative
                    @test actual.radial.value == expected.radial.value
                    @test actual.radial.derivative == expected.radial.derivative
                end
            end
        end
    end

    @testset "Argument Validation" begin
        @test_throws ErrorException smn(0, 0, 1.0, [0.0]; precision=:bad)
        @test_throws ErrorException rmn(0, 0, 1.0, [1.1]; precision=:bad)
        @test_throws ErrorException eigenvalue(0, 0, 1.0; precision=:bad)
        @test_throws ErrorException jacobian_eigen(0, 0, 1.0; precision=:bad)
        @test_throws ErrorException smn(0, 0, 1.0, [0.0]; spheroid=:bad)
        @test_throws ErrorException rmn(0, 0, 1.0, [1.1]; spheroid=:bad)
        @test_throws ErrorException eigenvalue(0, 0, 1.0; spheroid=:bad)
        @test_throws ErrorException jacobian_eigen(0, 0, 1.0; h=0.0)
        @test_throws ErrorException jacobian_eigen(0, 0, 1.0; rtol=0.0)
        @test_throws ErrorException jacobian_eigen(0, 0, 1.0; atol=0.0)
        @test_throws ErrorException find_c_for_eigenvalue(0, 0, 1.0; bracket=(1.0, 0.0))
        @test_throws ErrorException find_c_for_eigenvalue(0, 0, 1.0; bracket=(0.0, 1.0), atol=0.0)
        @test_throws ErrorException find_c_for_eigenvalue(0, 0, 1.0; bracket=(0.0, 1.0), rtol=0.0)
        @test_throws ErrorException find_c_for_eigenvalue(0, 0, 1.0; bracket=(0.0, 1.0), maxiter=0)
        @test_throws ErrorException accuracy(0, 0, 1.0, [1.1]; target=:bad)
        @test_throws ErrorException eigenvalue_sweep(0, 1, Float64[])
        @test_throws ErrorException eigenvalue_sweep(0, 1, [1.0, 1.0])
        @test_throws ErrorException eigenvalue_sweep(0, 1, [0.0, 1.0, 0.5])
        @test_throws ErrorException eigenvalue_sweep(0, 1, [0.0, 1.0]; branch_window=-1)
    end

    @testset "Scalar smn Convenience Overload" begin
        vector_result = smn(0, 2, 0.0, [0.5]; spheroid=:prolate, precision=:double)
        scalar_result = smn(0, 2, 0.0, 0.5; spheroid=:prolate, precision=:double)

        @test scalar_result.value == vector_result.value
        @test scalar_result.derivative == vector_result.derivative
    end

    @testset "Eigenvalue Continuation Branch Lock" begin
        # Synthetic evaluator that swaps n=2 and n=3 branches at c=1.0 only.
        # The continuation predictor-corrector should stay on the original smooth
        # curve by locally selecting the best branch candidate.
        synthetic_eval = function (m, n, c)
            _ = m
            base = 100.0 * n + 10.0 * c
            if isapprox(c, 1.0; atol=0.0, rtol=0.0)
                if n == 2
                    return 100.0 * 3 + 10.0 * c
                elseif n == 3
                    return 100.0 * 2 + 10.0 * c
                end
            end
            return base
        end

        cgrid = [0.0, 1.0, 2.0, 3.0]
        raw = eigenvalue_sweep(0, 2, cgrid;
                               branch_lock=false,
                               use_jacobian_predictor=false,
                               evaluator=synthetic_eval)
        locked = eigenvalue_sweep(0, 2, cgrid;
                                  branch_lock=true,
                                  branch_window=1,
                                  use_jacobian_predictor=false,
                                  evaluator=synthetic_eval)

        @test raw.lambda == [200.0, 310.0, 220.0, 230.0]
        @test locked.lambda == [200.0, 210.0, 220.0, 230.0]
        @test raw.selected_n == [2, 2, 2, 2]
        @test locked.selected_n == [2, 3, 2, 2]
        @test locked.switched_branch == [false, true, true, false]
    end

    @testset "Inverse characteristic values" begin
        for (spheroid,use_jacobian) in ((:prolate,true),(:oblate,false))
            target = eigenvalue(0,1,2.0;spheroid)
            root = find_c_for_eigenvalue(0,1,target;spheroid,use_jacobian,bracket=(1.0,3.0))
            @test root.converged
            @test root.c ≈ 2.0 atol=1e-8
            @test abs(root.residual)<1e-8
        end
    end

    @testset "Inverse characteristic endpoints and failures" begin
        target = eigenvalue(0, 1, 2.0)
        for bracket in ((2.0, 3.0), (1.0, 2.0))
            root = find_c_for_eigenvalue(0, 1, target; bracket)
            @test root.converged && root.c == 2 && root.iterations == 0
            @test root.method === :endpoint && root.residual == 0
        end
        failed = find_c_for_eigenvalue(0, 1, target; bracket=(1.0, 4.0), maxiter=1,
                                      atol=1e-15, rtol=1e-15)
        @test !failed.converged && failed.method === :maxiter && failed.iterations == 1
        @test failed.residual ≈ eigenvalue(0, 1, failed.c)-target
        @test failed.bracket[1] <= failed.c <= failed.bracket[2]
        @test_throws r"straddle" find_c_for_eigenvalue(0, 1, target; bracket=(3.0, 4.0))
        @test_throws r"non-finite residual" find_c_for_eigenvalue(0, 1, NaN; bracket=(1.0, 4.0))
        @test_throws ErrorException find_c_for_eigenvalue(0, 1, target; bracket=(1, 4), spheroid=:bad)
        @test_throws ArgumentError find_c_for_eigenvalue(0, 1, target; bracket=(1, 4), form=:log)
        # The secant proposal for x^2-1/4 on [0,1] is x=1/4. Simulate a
        # failed evaluation there and require recovery at the midpoint root.
        samples = Float64[]
        residual(x) = (push!(samples, x); x == 0.25 ? NaN : x^2-0.25)
        unexpected_jacobian(x) = error("derivative-free inversion called its Jacobian")
        recovered = SpheroidalWaves._find_separation_root(residual, unexpected_jacobian,
            0.25, 0., 1., 1e-12, 1e-10, 10, false)
        @test 0.25 in samples
        @test recovered.converged && recovered.c == 0.5 && recovered.residual == 0
        @test recovered.method === :bisection && recovered.iterations == 1
    end

    @testset "Continuation grids and unavailable candidates" begin
        @test_throws ErrorException eigenvalue_sweep(0, 1, [0.0]; spheroid=:bad)
        @test_throws ErrorException eigenvalue_sweep(0, 1, [0.0]; evaluator=(m,n,c) -> NaN)
        @test_throws ErrorException eigenvalue_sweep(0, 1, [0.0]; evaluator=(m,n,c) -> 1im)
        @test_throws ErrorException eigenvalue_sweep(0, 1, [0., 1.];
            evaluator=(m,n,c) -> c == 0 ? 2. : NaN, use_jacobian_predictor=false)
        @test_throws ErrorException eigenvalue_sweep(0, 1, [0., 1.];
            evaluator=(m,n,c) -> c == 0 ? 2. : 1im, use_jacobian_predictor=false)
        single = eigenvalue_sweep(0, 1, [0.0])
        @test single.lambda == [2.] && single.selected_n == [1]
        @test single.switched_branch == [false]
        for precision in (:double, :quad)
            grid = [1., 0.75, 0.5]
            sweep = eigenvalue_sweep(0, 1, grid; precision)
            @test sweep.lambda ≈ [eigenvalue(0, 1, c; precision) for c in grid]
            @test sweep.selected_n == [1, 1, 1]
            path = [1.0+0.0im, 1.0+0.1im, 1.0+0.1im, 1.0+0.0im]
            for branch_lock in (false, true)
                sweep = eigenvalue_sweep(0, 1, path; precision, branch_lock)
                @test sweep.lambda ≈ [eigenvalue(0, 1, c; precision) for c in path]
                @test sweep.c == path && sweep.selected_n == ones(Int, 4)
                @test !any(sweep.switched_branch)
            end
        end
        for (m,n,path,kwargs) in ((-1,0,[1im],(;)), (0,0,ComplexF64[],(;)),
                                  (0,0,[Inf+0im],(;)), (0,0,[1im],(;spheroid=:bad)),
                                  (0,0,[1im],(;branch_window=-1)))
            @test_throws ErrorException eigenvalue_sweep(m,n,path;kwargs...)
        end
        @test_throws ArgumentError eigenvalue_sweep(0, 0, [1im]; evaluator=identity)
    end

    @testset "Inline angular and radial references" begin
            angular_reference_cases = [
                ((0, 0, 0.0, 0.3), 1.0, 0.0),
                ((0, 1, 0.0, 0.3), 0.3, 1.0),
                ((1, 1, 0.3, 0.0), -1.001795214382128, 0.0),
                ((1, 1, 0.0, 0.3), -sqrt(0.91), 0.3 / sqrt(0.91)),
                ((1, 2, 0.5, 0.25), -0.7309295344004918, -2.7222906159495754),
                ((2, 3, 1.0, 0.5), 5.65036805385163, 3.4543263211124438),
                ((2, 5, 10.0, 0.6), 13.831303334796436, 21.40430253808833),
            ]
            for ((m, n, c, eta), expected_value, expected_derivative) in angular_reference_cases
                result = smn(m, n, c, eta; spheroid=:prolate, precision=:double)
                @test isapprox(result.value[1], expected_value; atol=1e-11, rtol=1e-11)
                @test isapprox(result.derivative[1], expected_derivative; atol=1e-11, rtol=1e-11)
            end

            x_prolate = [1.1]
            prolate_r1_l0 = [-6.32691894914518e-3]
            prolate_r1d_l0 = [-1.47014223025242e0]
            prolate_r2_l0 = [3.10920865064746e-3]
            prolate_r2d_l0 = [-3.04074463797897e0]
            prolate_r1_l1 = [-4.46044642130715e-3]
            prolate_r1d_l1 = [-2.59985144640758e0]
            prolate_r2_l1 = [5.47799916812461e-3]
            prolate_r2d_l1 = [-2.14497358451685e0]

            rp1_l0 = rmn(0, 0, 200.0, x_prolate; spheroid=:prolate, precision=:double, kind=1)
            rp2_l0 = rmn(0, 0, 200.0, x_prolate; spheroid=:prolate, precision=:double, kind=2)
            rp1_l1 = rmn(0, 1, 200.0, x_prolate; spheroid=:prolate, precision=:double, kind=1)
            rp2_l1 = rmn(0, 1, 200.0, x_prolate; spheroid=:prolate, precision=:double, kind=2)

            assert_allclose(real.(rp1_l0.value), prolate_r1_l0; atol=1e-12)
            assert_allclose(imag.(rp1_l0.value), [0.0]; atol=1e-12)
            assert_allclose(real.(rp1_l0.derivative), prolate_r1d_l0; atol=1e-12)
            assert_allclose(imag.(rp1_l0.derivative), [0.0]; atol=1e-12)
            assert_allclose(real.(rp2_l0.value), prolate_r2_l0; atol=1e-12)
            assert_allclose(imag.(rp2_l0.value), [0.0]; atol=1e-12)
            assert_allclose(real.(rp2_l0.derivative), prolate_r2d_l0; atol=1e-12)
            assert_allclose(imag.(rp2_l0.derivative), [0.0]; atol=1e-12)

            assert_allclose(real.(rp1_l1.value), prolate_r1_l1; atol=1e-12)
            assert_allclose(imag.(rp1_l1.value), [0.0]; atol=1e-12)
            assert_allclose(real.(rp1_l1.derivative), prolate_r1d_l1; atol=1e-12)
            assert_allclose(imag.(rp1_l1.derivative), [0.0]; atol=1e-12)
            assert_allclose(real.(rp2_l1.value), prolate_r2_l1; atol=1e-12)
            assert_allclose(imag.(rp2_l1.value), [0.0]; atol=1e-12)
            assert_allclose(real.(rp2_l1.derivative), prolate_r2d_l1; atol=1e-12)
            assert_allclose(imag.(rp2_l1.derivative), [0.0]; atol=1e-12)

            # Independent radial reference values.
            radial_kind1_reference_cases = [
                ((0, 0, 0.5, 1.5), 0.9355129525869241),
                ((0, 1, 1.0, 2.0), 0.45603690333332372100494242600202824193682868459324),
                ((2, 3, 2.0, 1.8), 0.1681696391119482),
            ]
            for ((m, n, c, x), expected_value) in radial_kind1_reference_cases
                result = rmn(m, n, c, [x]; spheroid=:prolate, precision=:double, kind=1)
                @test isapprox(real(result.value[1]), expected_value; atol=1e-11, rtol=1e-11)
                @test isapprox(imag(result.value[1]), 0.0; atol=1e-12, rtol=0.0)
            end

            radial_kind2_reference_cases = [
                ((1, 2, 0.5, 1.3), -24.89681123745018),
            ]
            for ((m, n, c, x), expected_value) in radial_kind2_reference_cases
                result = rmn(m, n, c, [x]; spheroid=:prolate, precision=:double, kind=2)
                @test isapprox(real(result.value[1]), expected_value; atol=1e-10, rtol=1e-11)
                @test isapprox(imag(result.value[1]), 0.0; atol=1e-12, rtol=0.0)
            end

            # Published upstream benchmark from Oblate_swf sample files:
            # oblfcndat.txt -> c = 500, x = 0.2, m = 0, eta = 0:0.2:1.0
            # oblfort20.txt gives literal fort.20 outputs.
            x_oblate = [0.2]
            oblate_r1_l0 = [1.46470793668965e-3]
            oblate_r1d_l0 = [6.51950929637809e-1]
            oblate_r2_l0 = [-1.30698205147780e-3]
            oblate_r2d_l0 = [7.31196119559882e-1]
            oblate_r1_l1 = [-1.30698205147780e-3]
            oblate_r1d_l1 = [7.31196119559882e-1]
            oblate_r2_l1 = [-1.46470793668965e-3]
            oblate_r2d_l1 = [-6.51950929637809e-1]

            ro1_l0 = rmn(0, 0, 500.0, x_oblate; spheroid=:oblate, precision=:double, kind=1)
            ro2_l0 = rmn(0, 0, 500.0, x_oblate; spheroid=:oblate, precision=:double, kind=2)
            ro1_l1 = rmn(0, 1, 500.0, x_oblate; spheroid=:oblate, precision=:double, kind=1)
            ro2_l1 = rmn(0, 1, 500.0, x_oblate; spheroid=:oblate, precision=:double, kind=2)

            assert_allclose(real.(ro1_l0.value), oblate_r1_l0; atol=1e-12)
            assert_allclose(imag.(ro1_l0.value), [0.0]; atol=1e-12)
            assert_allclose(real.(ro1_l0.derivative), oblate_r1d_l0; atol=1e-12)
            assert_allclose(imag.(ro1_l0.derivative), [0.0]; atol=1e-12)
            assert_allclose(real.(ro2_l0.value), oblate_r2_l0; atol=1e-12)
            assert_allclose(imag.(ro2_l0.value), [0.0]; atol=1e-12)
            assert_allclose(real.(ro2_l0.derivative), oblate_r2d_l0; atol=1e-12)
            assert_allclose(imag.(ro2_l0.derivative), [0.0]; atol=1e-12)

            assert_allclose(real.(ro1_l1.value), oblate_r1_l1; atol=1e-12)
            assert_allclose(imag.(ro1_l1.value), [0.0]; atol=1e-12)
            assert_allclose(real.(ro1_l1.derivative), oblate_r1d_l1; atol=1e-12)
            assert_allclose(imag.(ro1_l1.derivative), [0.0]; atol=1e-12)
            assert_allclose(real.(ro2_l1.value), oblate_r2_l1; atol=1e-12)
            assert_allclose(imag.(ro2_l1.value), [0.0]; atol=1e-12)
            assert_allclose(real.(ro2_l1.derivative), oblate_r2d_l1; atol=1e-12)
            assert_allclose(imag.(ro2_l1.derivative), [0.0]; atol=1e-12)

            @test_throws ErrorException smn(0, 1, 1.0, [1.2]; spheroid=:prolate, precision=:double)
            @test_throws ErrorException rmn(0, 1, 1.0, [0.9]; spheroid=:prolate, precision=:double, kind=1)
            @test_throws ErrorException rmn(0, 1, 1.0, [-0.1]; spheroid=:oblate, precision=:double, kind=1)
            @test_throws ErrorException rmn(0, 1, 1.0 + 0.1im, [0.9]; spheroid=:prolate, precision=:double, kind=1)
            @test_throws ErrorException rmn(0, 1, 1.0 + 0.1im, [-0.1]; spheroid=:oblate, precision=:double, kind=1)

    end
end

include("degree_ranges.jl")
include("degree_batch.jl")
include("angular_conventions.jl")

include("radial_domain.jl")
include("quad_precision.jl")

include("accuracy.jl")
include("numerical_regressions.jl")
include("angular_precision.jl")
include("parameter_derivatives.jl")
include("angular_second_kind.jl")
include("oblate_precision.jl")
include("angular_parameter_derivatives.jl")
include("radial_parameter_derivatives.jl")
include("numerical_domain.jl")
include("boundary_diagnostics.jl")
include("radial_hankel.jl")
include("integral_eigenvalues.jl")
include("integral_sensitivities.jl")
