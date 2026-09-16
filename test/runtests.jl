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

    @testset "Inline angular and radial references" begin
            wolfram_smn_cases = [
                ((0, 0, 0.0, 0.3), 1.0, 0.0),
                ((0, 1, 0.0, 0.3), 0.3, 1.0),
                ((1, 1, 0.3, 0.0), -1.001795214382128, 0.0),
                ((1, 1, 0.0, 0.3), -sqrt(0.91), 0.3 / sqrt(0.91)),
                ((1, 2, 0.5, 0.25), -0.7309295344004918, -2.7222906159495754),
                ((2, 3, 1.0, 0.5), 5.65036805385163, 3.4543263211124438),
                ((2, 5, 10.0, 0.6), 13.831303334796436, 21.40430253808833),
            ]
            for ((m, n, c, eta_wa), expected_value, expected_derivative) in wolfram_smn_cases
                result = smn(m, n, c, eta_wa; spheroid=:prolate, precision=:double)
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

            # Direct Wolfram|Alpha radial benchmarks.
            wolfram_rmn_kind1_cases = [
                ((0, 0, 0.5, 1.5), 0.9355129525869241),
                ((0, 1, 1.0, 2.0), 0.45603690333332372100494242600202824193682868459324),
                ((2, 3, 2.0, 1.8), 0.1681696391119482),
            ]
            for ((m, n, c, x_wa), expected_value) in wolfram_rmn_kind1_cases
                result = rmn(m, n, c, [x_wa]; spheroid=:prolate, precision=:double, kind=1)
                @test isapprox(real(result.value[1]), expected_value; atol=1e-11, rtol=1e-11)
                @test isapprox(imag(result.value[1]), 0.0; atol=1e-12, rtol=0.0)
            end

            wolfram_rmn_kind2_cases = [
                ((1, 2, 0.5, 1.3), -24.89681123745018),
            ]
            for ((m, n, c, x_wa), expected_value) in wolfram_rmn_kind2_cases
                result = rmn(m, n, c, [x_wa]; spheroid=:prolate, precision=:double, kind=2)
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
