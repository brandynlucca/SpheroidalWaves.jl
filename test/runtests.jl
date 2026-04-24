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

    @testset "Double Precision Numerical Paths" begin
        if !has_backend(:double)
            @info "Skipping double-precision numerical tests: backend library not available."
        else
            # Exact benchmark points (Wolfram-compatible spherical limit c = 0).
            eta = [-0.5, 0.25, 0.75]
            expected_p0 = [1.0, 1.0, 1.0]
            expected_p0d = [0.0, 0.0, 0.0]
            expected_p1 = [-0.5, 0.25, 0.75]
            expected_p1d = [1.0, 1.0, 1.0]
            expected_p2 = [-0.125, -0.40625, 0.34375]
            expected_p2d = [-1.5, 0.75, 2.25]

            s00 = smn(0, 0, 0.0, eta; spheroid=:prolate, precision=:double)
            assert_allclose(s00.value, expected_p0; atol=1e-10)
            assert_allclose(s00.derivative, expected_p0d; atol=1e-10)

            s01 = smn(0, 1, 0.0, eta; spheroid=:prolate, precision=:double)
            assert_allclose(s01.value, expected_p1; atol=1e-10)
            assert_allclose(s01.derivative, expected_p1d; atol=1e-10)

            s02 = smn(0, 2, 0.0, eta; spheroid=:prolate, precision=:double)
            assert_allclose(s02.value, expected_p2; atol=1e-9)
            assert_allclose(s02.derivative, expected_p2d; atol=1e-9)

            so0 = smn(0, 0, 0.0, eta; spheroid=:oblate, precision=:double)
            assert_allclose(so0.value, expected_p0; atol=1e-10)
            assert_allclose(so0.derivative, expected_p0d; atol=1e-10)

            so1 = smn(0, 1, 0.0, eta; spheroid=:oblate, precision=:double)
            assert_allclose(so1.value, expected_p1; atol=1e-10)
            assert_allclose(so1.derivative, expected_p1d; atol=1e-10)

            # Exact prolate radial benchmarks at c = 0, m = n = 0:
            # kind 1: P_0(x)=1, d/dx=0
            # kind 2: Q_0(x)=0.5*log((x+1)/(x-1)), d/dx=-1/(x^2-1)
            x = [1.2, 1.35]
            q0 = [1.1989476363991853, 0.9521187263273725]
            q0d = [-2.2727272727272730, -1.2158054711246200]

            r1 = rmn(0, 0, 0.0, x; spheroid=:prolate, precision=:double, kind=1)
            r2 = rmn(0, 0, 0.0, x; spheroid=:prolate, precision=:double, kind=2)
            r3 = rmn(0, 0, 0.0, x; spheroid=:prolate, precision=:double, kind=3)
            r4 = rmn(0, 0, 0.0, x; spheroid=:prolate, precision=:double, kind=4)
            w = radial_wronskian(0, 0, 0.0, x; spheroid=:prolate, precision=:double)

            assert_allclose(real.(r1.value), [1.0, 1.0]; atol=1e-10)
            assert_allclose(imag.(r1.value), [0.0, 0.0]; atol=1e-10)
            assert_allclose(real.(r1.derivative), [0.0, 0.0]; atol=1e-10)
            assert_allclose(imag.(r1.derivative), [0.0, 0.0]; atol=1e-10)

            assert_allclose(real.(r2.value), q0; atol=5e-9)
            assert_allclose(imag.(r2.value), [0.0, 0.0]; atol=1e-10)
            assert_allclose(real.(r2.derivative), q0d; atol=5e-9)
            assert_allclose(imag.(r2.derivative), [0.0, 0.0]; atol=1e-10)

            assert_allclose(real.(r3.value), [1.0, 1.0]; atol=1e-9)
            assert_allclose(imag.(r3.value), q0; atol=5e-9)
            assert_allclose(real.(r3.derivative), [0.0, 0.0]; atol=1e-9)
            assert_allclose(imag.(r3.derivative), q0d; atol=5e-9)

            assert_allclose(real.(r4.value), [1.0, 1.0]; atol=1e-9)
            assert_allclose(imag.(r4.value), -q0; atol=5e-9)
            assert_allclose(real.(r4.derivative), [0.0, 0.0]; atol=1e-9)
            assert_allclose(imag.(r4.derivative), -q0d; atol=5e-9)

            # Exact c=0 Wronskian identity: P_n Q_n' - P_n' Q_n = -1 / (x^2 - 1)
            w_expected = ComplexF64.(-1.0 ./ (x .^ 2 .- 1.0))
            assert_allclose(real.(w), real.(w_expected); atol=5e-9)
            assert_allclose(imag.(w), imag.(w_expected); atol=5e-9)

            # Complex-family c=0 benchmark equivalence to real Legendre limit.
            scp = smn(0, 2, 0.0 + 0.0im, eta; spheroid=:prolate, precision=:double)
            assert_allclose(real.(scp.value), expected_p2; atol=1e-9)
            assert_allclose(imag.(scp.value), [0.0, 0.0, 0.0]; atol=1e-10)
            assert_allclose(real.(scp.derivative), expected_p2d; atol=1e-9)
            assert_allclose(imag.(scp.derivative), [0.0, 0.0, 0.0]; atol=1e-10)

            sco = smn(0, 2, 0.0 + 0.0im, eta; spheroid=:oblate, precision=:double)
            assert_allclose(real.(sco.value), expected_p2; atol=1e-9)
            assert_allclose(imag.(sco.value), [0.0, 0.0, 0.0]; atol=1e-10)
            assert_allclose(real.(sco.derivative), expected_p2d; atol=1e-9)
            assert_allclose(imag.(sco.derivative), [0.0, 0.0, 0.0]; atol=1e-10)

            # Published upstream benchmark from Prolate_swf sample files:
            # profcndat.txt -> c = 200, x = 1.1, m = 0, eta = 0:0.2:1.0
            # profort20.txt / profort30.txt give literal fort.20 / fort.30 outputs.
            published_eta = [0.0, 0.2, 0.4]
            published_s00 = [2.82335544836906e0, 5.04086281350719e-2, 1.68883998486357e-7]
            published_s00d = [0.0, -2.05006402662351e0, -1.46818835366416e-5]

            sp00 = smn(0, 0, 200.0, published_eta; spheroid=:prolate, precision=:double, normalize=true)
            assert_allclose(sp00.value, published_s00; atol=1e-12)
            assert_allclose(sp00.derivative, published_s00d; atol=1e-12)

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

            # Accuracy API direct checks.
            acc_ang_c0 = accuracy(0, 2, 0.0, eta; spheroid=:prolate, precision=:double, target=:angular)
            acc_rad_c0 = accuracy(0, 0, 0.0, x; spheroid=:prolate, precision=:double, target=:radial, kind=1)
            @test acc_ang_c0 == fill(15, length(eta))
            @test acc_rad_c0 == fill(15, length(x))

            acc_ang_real = accuracy(0, 1, 200.0, [0.0, 0.2, 0.4]; spheroid=:prolate, precision=:double, target=:angular)
            acc_rad_real = accuracy(0, 1, 200.0, [1.1, 1.2]; spheroid=:prolate, precision=:double, target=:radial, kind=1)
            acc_ang_cplx = accuracy(0, 1, 200.0 + 0.1im, [0.0, 0.2, 0.4]; spheroid=:oblate, precision=:double, target=:angular)
            acc_rad_cplx = accuracy(0, 1, 200.0 + 0.1im, [1.2, 1.4]; spheroid=:oblate, precision=:double, target=:radial, kind=3)

            @test length(acc_ang_real) == 3 && all(>=(0), acc_ang_real)
            @test length(acc_rad_real) == 2 && all(>=(0), acc_rad_real)
            @test length(acc_ang_cplx) == 3 && all(>=(0), acc_ang_cplx)
            @test length(acc_rad_cplx) == 2 && all(>=(0), acc_rad_cplx)

            # Eigenvalue API: real/complex and prolate/oblate symbol routing.
            @test eigenvalue(0, 2, 0.0; spheroid=:prolate, precision=:double) == 6.0
            @test eigenvalue(0, 2, 0.0; spheroid=:oblate, precision=:double) == 6.0
            @test eigenvalue(0, 2, 0.0 + 0.0im; spheroid=:prolate, precision=:double) == 6.0 + 0.0im
            @test eigenvalue(0, 2, 0.0 + 0.0im; spheroid=:oblate, precision=:double) == 6.0 + 0.0im

            λp = eigenvalue(0, 1, 200.0; spheroid=:prolate, precision=:double)
            λo = eigenvalue(0, 1, 200.0; spheroid=:oblate, precision=:double)
            λcp = eigenvalue(0, 1, 200.0 + 0.1im; spheroid=:prolate, precision=:double)
            λco = eigenvalue(0, 1, 200.0 + 0.1im; spheroid=:oblate, precision=:double)
            @test isfinite(λp)
            @test isfinite(λo)
            @test isfinite(real(λcp)) && isfinite(imag(λcp))
            @test isfinite(real(λco)) && isfinite(imag(λco))

            # Nonzero-c cross-family consistency at purely real c.
            λpr = eigenvalue(0, 1, 200.0; spheroid=:prolate, precision=:double)
            λpc = eigenvalue(0, 1, 200.0 + 0.0im; spheroid=:prolate, precision=:double)
            λor = eigenvalue(0, 1, 200.0; spheroid=:oblate, precision=:double)
            λoc = eigenvalue(0, 1, 200.0 + 0.0im; spheroid=:oblate, precision=:double)
            @test isapprox(real(λpc), λpr; atol=1e-10, rtol=1e-10)
            @test isapprox(imag(λpc), 0.0; atol=1e-10, rtol=0.0)
            @test isapprox(real(λoc), λor; atol=1e-10, rtol=1e-10)
            @test isapprox(imag(λoc), 0.0; atol=1e-10, rtol=0.0)

            # Root finding: recover c from a target eigenvalue.
            c_true_p = 80.0
            λ_target_p = eigenvalue(0, 1, c_true_p; spheroid=:prolate, precision=:double)
            root_p = find_c_for_eigenvalue(0, 1, λ_target_p;
                                           bracket=(60.0, 100.0),
                                           spheroid=:prolate,
                                           precision=:double)
            @test root_p.converged
            @test isapprox(root_p.c, c_true_p; atol=1e-6, rtol=1e-8)
            @test abs(root_p.residual) <= 1e-6

            c_true_o = 60.0
            λ_target_o = eigenvalue(0, 1, c_true_o; spheroid=:oblate, precision=:double)
            root_o = find_c_for_eigenvalue(0, 1, λ_target_o;
                                           bracket=(40.0, 80.0),
                                           spheroid=:oblate,
                                           precision=:double,
                                           use_jacobian=false)
            @test root_o.converged
            @test isapprox(root_o.c, c_true_o; atol=1e-6, rtol=1e-8)

            @test_throws ErrorException find_c_for_eigenvalue(0, 1, 1e12;
                                                               bracket=(0.0, 1.0),
                                                               spheroid=:prolate,
                                                               precision=:double)

            # Jacobians: compare against independent finite differences.
            h = 1e-6
            λjac = jacobian_eigen(0, 1, 200.0; spheroid=:prolate, precision=:double, h=h, adaptive=false)
            λfd = (eigenvalue(0, 1, 200.0 + h; spheroid=:prolate, precision=:double) -
                   eigenvalue(0, 1, 200.0 - h; spheroid=:prolate, precision=:double)) / (2h)
            @test isapprox(λjac, λfd; atol=1e-8, rtol=1e-8)

                λjac_meta = jacobian_eigen(0, 1, 200.0; spheroid=:prolate, precision=:double, h=h, with_metadata=true, adaptive=false)
                 @test isapprox(λjac_meta.derivative, λjac; atol=1e-8, rtol=1e-8)
                 @test λjac_meta.metadata.step_used > 0.0
                 @test λjac_meta.metadata.step_used <= h
                 @test λjac_meta.metadata.relative_change_when_halving_step >= 0.0
                 @test λjac_meta.metadata.finite_flag
                 @test λjac_meta.metadata.conditioning_flag in (:good, :warning, :poor)
                 @test λjac_meta.metadata.suggested_action in (:accept, :retry_smaller_h, :use_quad)

            λcjac = jacobian_eigen(0, 1, 200.0 + 0.1im; spheroid=:oblate, precision=:double, h=h, adaptive=false)
            λcfd_re = (eigenvalue(0, 1, 200.0 + h + 0.1im; spheroid=:oblate, precision=:double) -
                       eigenvalue(0, 1, 200.0 - h + 0.1im; spheroid=:oblate, precision=:double)) / (2h)
            λcfd_im = (eigenvalue(0, 1, 200.0 + (0.1 + h)im; spheroid=:oblate, precision=:double) -
                       eigenvalue(0, 1, 200.0 + (0.1 - h)im; spheroid=:oblate, precision=:double)) / (2h)
            @test isapprox(λcjac.d_dcreal, λcfd_re; atol=1e-8, rtol=1e-8)
            @test isapprox(λcjac.d_dcimag, λcfd_im; atol=1e-8, rtol=1e-8)

                λcjac_meta = jacobian_eigen(0, 1, 200.0 + 0.1im; spheroid=:oblate, precision=:double, h=h, with_metadata=true, adaptive=false)
                 @test isapprox(λcjac_meta.d_dcreal, λcjac.d_dcreal; atol=1e-8, rtol=1e-8)
                 @test isapprox(λcjac_meta.d_dcimag, λcjac.d_dcimag; atol=1e-8, rtol=1e-8)
                 @test λcjac_meta.metadata_dcreal.finite_flag
                 @test λcjac_meta.metadata_dcimag.finite_flag
                 @test λcjac_meta.metadata_dcreal.conditioning_flag in (:good, :warning, :poor)
                 @test λcjac_meta.metadata_dcimag.conditioning_flag in (:good, :warning, :poor)

            eta_j = [0.0, 0.2]
            sjac = jacobian_smn(0, 1, 200.0, eta_j; spheroid=:prolate, precision=:double, h=h, adaptive=false)
            sp = smn(0, 1, 200.0 + h, eta_j; spheroid=:prolate, precision=:double)
            sm = smn(0, 1, 200.0 - h, eta_j; spheroid=:prolate, precision=:double)
            sfd_val = (sp.value .- sm.value) ./ (2h)
            sfd_der = (sp.derivative .- sm.derivative) ./ (2h)
            assert_allclose(sjac.dvalue_dc, sfd_val; atol=1e-8)
            assert_allclose(sjac.dderivative_dc, sfd_der; atol=1e-8)

                sjac_meta = jacobian_smn(0, 1, 200.0, eta_j; spheroid=:prolate, precision=:double, h=h, with_metadata=true, adaptive=false)
                 assert_allclose(sjac_meta.dvalue_dc, sjac.dvalue_dc; atol=1e-8)
                 assert_allclose(sjac_meta.dderivative_dc, sjac.dderivative_dc; atol=1e-8)
                 @test sjac_meta.metadata_value.finite_flag
                 @test sjac_meta.metadata_derivative.finite_flag

            x_j = [1.1, 1.2]
            rjac = jacobian_rmn(0, 1, 200.0, x_j; spheroid=:prolate, precision=:double, kind=1, h=h, adaptive=false)
            rp = rmn(0, 1, 200.0 + h, x_j; spheroid=:prolate, precision=:double, kind=1)
            rm = rmn(0, 1, 200.0 - h, x_j; spheroid=:prolate, precision=:double, kind=1)
            rfd_val = (rp.value .- rm.value) ./ (2h)
            rfd_der = (rp.derivative .- rm.derivative) ./ (2h)
            assert_allclose(real.(rjac.dvalue_dc), real.(rfd_val); atol=1e-8)
            assert_allclose(imag.(rjac.dvalue_dc), imag.(rfd_val); atol=1e-8)
            assert_allclose(real.(rjac.dderivative_dc), real.(rfd_der); atol=1e-8)
            assert_allclose(imag.(rjac.dderivative_dc), imag.(rfd_der); atol=1e-8)

            rjac_meta = jacobian_rmn(0, 1, 200.0, x_j; spheroid=:prolate, precision=:double, kind=1, h=h, with_metadata=true, adaptive=false)
            assert_allclose(real.(rjac_meta.dvalue_dc), real.(rjac.dvalue_dc); atol=1e-8)
            assert_allclose(imag.(rjac_meta.dvalue_dc), imag.(rjac.dvalue_dc); atol=1e-8)
            @test rjac_meta.metadata_value.finite_flag
            @test rjac_meta.metadata_derivative.finite_flag
        end
    end

    @testset "Quad Precision Numerical Paths" begin
        if !has_backend(:quad)
            @info "Skipping quad-precision numerical tests: backend library not available."
        else
            eta = [-0.5, 0.25, 0.75]
            expected_p2 = [-0.125, -0.40625, 0.34375]
            expected_p2d = [-1.5, 0.75, 2.25]

            sq = smn(0, 2, 0.0, eta; spheroid=:prolate, precision=:quad)
            assert_allclose(sq.value, expected_p2; atol=1e-9)
            assert_allclose(sq.derivative, expected_p2d; atol=1e-9)

            sqc = smn(0, 2, 0.0 + 0.0im, eta; spheroid=:oblate, precision=:quad)
            assert_allclose(real.(sqc.value), expected_p2; atol=1e-9)
            assert_allclose(imag.(sqc.value), [0.0, 0.0, 0.0]; atol=1e-10)

            @test eigenvalue(0, 2, 0.0; spheroid=:prolate, precision=:quad) == 6.0
            @test eigenvalue(0, 2, 0.0 + 0.0im; spheroid=:oblate, precision=:quad) == 6.0 + 0.0im

            accq_ang_c0 = accuracy(0, 2, 0.0, eta; spheroid=:prolate, precision=:quad, target=:angular)
            accq_rad_c0 = accuracy(0, 0, 0.0, [1.2, 1.35]; spheroid=:prolate, precision=:quad, target=:radial, kind=1)
            @test accq_ang_c0 == fill(30, length(eta))
            @test accq_rad_c0 == fill(30, 2)

            wq = radial_wronskian(0, 0, 0.0, [1.2, 1.35]; spheroid=:prolate, precision=:quad)
            wq_expected = ComplexF64.(-1.0 ./ ([1.2, 1.35] .^ 2 .- 1.0))
            assert_allclose(real.(wq), real.(wq_expected); atol=5e-9)
            assert_allclose(imag.(wq), imag.(wq_expected); atol=5e-9)

            if has_backend(:double)
                λd = eigenvalue(0, 1, 200.0; spheroid=:prolate, precision=:double)
                λq = eigenvalue(0, 1, 200.0; spheroid=:prolate, precision=:quad)
                @test isapprox(λq, λd; rtol=1e-6, atol=1e-6)

                λcd = eigenvalue(0, 1, 200.0 + 0.1im; spheroid=:oblate, precision=:double)
                λcq = eigenvalue(0, 1, 200.0 + 0.1im; spheroid=:oblate, precision=:quad)
                @test isapprox(real(λcq), real(λcd); rtol=1e-6, atol=1e-6)
                @test isapprox(imag(λcq), imag(λcd); rtol=1e-6, atol=1e-6)

                hq = 1e-6
                jq = jacobian_eigen(0, 1, 200.0; spheroid=:prolate, precision=:quad, h=hq, adaptive=false)
                jd = jacobian_eigen(0, 1, 200.0; spheroid=:prolate, precision=:double, h=hq, adaptive=false)
                @test isapprox(jq, jd; rtol=1e-4, atol=1e-4)

                c_true_q = 80.0
                λ_target_q = eigenvalue(0, 1, c_true_q; spheroid=:prolate, precision=:quad)
                root_q = find_c_for_eigenvalue(0, 1, λ_target_q;
                                               bracket=(60.0, 100.0),
                                               spheroid=:prolate,
                                               precision=:quad)
                @test root_q.converged
                @test isapprox(root_q.c, c_true_q; atol=1e-6, rtol=1e-8)
            end
        end
    end
end

