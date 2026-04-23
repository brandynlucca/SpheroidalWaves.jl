using SpheroidalWaveFunctions
using Test

@testset "SpheroidalWaveFunctions.jl" begin
    assert_allclose(actual, expected; atol=1e-10, rtol=0.0) = begin
        @test length(actual) == length(expected)
        @test all(isapprox.(actual, expected; atol=atol, rtol=rtol))
    end

    has_backend(precision::Symbol) = begin
        lib = SpheroidalWaveFunctions.backend_library(precision=precision)
        lib !== nothing && isfile(lib)
    end

    @testset "Public API Surface" begin
        exported = names(SpheroidalWaveFunctions)
        @test :smn in exported
        @test :rmn in exported
        @test !(:set_backend_library! in exported)
        @test !(:backend_library in exported)
    end

    @testset "Argument Validation" begin
        @test_throws ErrorException smn(0, 0, 1.0, [0.0]; precision=:bad)
        @test_throws ErrorException rmn(0, 0, 1.0, [1.1]; precision=:bad)
        @test_throws ErrorException smn(0, 0, 1.0, [0.0]; spheroid=:bad)
        @test_throws ErrorException rmn(0, 0, 1.0, [1.1]; spheroid=:bad)
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
        end
    end
end
