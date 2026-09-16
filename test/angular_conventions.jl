using SpheroidalWaves, Test

@testset "Angular spherical limit and Condon–Shortley phase" begin
    # Closed-form Ferrers polynomials, independently of the implementation's
    # recurrence. Both normalization choices use the same phase.
    for precision in (:double, :quad), spheroid in (:prolate, :oblate),
        complex_c in (false, true), normalize in (false, true)
        T = precision === :quad ? BigFloat : Float64
        eta = T.([-0.5, 0.0, 0.25])
        u = 1 .- eta.^2
        c = complex_c ? complex(zero(T)) : zero(T)
        cases = (
            (0, 2, (3eta.^2 .- 1) ./ 2, 3eta),
            (1, 1, -sqrt.(u), eta ./ sqrt.(u)),
            (1, 2, -3eta .* sqrt.(u), -3 .* (1 .- 2eta.^2) ./ sqrt.(u)),
            (1, 3, -(3//2) .* (5eta.^2 .- 1) .* sqrt.(u), (3//2) .* eta .* (15eta.^2 .- 11) ./ sqrt.(u)),
            (1, 4, -(5//2) .* (7eta.^3 .- 3eta) .* sqrt.(u), (5//2) .* (28eta.^4 .- 27eta.^2 .+ 3) ./ sqrt.(u)),
            (2, 2, 3u, -6eta),
            (2, 3, 15eta .* u, 15 .* (1 .- 3eta.^2)),
            (3, 3, -15u .* sqrt.(u), 45eta .* sqrt.(u)),
        )
        for (m, n, values, derivatives) in cases
            norm_squared = T(2) / (2n + 1) * T(factorial(n + m)) / factorial(n - m)
            scale = normalize ? inv(sqrt(norm_squared)) : one(T)
            result = smn(m, n, c, eta; precision, spheroid, normalize)
            expected_type = complex_c ? Complex{T} : T
            @test eltype(result.value) == expected_type
            @test eltype(result.derivative) == expected_type
            @test result.value ≈ scale .* values rtol=64eps(T) atol=64eps(T)
            @test result.derivative ≈ scale .* derivatives rtol=64eps(T) atol=64eps(T)
        end
    end

    # Exact endpoint limits: m=1 really has divergent derivatives; the other
    # listed cases have finite limits and must not produce 0/0 NaNs.
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)
        s = smn(0, 2, 0.0, [-1.0, 1.0]; precision, spheroid)
        @test s.value == [1, 1]
        @test s.derivative == [-3, 3]
        @test smn(1, 1, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [-Inf, Inf]
        @test smn(1, 2, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [Inf, Inf]
        @test smn(2, 2, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [6, -6]
        @test smn(3, 3, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [0, 0]
        @test smn(2, 2, 0.0, [-1.0, 1.0]; precision, spheroid, normalize=true).derivative ≈ [6, -6] ./ sqrt(48 / 5)
    end

    setprecision(BigFloat, 256) do
        x = big"0.123456789012345678901234567890123456789"
        result = smn(1, 1, big"0", [x]; precision=:quad)
        @test abs(only(result.value) + sqrt(1 - x^2)) < big"1e-70"
        @test abs(only(result.derivative) - x / sqrt(1 - x^2)) < big"1e-70"
        # Unit normalization stays representable even when unnormalized P_m^m
        # would overflow double precision.
        m = 201
        expected = -sqrt(BigFloat(2m + 1) * binomial(big(2m), m) / big(2)^(2m + 1))
        normalized = smn(m, m, 0.0, 0.0; normalize=true)
        @test only(normalized.value) ≈ Float64(expected) rtol=1e-13
        @test only(normalized.derivative) == 0
    end

    @test_throws ErrorException smn(-1, 1, 0.0, 0.3)
    @test_throws ErrorException smn(2, 1, 0.0, 0.3)
    @test_throws ErrorException smn(1, 1, 0.0, 1.1)
    @test_throws ErrorException smn(1, 1, 0.0, NaN)
end

@testset "Angular phase across native backends and degree ranges" begin
    for precision in (:double, :quad)
        lib = SpheroidalWaves.backend_library(; precision)
        if lib === nothing || !isfile(lib)
            @info "Skipping native angular phase tests: backend unavailable." precision
            continue
        end
        # Independent reference for m=1, n=1, c=0.3, eta=0.
        value = only(smn(1, 1, 0.3, 0.0; precision).value)
        @test value ≈ -1.001795214382128 rtol=1e-14

        for spheroid in (:prolate, :oblate), complex_c in (false, true), normalize in (false, true)
            c = complex_c ? 0.3 + 0.0im : 0.3
            eta = [-0.25, 0.0, 0.25]
            for m in (1, 2)
                block = smn(m, m:(m + 3), c, eta; precision, spheroid, normalize)
                for (j, n) in enumerate(m:(m + 3))
                    single = smn(m, n, c, eta; precision, spheroid, normalize)
                    @test block.value[:, j] ≈ single.value rtol=1e-12 atol=1e-14
                    @test block.derivative[:, j] ≈ single.derivative rtol=1e-12 atol=1e-14
                    # DLMF 30.4(i), the phase specified after Eq. 30.4.1.
                    if iseven(n - m)
                        @test sign(real(block.value[2, j])) == (-1)^((n + m) ÷ 2)
                    else
                        @test sign(real(block.derivative[2, j])) == (-1)^((n + m - 1) ÷ 2)
                    end
                end
            end
            # Continuity to the closed-form spherical limit, independently of
            # whether the backend uses real or complex arithmetic.
            small_c = complex_c ? 1e-3 + 0.0im : 1e-3
            near = smn(1, 2, small_c, eta; precision, spheroid, normalize)
            exact = smn(1, 2, zero(small_c), eta; precision, spheroid, normalize)
            @test near.value ≈ exact.value rtol=1e-6 atol=1e-12
            @test near.derivative ≈ exact.derivative rtol=1e-6 atol=1e-12
        end

        # Approach zero through both sides of arg(c)=pi/4, where the old
        # complex oblate normalization selected opposite signs for n=m+2.
        for spheroid in (:prolate, :oblate), normalize in (false, true),
            c in (1e-3 + 1e-4im, 1e-4 + 1e-3im), n in (3, 4)
            eta = [-0.25, 0.25] # Phase probes are independent of this batch.
            near = smn(1, n, c, eta; precision, spheroid, normalize)
            exact = smn(1, n, zero(c), eta; precision, spheroid, normalize)
            @test near.value ≈ exact.value rtol=1e-6 atol=1e-12
            @test near.derivative ≈ exact.derivative rtol=1e-6 atol=1e-12
            with_anchor = smn(1, n, c, [0.0; eta]; precision, spheroid, normalize)
            @test near.value ≈ with_anchor.value[2:end] rtol=1e-13
            @test near.derivative ≈ with_anchor.derivative[2:end] rtol=1e-13
        end

        # The phase must propagate to the parameter Jacobian: the negative
        # S_11(c, 0) decreases locally around this reference point.
        jac = jacobian_smn(1, 1, 0.3, [0.0]; precision, h=1e-4, adaptive=false)
        @test only(jac.dvalue_dc) < 0
        @test only(jac.dderivative_dc) == 0
    end
end
