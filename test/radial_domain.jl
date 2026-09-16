using SpheroidalWaves, Test

@testset "Radial functions exclude exact zero parameter" begin
    # Reject before calling a backend, for every radial kind and input family.
    # The former prolate m=0 fallback returned unrelated Legendre P/Q values.
    for spheroid in (:prolate, :oblate), precision in (:double, :quad),
        c in (0, 0.0, -0.0, 0.0 + 0.0im, big"0.0", complex(big"0.0", big"0.0"))
        for kind in 1:4
            @test_throws DomainError rmn(0, 1, c, [2.0]; spheroid, precision, kind)
            @test_throws DomainError accuracy(0, 1, c, [2.0]; spheroid, precision, kind, target=:radial)
        end
        @test_throws DomainError rmn(1, 1, c, [2.0]; spheroid, precision)
        for (m, n) in ((0, 0), (0, 2), (1, 2))
            @test_throws DomainError rmn(m, n, c, [2.0]; spheroid, precision)
        end
        @test_throws DomainError rmn(0, 0:2, c, [2.0]; spheroid, precision)
        @test_throws DomainError radial_wronskian(0, 1, c, [2.0]; spheroid, precision)
        # A finite-difference derivative must reject an undefined center too.
        @test_throws DomainError jacobian_rmn(0, 1, c, [2.0]; spheroid, precision)
    end

    # A singular-coordinate error must not conceal the parameter-domain error.
    err = try
        rmn(0, 0, 0.0, [1.0])
    catch e
        e
    end
    @test err isa DomainError
    @test occursin("c != 0", sprint(showerror, err))
end

@testset "Small nonzero radial parameters remain supported" begin
    # Independent signed reference values, with decimal inputs preserved
    # for the high-precision second-kind comparisons.
    for precision in (:double, :quad)
        lib = SpheroidalWaves.backend_library(; precision)
        if lib === nothing || !isfile(lib)
            @info "Skipping nonzero radial references: backend unavailable." precision
            continue
        end
        first = rmn(0, 1, big"0.1", [2]; precision, kind=1)
        @test only(first.value) ≈ 0.06642698122211762 rtol=5e-15
        tolerance = precision === :quad ? big"1e-27" : big"5e-15"
        for (c, expected) in [
            (big"0.01", big"-2958.897611846449350102731791001293876590218151081837314"),
            (big"0.001", big"-295837.3950030077336871204221250744159129119120049634542")]
            second = rmn(0, 1, c, [2]; precision, kind=2)
            @test only(second.value) ≈ expected rtol=tolerance
        end
    end
end
