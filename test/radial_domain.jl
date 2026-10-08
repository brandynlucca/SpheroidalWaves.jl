using SpheroidalWaves, Test

@testset "Mixed radial boundary and interior batches" begin
    SW = SpheroidalWaves
    for precision in (:double, :quad), kind in 2:4
        tolerance = precision === :quad ? big"1e-28" : 1e-12
        for scaled in (false, true), points in ([1.], [1., 2.])
            result = rmn(0, 1, 1.25, points; precision, kind, scaled, second_derivative=true)
            unpack(v) = scaled ? v.mantissa .* BigFloat(10).^v.exponent : v
            @test isnan(first(unpack(result.value)))
            @test isnan(first(unpack(result.derivative)))
            @test isnan(first(unpack(result.second_derivative)))
            if length(points) == 2
                interior = rmn(0, 1, 1.25, 2; precision, kind, second_derivative=true)
                for field in keys(interior)
                    @test unpack(getproperty(result, field))[2] ≈ only(getproperty(interior, field)) rtol=tolerance
                end
            end
        end
        if precision === :quad
            native = SW._call_real_rmn(:psms, 0, 1, 1.25, [1., 2.]; precision, kind, with_accuracy=true)
            @test native.accuracy[1] == -1
            @test native.accuracy[2] == only(SW._call_real_rmn(:psms, 0, 1, 1.25, [2.]; precision, kind, with_accuracy=true).accuracy)
        end
    end
end

@testset "Complex native and analytic radial transfers" begin
    SW = SpheroidalWaves
    for spheroid in (:prolate, :oblate), precision in (:double, :quad)
        prefix = spheroid === :prolate ? :cprolate : :coblate
        tolerance = precision === :quad ? big"1e-27" : 1e-11
        for c in (10im, 1.25+0.1im), kind in 3:4
            direct = rmn(0, 1, c, 2; spheroid, precision, kind)
            native = SW._call_complex_rmn_raw(prefix, 0, 1, c, [2]; precision, kind)
            @test native.value ≈ direct.value rtol=tolerance
            @test native.derivative ≈ direct.derivative rtol=tolerance
            if precision === :quad && !iszero(real(c))
                quad = SW._complex_quad_radial(prefix, 0, 1, c, [2], kind)
                @test quad.value ≈ direct.value rtol=tolerance
                @test quad.derivative ≈ direct.derivative rtol=tolerance
            end
            if iszero(real(c))
                scaled = SW._scaled_native_values(0, 1, c, [2], spheroid, precision, :radial, kind)
                @test scaled.value ≈ direct.value rtol=tolerance
            end
        end
    end
    @test_throws DomainError SW._scaled_native_values(0, 1, 1.25+0.1im, [1], :prolate, :double, :radial, 1)
end

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
