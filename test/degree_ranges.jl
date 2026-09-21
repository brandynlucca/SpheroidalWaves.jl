using SpheroidalWaves, Test, Libdl
@testset "Degree ranges preserve wave functions" begin
    # Exercise multiple coordinates, a range starting above m, both norms,
    # angular endpoints, and all four radial kinds.
    for (m, degrees, c) in ((0, 0:6, 3.0), (2, 4:9, 11.0))
        eta = [-1.0, -0.7, 0.0, 0.3, 1.0]
        x = [1.02, 1.3]
        for normalize in (false, true)
            block = smn(m, degrees, c, eta; precision=:quad, normalize)
            @test size(block.value) == (length(eta), length(degrees))
            @test eltype(block.value) == BigFloat
            for (i, n) in enumerate(degrees)
                single = smn(m, n, c, eta; precision=:quad, normalize)
                @test block.value[:, i] ≈ single.value rtol=1e-24
                @test block.derivative[:, i] ≈ single.derivative rtol=1e-24
            end
        end
        for kind in 1:4
            block = rmn(m, degrees, c, x; precision=:quad, kind)
            @test eltype(block.value) == Complex{BigFloat}
            for (i, n) in enumerate(degrees)
                single = rmn(m, n, c, x; precision=:quad, kind)
                @test block.value[:, i] ≈ single.value rtol=1e-24
                @test block.derivative[:, i] ≈ single.derivative rtol=1e-24
            end
        end
        r1 = rmn(m, degrees, c, x; precision=:quad, kind=1)
        r2 = rmn(m, degrees, c, x; precision=:quad, kind=2)
        wronskian = r1.value .* r2.derivative .- r1.derivative .* r2.value
        @test wronskian ≈ repeat(1 ./ (c .* (x.^2 .- 1)), 1, length(degrees)) rtol=1e-12
    end
    for spheroid in (:prolate, :oblate), precision in (:double, :quad)
        a = smn(1, 2:3, 2.0, [0.3]; spheroid, precision)
        r = rmn(1, 2:3, 2.0, [1.2]; spheroid, precision)
        for (i,n) in enumerate(2:3)
            @test a.value[:, i] ≈ smn(1, n, 2.0, [0.3]; spheroid, precision).value rtol=(precision === :quad ? 1e-24 : 1e-12)
            @test r.value[:, i] ≈ rmn(1, n, 2.0, [1.2]; spheroid, precision).value rtol=(precision === :quad ? 1e-24 : 1e-12)
        end
    end
    @test smn(0, 0:2, 0.0, 0.2).value == hcat((smn(0,n,0.0,0.2).value for n in 0:2)...)
    @test_throws ArgumentError smn(2, 1:4, 1.0, [0.0])
    @test_throws ArgumentError smn(0, 2:1, 1.0, [0.0])
    @test_throws ArgumentError smn(0, 0:2, Inf, [0.0])
    @test_throws ArgumentError smn(0, 0:2, 1.0, [1.1])
    @test_throws ArgumentError rmn(0, 0:2, 1.0, Float64[])
    @test_throws ArgumentError rmn(0, 0:2, 1.0, [1.1]; kind=5)
end
