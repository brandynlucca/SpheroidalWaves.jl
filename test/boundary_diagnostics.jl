using SpheroidalWaves,Test

@testset "Large bandwidth tolerance and second-kind residuals" begin
    # At c=375, exp(-2c) underflows in Float64 even though the requested
    # coefficient tolerance is meant to be formed at the guarded precision.
    # Wider-input reference values are inline to avoid rebuilding the large
    # quad expansions during every normal test run.
    for (spheroid,value,ve,derivative,de) in (
        (:prolate,big"1.5256488534902150635105594465576701032672",-7,
                  big"-1.7954837073168552408315797203151455810388",-5),
        (:oblate,big"4.1896917085417494071462469474848952172303",-113,
                 big"1.5679082330112565415126277769813336590425",-110))
        ordinary=smn(0,0,375.,big".3";spheroid,precision=:double,scaled=true)
        @test only(ordinary.value.mantissa)≈value rtol=1e-13
        @test only(ordinary.value.exponent)==ve
        @test only(ordinary.derivative.mantissa)≈derivative rtol=1e-13
        @test only(ordinary.derivative.exponent)==de
    end
    for spheroid in (:prolate,:oblate),kind in (1,2)
        diagnostic=SpheroidalWaves._spheroidal_residual(0,0,big"1",big".3";
            spheroid,kind,precision=:quad,h=big"1e-5")
        actual=smn(0,0,big"1",big".3";spheroid,kind,precision=:quad,second_derivative=true)
        @test diagnostic.derivative≈actual.derivative rtol=big"1e-20"
        @test diagnostic.second_derivative≈actual.second_derivative rtol=big"1e-20"
        @test maximum(diagnostic.relative_residual)<big"1e-20"
    end
end
