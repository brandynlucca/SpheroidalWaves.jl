using SpheroidalWaves, Test

@testset "Small values and imaginary bandwidth" begin
    setprecision(BigFloat,256) do
        # Independent Legendre eigensolves at 512/640 bits and 192/256 terms
        # agreed through the digits retained here; neither used the backend.
        expected = big"7.1783988291048889209520759337759276912940375371729226688757616389946757807e-49"
        slope = big"-2.9420866932171962019988954259535047772049356026568309326530906613819548479e-46"
        # Independently summed outgoing spherical Hankel waves, initialized
        # with exp(i*z), at 100/130 decimal digits and 96/120 terms.
        references = (
            (:prolate,big"-1.964709742392256121324629526258086287074741262368975984e-10",
                      big"2.152367602399809247697497954705026813163816880885369716e-9"),
            (:oblate, big"-1.140537841478724748209688801620624775874000182876878836e-11",
                      big"1.082403572638871876979451050094353588649533756952487463e-10"),
        )
        unpack(r) = r.mantissa .* BigFloat(10) .^ r.exponent
        for precision in (:double,:quad)
            lib = SpheroidalWaves.backend_library(;precision)
            if lib === nothing || !isfile(lib)
                @info "Skipping numerical domain regressions: backend unavailable" precision
                continue
            end
            tolerance = precision===:quad ? big"1e-28" : big"1e-12"
            s = smn(0,0,200,big"0.9";precision,scaled=true,logderivative=true,second_derivative=true)
            @test only(unpack(s.value)) ≈ expected rtol=tolerance atol=0
            @test only(unpack(s.derivative)) ≈ slope rtol=tolerance atol=0
            @test only(s.logderivative) ≈ slope/expected rtol=tolerance
            @test all(isfinite,s.second_derivative.mantissa)
            degrees = smn(0,0:1,200,big"0.9";precision)
            @test degrees.value[1,1] ≈ expected rtol=tolerance atol=0
            @test accuracy(0,0,200,[big"0.9"];precision,target=:angular) == [-1]

            for (spheroid,value,derivative) in references
                opposite = spheroid===:prolate ? :oblate : :prolate
                @test eigenvalue(0,0,10im;spheroid,precision) ==
                      eigenvalue(0,0,10;spheroid=opposite,precision)
                @test smn(0,0,10im,big"0.3";spheroid,precision).value ==
                      smn(0,0,10,big"0.3";spheroid=opposite,precision).value
                for (c,kind) in ((10im,3),(-10im,4))
                    r = rmn(0,0,c,2;spheroid,precision,kind,scaled=true,logderivative=true)
                    @test only(unpack(r.value)) ≈ value rtol=tolerance atol=0
                    @test only(unpack(r.derivative)) ≈ derivative rtol=tolerance atol=0
                    @test only(r.logderivative) ≈ derivative/value rtol=tolerance
                end
            end
        end

        libraries = [SpheroidalWaves.backend_library(;precision) for precision in (:double,:quad)]
        all(lib -> lib !== nothing && isfile(lib),libraries) || return
        # A double output request must not move a supplied interior BigFloat
        # coordinate onto a singular endpoint in the analytic evaluators.
        for (spheroid,kind,c) in ((:oblate,1,5),(:prolate,2,1))
            results = [smn(1,1,c,1-big"1e-40";spheroid,kind,precision,
                          scaled=true,logderivative=true) for precision in (:double,:quad)]
            for r in results
                @test all(isfinite,r.value.mantissa)
                @test all(isfinite,r.derivative.mantissa)
                @test all(isfinite,r.logderivative)
            end
            @test unpack(results[1].value) ≈ unpack(results[2].value) rtol=1e-12
            @test results[1].logderivative ≈ results[2].logderivative rtol=1e-12
        end
        r = jacobian_rmn(1,2,1.25,[1+big"1e-40"];kind=2,precision=:double)
        # Refined differentiated-equation reference at the same coordinate.
        @test only(r.dvalue_dc) ≈ big"6.61972135479573379523347505281974739772e20" rtol=1e-12
        @test only(r.dderivative_dc) ≈ big"-3.30986067739786689761673752640987369884e60" rtol=1e-12
    end
end
