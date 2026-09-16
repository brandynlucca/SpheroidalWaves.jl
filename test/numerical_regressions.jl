using SpheroidalWaves, Test

@testset "Focused numerical regressions" begin
    for precision in (:double,:quad)
        lib = SpheroidalWaves.backend_library(;precision)
        (lib === nothing || !isfile(lib)) && continue
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-27" : big"2e-11"

        # Independent spectral value, embedded here so CI need not reconstruct
        # a high-precision reference or execute the extensive path matrices.
        c = Complex{T}(2+3im)
        expected_lambda = complex(big"1.8526077784286968122645450083859036836748134683954",
                                  big"3.1017112824466201071716399658076138048679053435607")
        expected_s = complex(big"1.4672095477135660718169187632852425554599094354553",
                             big"0.02551066455758449011154602903404097070828353999752")
        @test eigenvalue(0,0,c;precision) ≈ expected_lambda rtol=tolerance
        @test only(smn(0,0,c,T(3)/10;precision).value) ≈ expected_s rtol=tolerance
        @test only(radial_wronskian(0,0,c,T(2);precision,form=:normalized)) ≈ 1 rtol=100tolerance

        for spheroid in (:prolate,:oblate)
            parameter = T(5)/4
            r = rmn(1,2,parameter,T(2);precision,spheroid,scaled=true,logderivative=true)
            ordinary = rmn(1,2,parameter,T(2);precision,spheroid)
            @test r.value.mantissa.*BigFloat(10).^r.value.exponent ≈ ordinary.value rtol=tolerance
            @test r.logderivative ≈ ordinary.derivative./ordinary.value rtol=tolerance
            @test only(radial_wronskian(1,2,parameter,T(2);precision,spheroid,form=:error)) < tolerance
        end

        # Exact spherical polynomials independently check coordinate derivatives
        # and zeros, without solving a numerical reference eigenproblem.
        s = smn(0,3,zero(T),T(1)/3;precision,second_derivative=true)
        @test only(s.second_derivative) ≈ T(5) rtol=tolerance
        roots = SpheroidalWaves.angular_zeros(0,3,zero(T);precision)
        a = sqrt(T(3)/5)
        @test roots ≈ T[-a,0,a] rtol=tolerance atol=tolerance
    end

    # Original exact-decimal radial boundary regressions, kept in the normal
    # suite even though the extensive boundary integration tests are local.
    lib = SpheroidalWaves.backend_library(;precision=:quad)
    if lib !== nothing && isfile(lib)
        for (m,n,value) in ((0,1,big"-28.19279152071340627099348388765367655558239286251609457"),
                            (1,2,big"-2846358.502301398239201169080317823638197684953289129961"))
            @test only(rmn(m,n,big"1.25",big"1.000000000001";precision=:quad,kind=2).value) ≈ value rtol=big"1e-27"
        end
    end

    s = smn(200,200,0.0,0.3;scaled=true,logderivative=true)
    @test only(s.value.exponent) > 308
    @test only(s.logderivative) ≈ -200*0.3/(1-0.3^2) rtol=2e-13
end
