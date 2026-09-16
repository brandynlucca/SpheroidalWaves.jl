using SpheroidalWaves,Test

@testset "Direct decaying Hankel waves and sensitivities" begin
    # Independent dense spectral perturbation sums and explicit finite Hankel
    # polynomials. Refinements at 100/140 digits and 96/128 terms agree to
    # more than 38 relative digits. References below use the larger calculation.
    references = (
        (:prolate,
         big"-211.8274401177637813110808080716533763897461301849761356",
         [big"-1.526822014336009454118764598602475850723497336614650773e-18",
          big"-5.158392139587137360274175832276312890054448985327741949e-351"],
         [big"3.380652574194445468410271922948619373768973584783163882e-17",
          big"1.033120395097099409107000407859459028442570529122093247e-349"],
         [big"-3.112581684701817299881484821959787495950722731849506588e-18",
          big"-2.065915171886190073802838431717392526186141849496115023e-349"],
         [big"6.737752316256361296950430799084451105560621886070181174e-17",
          big"4.132446919131259074862177708503548146987135202892203722e-348"]),
        (:oblate,
         big"100.7131752364551773694329255610312012453908688656851678",
         [big"-2.930662508835181925646146911699722400696736591728175713e-21",
          big"-3.802038016847395155618029612498744199455358714591274911e-351"],
         [big"5.531297121407866376425594049147470245448159550908816509e-20",
          big"7.611800331549399620260330117854328694240686039155645507e-350"],
         [big"-6.698498171318505423601671530361688643448815958971790325e-21",
          big"-1.523192881767901826737513419651278102789492009746770819e-349"],
         [big"1.238026299728219560347271417121012291768675071211830741e-19",
          big"3.045679466260037154721896946195389166838334005842770057e-348"]),
    )
    unpack(r)=r.mantissa.*BigFloat(10).^r.exponent
    for precision in (:double,:quad)
        lib=SpheroidalWaves.backend_library(;precision)
        (lib===nothing || !isfile(lib)) && continue
        T=precision===:quad ? BigFloat : Float64
        tolerance=precision===:quad ? big"1e-28" : big"1e-12"
        xs=T[2,40]
        for (spheroid,lambda,v,d,vc,dc) in references,kind in (3,4)
            sign=kind==3 ? 1 : -1
            c=complex(zero(T),sign*T(20))
            r=rmn(2,4,c,xs;spheroid,precision,kind,scaled=true,logderivative=true,second_derivative=true)
            @test eltype(r.value.mantissa)==Complex{T}
            @test r.value.exponent[2]==-351
            for i in eachindex(xs)
                @test unpack(r.value)[i]≈v[i] rtol=tolerance atol=0
                @test unpack(r.derivative)[i]≈d[i] rtol=tolerance atol=0
                @test r.logderivative[i]≈d[i]/v[i] rtol=tolerance
                x=BigFloat(xs[i])
                sigma=spheroid===:prolate ? 1 : -1
                second=-(2x*d[i]+(-400x^2-lambda-sigma*4/(x^2-sigma))*v[i])/(x^2-sigma)
                @test unpack(r.second_derivative)[i]≈second rtol=tolerance atol=0
            end
            j=jacobian_rmn(2,4,c,xs;spheroid,precision,kind)
            # Double outputs underflow at x=40; quad derivatives preserve them.
            for i in (precision===:quad ? (1:2) : (1:1))
                @test j.dvalue_dcreal[i]≈sign*im*vc[i] rtol=tolerance atol=0
                @test j.dderivative_dcreal[i]≈sign*im*dc[i] rtol=tolerance atol=0
            end
        end
    end
end
