using SpheroidalWaves,Test

@testset "Angular parameter derivatives for both kinds" begin
    # Independent finite differences of the full DLMF Ferrers-Q/P series,
    # including the negative-degree extension. Refining 95 to 115 digits,
    # 32 to 40 terms, and h=1e-10 to 1e-12 retains over 35 digits.
    # All reference values are inline.
    references = (
        (0,0,"1","0",:prolate,
         ("-0.1313518835856766918752924659743676345239","0"),
         ("-0.4752047408624118783574188646391051313358","0")),
        (1,1,"1","0",:prolate,
         ("-0.8434977614096064722355337593609288960688","0"),
         ("-2.870788837566549191218142759113917841296","0")),
        (1,2,"1.25","0",:prolate,
         ("-0.08017374880523985786470200285082474697146","0"),
         ("-0.5376226274474168135618373699166373460832","0")),
        (1,1,"1.25","0.2",:prolate,
         ("-0.8480663794480312512783796078823772145376","-1.906919009723626806024548779655968820397"),
         ("-2.885445842617955380528164280862569194937","-6.507038442270115185753488599328174746823")),
        (1,2,"1.25","0.2",:oblate,
         ("0.2179613780587120007502795149300736563302","0.05631418219550209522184540113541395714847"),
         ("0.4003472107347302543179468428970512511774","0.06424407974837104690466510092410146436666")),
    )
    for precision in (:double,:quad)
        lib=SpheroidalWaves.backend_library(;precision)
        (lib===nothing || !isfile(lib)) && continue
        T=precision===:quad ? BigFloat : Float64
        tolerance=precision===:quad ? big"1e-28" : big"2e-12"
        x=T(3)/10
        for (m,n,cr,ci,spheroid,v,d) in references
            c=ci=="0" ? parse(T,cr) : complex(parse(T,cr),parse(T,ci))
            j=jacobian_smn(m,n,c,[-x,x];precision,spheroid,kind=2,with_metadata=true)
            value=c isa Real ? j.dvalue_dc : j.dvalue_dcreal
            derivative=c isa Real ? j.dderivative_dc : j.dderivative_dcreal
            metadata=c isa Real ? j.metadata_value : j.metadata_value_dcreal
            @test value[2]≈complex(parse(BigFloat,v[1]),parse(BigFloat,v[2])) rtol=tolerance
            @test derivative[2]≈complex(parse(BigFloat,d[1]),parse(BigFloat,d[2])) rtol=tolerance
            @test value[1]==(-1)^(n-m+1)*value[2]
            @test derivative[1]==(-1)^(n-m)*derivative[2]
            @test eltype(value)==(c isa Real ? T : Complex{T})
            @test metadata.method===:differentiated_equation && metadata.step_used===nothing
            if !(c isa Real)
                @test j.dvalue_dcimag==im.*value
                @test j.dderivative_dcimag==im.*derivative
            end
        end
        for spheroid in (:prolate,:oblate)
            j=jacobian_smn(1,2,zero(T),[-x,x];precision,spheroid,kind=2)
            @test all(iszero,j.dvalue_dc) && all(iszero,j.dderivative_dc)
        end
        j=jacobian_smn(0,0,T(1),T[-1,1];precision,kind=2,with_metadata=true)
        @test all(isnan,j.dvalue_dc) && all(isnan,j.dderivative_dc)
        @test !j.metadata_value.finite_flag && j.metadata_value.conditioning_flag===:poor

        # Independent symmetric-matrix eigenvectors and differentiated Ferrers
        # sums at c=50. A tiny absolute tolerance would hide cancellation here.
        p=jacobian_smn(0,0,T(50),T[0,x];precision,spheroid=:oblate)
        @test p.dvalue_dc[1]≈big"-7.638376215356727332038888034783750088421157894e-21" rtol=tolerance atol=0
        @test iszero(p.dderivative_dc[1])
        @test p.dvalue_dc[2]≈big"-6.677179171731690380463396962857462465540847081e-15" rtol=tolerance atol=0
        @test p.dderivative_dc[2]≈big"-3.190030642668206332312178067136090394527345440e-13" rtol=tolerance atol=0
        unit=jacobian_smn(0,0,T(50),T[0,x];precision,spheroid=:oblate,normalize=true)
        @test unit.dvalue_dc≈p.dvalue_dc/sqrt(T(2)) rtol=tolerance
        odd=jacobian_smn(1,2,T(30),T[0,x];precision,spheroid=:oblate)
        @test iszero(odd.dvalue_dc[1])
        @test odd.dderivative_dc[1]≈big"8.981151898474642144202904101731327439038e-10" rtol=tolerance atol=0
        @test odd.dvalue_dc[2]≈big"5.213776262245525508155836340940084743715e-8" rtol=tolerance atol=0
        @test odd.dderivative_dc[2]≈big"1.387295596972673244090605137815319778311e-6" rtol=tolerance atol=0
    end
    lib=SpheroidalWaves.backend_library(;precision=:double)
    if lib!==nothing && isfile(lib)
        c=1.25+0.2im
        a=jacobian_smn(1,2,c,[0.3];kind=2,spheroid=:oblate)
        f=jacobian_smn(1,2,c,[0.3];kind=2,spheroid=:oblate,h=1e-4,adaptive=false)
        @test f.dvalue_dcreal≈a.dvalue_dcreal rtol=1e-7
        @test f.dderivative_dcimag≈a.dderivative_dcimag rtol=1e-7
    end
    @test_throws ArgumentError jacobian_smn(0,0,1.0,[0.3];kind=3)
    @test_throws ArgumentError jacobian_smn(0,0,1.0,[0.3];kind=2,normalize=true)
    @test_throws ArgumentError jacobian_smn(0,0,1.0,[0.3];kind=2,normalize=true,h=1e-4)
end
