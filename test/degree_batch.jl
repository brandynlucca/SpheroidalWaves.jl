using SpheroidalWaves,Test,Libdl

@testset "Shared real degree expansions" begin
    for precision in (:double,:quad),spheroid in (:prolate,:oblate)
        T=precision===:quad ? BigFloat : Float64
        tol=precision===:quad ? big"1e-25" : big"3e-12"
        c=T(1.25)
        points=T.([-.7,0,.3])
        # Start above m to check native degree indexing, and use odd m to
        # check the Condon-Shortley phase in both output formats.
        for (normalize,scaled) in ((false,false),(true,true))
            block=smn(1,2:3,c,points;spheroid,precision,normalize,scaled,logderivative=true,second_derivative=true)
            for (i,n) in enumerate(2:3)
                one=smn(1,n,c,points;spheroid,precision,normalize,scaled,logderivative=true,second_derivative=true)
                for field in (:value,:derivative,:second_derivative)
                    a,b=getproperty(block,field),getproperty(one,field)
                    if scaled
                        @test a.mantissa[:,i]≈b.mantissa rtol=tol atol=tol
                        @test a.exponent[:,i]==b.exponent
                    else
                        @test a[:,i]≈b rtol=tol atol=tol
                    end
                end
                @test all(isapprox(a,b;rtol=tol,atol=tol,nans=true) for (a,b) in zip(block.logderivative[:,i],one.logderivative))
            end
        end
        x=T.(spheroid===:prolate ? [1.2,2.] : [0.,.4,2.])
        for kind in 1:4,scaled in (false,true)
            block=rmn(1,2:4,c,x;spheroid,precision,kind,scaled,logderivative=true,second_derivative=true)
            for (i,n) in enumerate(2:4)
                one=rmn(1,n,c,x;spheroid,precision,kind,scaled,logderivative=true,second_derivative=true)
                for field in (:value,:derivative,:second_derivative)
                    a,b=getproperty(block,field),getproperty(one,field)
                    if scaled
                        @test a.mantissa[:,i]≈b.mantissa rtol=tol atol=tol
                        @test a.exponent[:,i]==b.exponent
                    else
                        @test a[:,i]≈b rtol=tol atol=tol
                    end
                end
                @test all(isapprox(a,b;rtol=tol,atol=tol,nans=true) for (a,b) in zip(block.logderivative[:,i],one.logderivative))
            end
        end
        # A mixed range keeps the refined low degrees and batches the rest.
        block=smn(0,0:3,T(3),T(.3);spheroid,precision,scaled=true)
        for (i,n) in enumerate(0:3)
            one=smn(0,n,T(3),T(.3);spheroid,precision,scaled=true)
            @test block.value.mantissa[:,i]≈one.value.mantissa rtol=tol
            @test block.value.exponent[:,i]==one.value.exponent
        end
        # Independent pair consistency, both geometries.
        r1=rmn(1,2:4,c,x;spheroid,precision)
        r2=rmn(1,2:4,c,x;spheroid,precision,kind=2)
        weight=c.*(x.^2 .+ (spheroid===:prolate ? -1 : 1))
        @test weight.*(r1.value.*r2.derivative-r1.derivative.*r2.value)≈ones(length(x),3) rtol=tol atol=tol
        @test size(smn(1,2:2,c,T(.3);spheroid,precision).value)==(1,1)
        symbol=precision===:quad ? :spheroidal_degrees_scaled_text : :spheroidal_degrees_scaled_double
        lib=SpheroidalWaves.backend_library(;precision)
        if Libdl.dlsym_e(SpheroidalWaves._require_backend_handle(lib),symbol)!=C_NULL
            @test SpheroidalWaves._shared_real_degree_range(1,2:4,c,points,spheroid,
                precision,:angular,0,true,false,false)!==nothing
        end
        # These values exceed Float64's exponent range. No reconstruction in
        # Float64 may occur before the scaled mantissas are returned.
        block=smn(200,200:201,T(1),T(.3);spheroid,precision,scaled=true)
        @test all(>(308),block.value.exponent)
        for (i,n) in enumerate(200:201)
            one=smn(200,n,T(1),T(.3);spheroid,precision,scaled=true)
            @test block.value.mantissa[:,i]≈one.value.mantissa rtol=tol
            @test block.value.exponent[:,i]==one.value.exponent
        end
    end
end
