using SpheroidalWaves,Test

@testset "Radial analytic parameter derivatives" begin
    for precision in (:double,:quad)
        lib=SpheroidalWaves.backend_library(;precision)
        (lib===nothing || !isfile(lib)) && continue
        T=precision===:quad ? BigFloat : Float64
        tolerance=precision===:quad ? big"1e-27" : big"2e-11"
        for spheroid in (:prolate,:oblate), complex_parameter in (false,true)
            c=complex_parameter ? complex(T(5)/4,T(1)/5) : T(5)/4
            xs=spheroid===:prolate ? T[T(5)/4,2] : T[0,2]
            a=jacobian_rmn(1,2,c,xs;precision,spheroid,kind=1)
            b=jacobian_rmn(1,2,c,xs;precision,spheroid,kind=2)
            av=complex_parameter ? a.dvalue_dcreal : a.dvalue_dc
            ad=complex_parameter ? a.dderivative_dcreal : a.dderivative_dc
            bv=complex_parameter ? b.dvalue_dcreal : b.dvalue_dc
            bd=complex_parameter ? b.dderivative_dcreal : b.dderivative_dc
            r=rmn(1,2,c,xs;precision,spheroid,kind=1)
            s=rmn(1,2,c,xs;precision,spheroid,kind=2)
            w=r.value.*s.derivative-r.derivative.*s.value
            wc=av.*s.derivative+r.value.*bd-ad.*s.value-r.derivative.*bv
            # Differentiate c*(x^2-sigma)*W=1. Native function values are
            # independent of the differentiated Bessel/ODE evaluator.
            sigma=spheroid===:prolate ? 1 : -1
            @test maximum(abs,(xs.^2 .-sigma).*(c.*wc+w))<tolerance
            @test eltype(av)==Complex{T}
            for kind in (3,4)
                h=jacobian_rmn(1,2,c,xs;precision,spheroid,kind)
                v=complex_parameter ? h.dvalue_dcreal : h.dvalue_dc
                d=complex_parameter ? h.dderivative_dcreal : h.dderivative_dc
                sign=kind==3 ? im : -im
                @test v≈av+sign*bv rtol=tolerance
                @test d≈ad+sign*bd rtol=tolerance
                if complex_parameter
                    @test h.dvalue_dcimag==im.*v
                    @test h.dderivative_dcimag==im.*d
                end
            end
        end
        for m in (0,1,2,3)
            endpoint=jacobian_rmn(m,m+1,T(5)/4,T[1];precision)
            if m==1
                @test isinf(only(endpoint.dderivative_dc))
            else
                @test all(isfinite,endpoint.dderivative_dc)
            end
            m>0 && @test iszero(only(endpoint.dvalue_dc))
        end
        for kind in (2,3,4)
            singular=jacobian_rmn(0,1,T(5)/4,T[1];precision,kind,with_metadata=true)
            @test all(isnan,singular.dvalue_dc) && !singular.metadata_value.finite_flag
        end
        @test_throws DomainError jacobian_rmn(0,0,complex(T(1),T(1)/10),T[1];precision)
    end
    lib=SpheroidalWaves.backend_library(;precision=:quad)
    if lib!==nothing && isfile(lib)
        # Independent 448-bit Bessel/Ferrers coefficients, coordinate Taylor
        # continuation, and a five-point parameter stencil. Increasing the
        # basis/order and reducing h preserved more than 34 relative digits.
        c=complex(big"1.25",big".2")
        j=jacobian_rmn(1,2,c,[1+big"1e-12"];precision=:quad,kind=2)
        @test only(j.dvalue_dcreal)≈complex(big"5.084528823480518360580906417683997069e6",big"-3.704671518492236309400812412367944905e6") rtol=big"1e-28"
        @test only(j.dderivative_dcreal)≈complex(big"-2.542264412107148688047455965219589695e18",big"1.852335759521177214548397423902695205e18") rtol=big"1e-28"
        j=jacobian_rmn(1,2,c,[big"0"];precision=:quad,kind=2,spheroid=:oblate)
        @test only(j.dvalue_dcreal)≈complex(big"14.46503453901280437825520071089769767",big"-10.50826690336147568389830712933445971") rtol=big"1e-28"
        @test only(j.dderivative_dcreal)≈complex(big"-33.24372188880857914258430590653508225",big"24.44921998275252128904198999440996139") rtol=big"1e-28"
    end
end
