using SpheroidalWaves, Test

# Inline references from an independent 384-bit Legendre/Bessel expansion,
# differentiated with a five-point stencil. Refining 80 to 104 terms and
# h=1e-10 to 1e-12 changed every entry by less than 6e-39.
# Columns: lambda_c, S_c, S_xc, R1_c, R1_xc, R2_c, R2_xc.
const PARAMETER_DERIVATIVE_REFERENCES = [
    (:prolate,false, [
        complex(big"1.04134447575229814513447473541903867714463073707905555004405055453489178829171256160062531382864018350827428130329657",big"0.0"),
        complex(big"-0.0513121504515224585547864221694417917086229280508057027823224811374064888107610905943325949778400243738225354799840598",big"0.0"),
        complex(big"-0.0575658726725869537453051243420078914757267442108460241358803927302013767511080852806641535064428438274206372840179939",big"0.0"),
        complex(big"0.225202829586489086277711787707580587517145972444825371960898112643335747142790172008508019785219671357233778208491892",big"0.0"),
        complex(big"-0.0574677229022064974251895607226830070295496940382235334109509760074957737728000579377542721166413156684474449533536187",big"0.0"),
        complex(big"0.965053128286810610792950316720844271850598093237587568971049974021578540542189397819009981142133922674903951516199622",big"0.0"),
        complex(big"-0.7727583239143221462076209645045089919037824652722906050373778862105538244098701402935949841548020145434250599821084",big"0.0")
    ]),
    (:prolate,true, [
        complex(big"1.04359943711343590382026656914914633285292878432125839614374439401105335048409863991369895958598356366369041401414584",big"0.157196763124144530054729220357014166616443689611064350846565580089742736842297101109988609225967740057036881815855036"),
        complex(big"-0.051374059203743321904356214048130137981277741390027742699126361489045036234095830833906897939066272834228764719748657",big"-0.00799524485846281266095524472196471787825856593286405743766086518820353032723136677237558993734113419635123580492066403"),
        complex(big"-0.0580693451944157968154497623794223500835953416656005339320663897126113428588502412870488081130962122309442099063050749",big"-0.00717462601977087477138085669263023122665140255591845868801550894655210011112597929037657979857066289284343103694373339"),
        complex(big"0.23900280214911622726143094214644909571552630388645093722060846658162529876080650465204987303823627398382163838751721",big"-0.0520278158634895885114689462880911001412441947663646847428004533376577164336355677294141634905025031203767460046104753"),
        complex(big"-0.0429876290638710047490999532799681492509420548677983760394304544660943487804155595632931624673739365869930556041303746",big"-0.165011374802337686049997603692306233229708143457092189571295904563024646288089598176738290139298655035294447632325405"),
        complex(big"0.819962497886195281842152285492907315439296972071253781271346991850044070122819731517392945348271626120293125689743724",big"-0.30857960248008913823369185877541675612365616595439106831375116236939512416881469021841607680728518524015380210764287"),
        complex(big"-0.4812748249017850223258836225434896437469022544238978382437244003765408669049958344396349179048128554215000863305096",big"0.67749915645313249278248258872414407676427105582308173957659665609005128915372944953358814746635425248413523558215406")
    ]),
    (:oblate,false, [
        complex(big"-1.1020397556358890209529722275273234534563466579988070568933440598884194385224558813396771231459084873276106807760664",big"0.0"),
        complex(big"0.052301839645565119154344545437399061626321929092102364486384750597748237215699377040372020396493832702225498789690102",big"0.0"),
        complex(big"0.0700370318259530140216816676338459456231704747351403922131318318817051717888280995157390131816019453018906605229173932",big"0.0"),
        complex(big"0.182379140037948545780052563704044405297010874601869896919234972000163855440067799999524728463386786621079098517201509",big"0.0"),
        complex(big"-0.182442797866573472205938310741076874946741569944314980041948015114967834064847494334218297915573502243755645479860059",big"0.0"),
        complex(big"0.811576586279751348598399089244799771495851663469075871319187066638468408910525701962764747179058511137577617348336963",big"0.0"),
        complex(big"-0.197141629400369066928382650697059973589438758128878068892882433456093995499288841439810781663270333421355829723326374",big"0.0")
    ]),
    (:oblate,true, [
        complex(big"-1.09965225173240791668049824207011897460348401928725911686665322611707288072499856485639272603924280741133972303952984",big"-0.18606011922462646957021183952753546041287161510801812134856661709410493418489202874356579184477289755279328483481348"),
        complex(big"0.0522902619603366409825151556852014227558615798783338955845832903085602504939332158216181642896617190795626440494683444",big"0.00846333518879194908077391262275843538269047309886958997467128701671686730953622773448194616987483177628353870547571647"),
        complex(big"0.0695954602453029766953491313183064703234467015526888620144377381412059824414424332750657187600860979842433543036895731",big"0.0130950312549669321758912970709879375947817581484013455034884105139895650322403759517704416930246655750343697449082054"),
        complex(big"0.198083825356482368551115851954835874302798360395832758704604564748739373400730693959931790449557974628645509842033137",big"-0.0991192915483226248046917486379512451694178361045442680411116692257512997021991480918368377639434900960370418870557795"),
        complex(big"-0.180392665365300671786774644288565223820173370958903221191142012253615217372571293912159784683957615236786998003517322",big"-0.17957197901357865112801002911627628655878540444926039874312582578385572571092011801851181910611653806634564096707391"),
        complex(big"0.731486050642069795088602125676088066633082983378585866912168652889547200093225798632982893738534158351533196143848252",big"-0.174011317689039288211180312454906176729583992252123102227005584061484689126251442047706565855473778924414748034863727"),
        complex(big"-0.036775120702678747829352452010564976992972020323148230069992566459775728847203267181574972248181846703760915738721502",big"0.268326045993420373017460298023387114464461064401848927353798883795352876603021042312179887962900970648765282352897022")
    ])
]

@testset "Parameter derivatives and reusable expansions" begin
    for precision in (:double,:quad)
        lib = SpheroidalWaves.backend_library(;precision)
        (lib === nothing || !isfile(lib)) && continue
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-28" : big"2e-12"
        radial_tolerance = precision === :quad ? big"1e-28" : big"5e-12"
        for (spheroid,complex_parameter,reference) in PARAMETER_DERIVATIVE_REFERENCES
            c = complex_parameter ? complex(T(5)/4,T(1)/5) : T(5)/4
            e = jacobian_eigen(1,2,c;precision,spheroid,with_metadata=true)
            @test (complex_parameter ? e.d_dcreal : e.derivative) ≈ reference[1] rtol=tolerance
            metadata = complex_parameter ? e.metadata_dcreal : e.metadata
            @test metadata.method === :coefficients
            @test metadata.step_used === nothing
            s = jacobian_smn(1,2,c,T[T(3)/10];precision,spheroid)
            sv = complex_parameter ? s.dvalue_dcreal : s.dvalue_dc
            sd = complex_parameter ? s.dderivative_dcreal : s.dderivative_dc
            @test only(sv) ≈ reference[2] rtol=tolerance
            @test only(sd) ≈ reference[3] rtol=tolerance
            unit = jacobian_smn(1,2,c,T[T(3)/10];precision,spheroid,normalize=true)
            @test (complex_parameter ? unit.dvalue_dcreal : unit.dvalue_dc)*sqrt(T(12)/5) ≈ sv rtol=tolerance
            if complex_parameter
                @test e.d_dcimag ≈ im*reference[1] rtol=tolerance
                @test only(s.dvalue_dcimag) ≈ im*reference[2] rtol=tolerance
            end
            for kind in (1,2)
                r = jacobian_rmn(1,2,c,T[2];precision,spheroid,kind,with_metadata=true)
                rv = complex_parameter ? r.dvalue_dcreal : r.dvalue_dc
                rd = complex_parameter ? r.dderivative_dcreal : r.dderivative_dc
                @test only(rv) ≈ reference[2kind+2] rtol=radial_tolerance
                @test only(rd) ≈ reference[2kind+3] rtol=radial_tolerance
                rm = complex_parameter ? r.metadata_value_dcreal : r.metadata_value
                @test rm.method === :differentiated_expansion
                @test rm.step_used === nothing
                if complex_parameter
                    @test only(r.dvalue_dcimag) ≈ im*reference[2kind+2] rtol=radial_tolerance
                    @test only(r.dderivative_dcimag) ≈ im*reference[2kind+3] rtol=radial_tolerance
                end
            end
        end

        # Perturbation of P_0: lambda_c/c -> 2sigma/3 and
        # S_c/c -> -2sigma*P_2/9. This remains informative when finite
        # differences of the function round to zero in the requested precision.
        for spheroid in (:prolate,:oblate), c in (T(0),T(1e-20),complex(T(1e-20),T(2e-20)))
            sigma = spheroid === :prolate ? 1 : -1
            e = jacobian_eigen(0,0,c;precision,spheroid)
            s = jacobian_smn(0,0,c,T[T(3)/10];precision,spheroid)
            v = c isa Real ? e : e.d_dcreal
            sv = only(c isa Real ? s.dvalue_dc : s.dvalue_dcreal)
            sd = only(c isa Real ? s.dderivative_dc : s.dderivative_dcreal)
            if iszero(c)
                @test iszero(v) && iszero(sv) && iszero(sd)
            else
                # Keep the O(c^2) value and coordinate slope as well as the
                # O(c) sensitivities; treating a tiny c as zero loses these.
                lambda = eigenvalue(0,0,c;precision,spheroid)
                wave = smn(0,0,c,T[T(3)/10];precision,spheroid)
                @test lambda/c^2 ≈ sigma*one(T)/3 rtol=tolerance
                @test only(wave.derivative)/c^2 ≈ -sigma*(T(3)/10)/3 rtol=tolerance
                @test v/c ≈ 2sigma*one(T)/3 rtol=tolerance
                @test sv/c ≈ sigma*(1-3*(T(3)/10)^2)/9 rtol=tolerance
                @test sd/c ≈ -2sigma*(T(3)/10)/3 rtol=tolerance
            end
        end

        # Reconstruct using spherical-limit functions, and ensure callers cannot
        # corrupt a cached expansion by modifying returned coefficient arrays.
        for spheroid in (:prolate,:oblate)
            c = T(5)/4
            expansion = SpheroidalWaves._angular_coefficients(1,2,c;precision,spheroid)
            @test expansion.converged
            @test expansion.quadrature_points == 0
            x = T[T(3)/10,T(7)/10]
            reconstructed = sum(d.*smn(1,l,zero(T),x;precision).value for (d,l) in zip(expansion.coefficients,expansion.degrees))
            @test reconstructed ≈ smn(1,2,c,x;precision,spheroid).value rtol=tolerance
            reconstructed_dc = sum(d.*smn(1,l,zero(T),x;precision).value for (d,l) in zip(expansion.dcoefficients_dc,expansion.degrees))
            @test reconstructed_dc ≈ jacobian_smn(1,2,c,x;precision,spheroid).dvalue_dc rtol=tolerance
            original = copy(expansion.coefficients)
            fill!(expansion.coefficients,zero(T))
            @test SpheroidalWaves._angular_coefficients(1,2,c;precision,spheroid).coefficients == original
        end
        small_c = precision === :quad ? big"1e-12" : 1e-6
        r = rmn(0,1,small_c,T[2];precision,kind=2)
        dr = jacobian_rmn(0,1,small_c,T[2];precision,kind=2,with_metadata=true)
        # The leading second-kind n=1 term is proportional to c^-2.
        @test only(small_c.*dr.dvalue_dc./r.value) ≈ -2 rtol=(precision === :quad ? big"1e-20" : 1e-9)
        @test dr.metadata_value.step_used === nothing
        endpoint = jacobian_smn(1,1,T(1.25),T[-1,1];precision,with_metadata=true)
        @test all(iszero,endpoint.dvalue_dc)
        @test all(isinf,endpoint.dderivative_dc)
        @test !endpoint.metadata_derivative.finite_flag
    end
    @test_throws ErrorException jacobian_eigen(0,0,1.0;h=Inf)
    @test_throws ArgumentError SpheroidalWaves._angular_coefficients(0,0,1.0;rtol=NaN)
end

@testset "Near-zero angular modes retain phase and scaling" begin
    for precision in (:double,:quad), spheroid in (:prolate,:oblate)
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-28" : 2e-13
        c,x = complex(T(1e-20),T(2e-20)),T(3)/10
        for (m,n) in ((1,1),(1,3),(2,3)), normalize in (false,true)
            spherical = smn(m,n,zero(c),x;precision,spheroid,normalize)
            wave = smn(m,n,c,x;precision,spheroid,normalize)
            @test wave.value ≈ spherical.value rtol=tolerance
            @test wave.derivative ≈ spherical.derivative rtol=tolerance
        end
        sigma = spheroid === :prolate ? 1 : -1
        scaled = smn(0,0,c,x;precision,spheroid,scaled=true,logderivative=true)
        slope = only(scaled.derivative.mantissa.*BigFloat(10).^scaled.derivative.exponent)
        @test slope/c^2 ≈ -sigma*x/3 rtol=tolerance
        @test only(scaled.logderivative)/c^2 ≈ -sigma*x/3 rtol=tolerance
        for parameter in (real(c),c)
            batch = smn(0,0:2,parameter,x;precision,spheroid)
            @test batch.value[:,1] ≈ smn(0,0,parameter,x;precision,spheroid).value rtol=tolerance
            @test batch.derivative[1,1]/parameter^2 ≈ -sigma*x/3 rtol=tolerance
        end
        @test accuracy(0,0,c,[x];precision,spheroid,target=:angular) == [-1]
    end
end
