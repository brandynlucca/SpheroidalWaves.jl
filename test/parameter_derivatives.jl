using SpheroidalWaves, Test

@testset "Jacobian diagnostics keyword" begin
    for c in (1.25, 1.25+0.125im),
        options in ((;), (; h = 1e-4, adaptive = false), (;
            order = 2, h = 1e-3, adaptive = false))

        for (jac, args) in ((jacobian_eigen, (0, 0, c)),
            (jacobian_smn, (0, 0, c, [0.25])),
            (jacobian_rmn, (0, 0, c, [2.0])),
            (jacobian_smn, (0, 0, c, [0.25+0.125im])),
            (jacobian_rmn, (0, 0, c, [2+0.125im])))
            ordinary = jac(args...; options...)
            reported = jac(args...; options..., diagnostics = true)
            @test isequal(ordinary, jac(args...; options..., diagnostics = false))
            @test any(name -> startswith(string(name), "metadata"), keys(reported))
            @test_throws MethodError jac(args...; options..., with_metadata = true)
            @test_throws MethodError jac(args...; options..., diagnostics = true, with_metadata = false)
        end
    end
end

@testset "Static radial sensitivities" begin
    offsets = collect(-4:4)
    first_weights, second_weights = SpheroidalWaves._difference_weights(offsets)
    for precision in (:double, :quad), spheroid in (:prolate, :oblate), kind in 1:2
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"2e-20" : 2e-9
        for c in (zero(T), complex(T(1.25), T(0.125)))
            h = big"0.001"
            wide = c isa Real ? BigFloat(c) : Complex{BigFloat}(c)
            samples = [rmn(1, 2, wide+k*h, [big"2"]; spheroid,
                           precision = :quad, kind, normalization = :static)
                       for k in offsets]
            first = jacobian_rmn(1, 2, c, [T(2)]; spheroid, precision, kind,
                normalization = :static, diagnostics = true)
            second = jacobian_rmn(1, 2, c, [T(2)]; spheroid, precision, kind,
                normalization = :static, order = 2, diagnostics = true)
            for (field, a, b) in ((:value, :dvalue_dc, :d2value_dc2),
                (:derivative, :dderivative_dc, :d2derivative_dc2))
                first_value = c isa Real ? getproperty(first, a) :
                              getproperty(first, field===:value ? :dvalue_dcreal :
                                                 :dderivative_dcreal)
                second_value = c isa Real ? getproperty(second, b) :
                               getproperty(second, field===:value ? :d2value_dcreal2 :
                                                   :d2derivative_dcreal2)
                values = [only(getproperty(r, field)) for r in samples]
                @test only(first_value)≈sum(BigFloat.(first_weights) .* values)/h rtol=tolerance atol=tolerance
                @test only(second_value)≈sum(BigFloat.(second_weights) .*
                                             (values .- values[5]))/h^2 rtol=tolerance atol=tolerance
            end
            @test second.metadata.finite_flag
            fd = jacobian_rmn(1, 2, c, [T(2)]; spheroid, precision,
                kind, normalization = :static, h = 1e-4)
            key = c isa Real ? :dvalue_dc : :dvalue_dcreal
            @test getproperty(fd, key)≈getproperty(first, key) rtol=1e-6 atol=1e-9
        end
        tiny = jacobian_rmn(1, 2, big"1e-40", [2]; spheroid,
            precision = :quad, kind, normalization = :static)
        limit = jacobian_rmn(1, 2, big"0", [2]; spheroid, precision = :quad,
            kind, normalization = :static, order = 2)
        @test tiny.dvalue_dc ./ big"1e-40" ≈ limit.d2value_dc2 rtol=big"1e-22"
        complex_zero = jacobian_rmn(1, 2, complex(zero(T)), [T(2)]; spheroid,
            precision, kind, normalization = :static, h = 1e-4)
        @test all(iszero, complex_zero.dvalue_dcreal)
    end
    for precision in (:double, :quad), c in (1.25, 1.25+0im, 1.25+0.125im)

        for (jac, x) in ((jacobian_smn, -1.0), (jacobian_smn, 1.0), (jacobian_rmn, 1.0))
            result = jac(1, 2, c, [x]; precision, order = 2, diagnostics = true)
            d = c isa Real ? result.d2derivative_dc2 : result.d2derivative_dcreal2
            @test isinf(only(d)) && !isnan(only(d))
            @test result.metadata.conditioning_flag === :singular
            if c isa Complex
                @test result.d2derivative_dcreal_dcimag == complex.(-imag.(d), real.(d))
                first = jac(1, 2, c, [x]; precision)
                @test !any(isnan, first.dderivative_dcimag)
            end
        end
        for cstatic in (zero(c), c)
            result = jacobian_rmn(
                1, 2, cstatic, [1]; precision, normalization = :static, order = 2)
            @test isinf(only(cstatic isa Real ? result.d2derivative_dc2 :
                             result.d2derivative_dcreal2))
            first = jacobian_rmn(1, 2, cstatic, [1]; precision, normalization = :static)
            value = cstatic isa Real ? first.dderivative_dc : first.dderivative_dcreal
            @test iszero(cstatic) ? all(iszero, value) : all(isinf, value)
        end
    end
    @test isnan(only(jacobian_rmn(0, 0, 0, [1]; normalization = :static, kind = 2).dvalue_dc))
end

@testset "Second parameter derivatives" begin
    # Nine function values supply an independent second-derivative stencil.
    offsets = collect(-4:4)
    weights = last(SpheroidalWaves._difference_weights(offsets))
    function reference(f, c, h)
        center = f(c)
        sum(BigFloat(w) .* (f(c+k*h) .- center) for (k, w) in zip(offsets, weights)) ./ h^2
    end
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T = precision === :quad ? BigFloat : Float64
        tol = precision === :quad ? big"5e-22" : 3e-10
        for c in (T(1.25), complex(T(1.25), T(0.125)))
            wide = c isa Real ? BigFloat(c) : Complex{BigFloat}(c)
            h = big"0.0001"
            e = jacobian_eigen(1, 2, c; spheroid, precision, order = 2, diagnostics = true)
            ev = c isa Real ? e.derivative : e.d2_dcreal2
            @test ev ≈
                  reference(z -> eigenvalue(1, 2, z; spheroid, precision = :quad), wide, h) rtol=tol
            @test e.metadata.method === :differenced_sensitivity && e.metadata.finite_flag
            @test e.metadata.step_used > 0
            @test ev isa (c isa Real ? T : Complex{T})
            @test jacobian_eigen(1, 2, c; spheroid, precision, order = 2) ==
                  (c isa Real ? ev :
                   (; d2_dcreal2 = ev, d2_dcreal_dcimag = im*ev, d2_dcimag2 = -ev))
            for (f, jac, point, kinds) in ((smn, jacobian_smn, T(0.3), (1, 2)),
                (rmn, jacobian_rmn, T(2), (1, 2, 3, 4)))
                for kind in kinds
                    j = jac(1, 2, c, [point]; spheroid, precision,
                        kind, order = 2, diagnostics = true)
                    v = c isa Real ? j.d2value_dc2 : j.d2value_dcreal2
                    d = c isa Real ? j.d2derivative_dc2 : j.d2derivative_dcreal2
                    if kind <= 2
                        values(z) = f(
                            1, 2, z, [BigFloat(point)]; spheroid, precision = :quad, kind)
                        @test v ≈ reference(z -> values(z).value, wide, h) rtol=tol
                        @test d ≈ reference(z -> values(z).derivative, wide, h) rtol=tol
                    else
                        a = jac(1, 2, c, [point]; spheroid, precision, kind = 1, order = 2)
                        b = jac(1, 2, c, [point]; spheroid, precision, kind = 2, order = 2)
                        av, bv = c isa Real ? (a.d2value_dc2, b.d2value_dc2) :
                                 (a.d2value_dcreal2, b.d2value_dcreal2)
                        @test v ≈ av+(kind==3 ? im : -im)*bv rtol=tol
                    end
                    @test j.metadata.finite_flag &&
                          j.metadata.method === :differenced_sensitivity
                    @test eltype(v) == (f===rmn || c isa Complex ? Complex{T} : T)
                    if c isa Complex
                        @test j.d2value_dcreal_dcimag == im .* v
                        @test j.d2value_dcimag2 == -v
                        @test j.d2derivative_dcreal_dcimag == im .* d
                        @test j.d2derivative_dcimag2 == -d
                    end
                end
            end
        end
        sigma = spheroid === :prolate ? 1 : -1
        for c in (zero(T), T(1e-20))
            @test jacobian_eigen(0, 0, c; spheroid, precision, order = 2) ≈ T(2sigma)/3 rtol=tol
            j = jacobian_smn(0, 0, c, T[0, 0.5]; spheroid, precision, order = 2)
            @test j.d2value_dc2 ≈ T(sigma) .* (1 .- 3 .* T[0, 0.5] .^ 2) ./ 9 rtol=tol
            @test j.d2derivative_dc2 ≈ -T(2sigma)/3 .* T[0, 0.5] rtol=tol
        end
    end
    @test_throws ArgumentError jacobian_eigen(0, 0, 1; order = 3)
    @test_throws ArgumentError jacobian_smn(0, 0, 1, [0.2]; order = 0)
    @test_throws ArgumentError jacobian_rmn(0, 0, 1, [2]; order = -1)
    @test_throws ArgumentError jacobian_rmn(0, 0, 1, [2]; order = 2, h = 0.5)
    @test_throws ArgumentError jacobian_eigen(0, 0, 1; order = 2, h = eps()/10)
    @test_throws r"h must be finite" jacobian_eigen(0, 0, 1; order = 2, h = -1)
    explicit = jacobian_eigen(
        0, 0, 1.0; order = 2, h = 1e-3, adaptive = false, diagnostics = true)
    @test explicit.metadata.step_used == 1e-3
    @test explicit.derivative ≈ jacobian_eigen(0, 0, 1.0; order = 2) rtol=1e-9
    bad = jacobian_smn(0, 0, 1.0, [-1.0, 1.0]; kind = 2, order = 2, diagnostics = true)
    @test all(isnan, bad.d2value_dc2) && !bad.metadata.finite_flag
    @test bad.metadata.conditioning_flag === :poor
    small = jacobian_rmn(0, 0, 1e-8, [2.0]; order = 2, diagnostics = true)
    @test small.metadata.step_used < 1e-8/4
    @test all(isfinite, small.d2value_dc2)
    for precision in (:double, :quad)
        unresolved = jacobian_rmn(
            0, 1, big"1e-50", [2.0]; precision, order = 2, diagnostics = true)
        @test unresolved.metadata.conditioning_flag === :poor
        @test unresolved.metadata.suggested_action ===
              (precision===:double ? :use_quad : :unresolved)
    end
    imaginary = jacobian_smn(1, 2, 1.25im, [0.3]; order = 2)
    opposite = jacobian_smn(1, 2, 1.25, [0.3]; spheroid = :oblate, order = 2)
    @test imaginary.d2value_dcreal2 ≈ -opposite.d2value_dc2 rtol=1e-9
    # Independent mixed real/imaginary differences of function values.
    c, h = complex(big"1.25", big"0.125"), big"0.0001"
    offsets, weights = (-2, -1, 1, 2), (1, -8, 8, -1)
    mixed = sum(a*b*only(smn(1, 2, c+(i+im*j)*h, [big"0.3"]; precision = :quad).value)
    for (i, a) in zip(offsets, weights), (j, b) in zip(offsets, weights))/(144h^2)
    @test only(jacobian_smn(1, 2, c, [big"0.3"]; precision = :quad, order = 2).d2value_dcreal_dcimag) ≈
          mixed rtol=1e-14
end

@testset "Second integral eigenvalue derivatives" begin
    for precision in (:double, :quad)
        T = precision === :quad ? BigFloat : Float64
        tol = precision === :quad ? big"1e-24" : 2e-9
        for n in 0:3
            @test jacobian_eigen(
                0, n, zero(T); operator = :concentration, precision, order = 2) == 0
            @test jacobian_eigen(0, n, zero(T); operator = :concentration,
                form = :complement, precision, order = 2) == 0
            @test jacobian_eigen(
                0, n, zero(T); operator = :fourier, precision, order = 2) ==
                  (n==0 ? -T(2)/9 : n==2 ? -T(8)/45 : zero(T))
            singular = jacobian_eigen(
                0, n, zero(T); operator = :concentration, form = :log,
                precision, order = 2, diagnostics = true)
            @test singular.derivative == -Inf &&
                  singular.metadata.conditioning_flag === :singular
            for operator in (:concentration, :fourier)
                c, h = big"1.25", big"0.0001"
                offsets = collect(-4:4)
                weights = last(SpheroidalWaves._difference_weights(offsets))
                f(z) = eigenvalue(0, n, z; operator, precision = :quad)
                expected = sum(BigFloat(w)*(f(c+k*h)-f(c))
                for (k, w) in zip(offsets, weights))/h^2
                a = jacobian_eigen(0, n, T(c); operator, precision, order = 2)
                @test a ≈ expected rtol=tol
                if operator === :concentration
                    @test jacobian_eigen(
                        0, n, T(c); operator, precision, form = :complement, order = 2) ≈ -a rtol=tol
                    first = jacobian_eigen(0, n, T(c); operator, precision)
                    value = eigenvalue(0, n, T(c); operator, precision)
                    @test jacobian_eigen(
                        0, n, T(c); operator, precision, form = :log, order = 2) ≈
                          a/value-(first/value)^2 rtol=tol
                end
            end
        end
    end
    @test_throws DomainError jacobian_eigen(
        0, 0, 0; operator = :concentration, form = :log, order = 2, h = 0.01)
    @test_throws ArgumentError jacobian_eigen(
        0, 0, 1; operator = :concentration, order = 2, h = 1)
    forward = jacobian_eigen(
        0, 0, 0.0; operator = :fourier, order = 2, h = 0.001, diagnostics = true)
    @test forward.derivative ≈ -2/9 rtol=1e-9
    @test forward.metadata.method === :differenced_sensitivity
    tiny = jacobian_eigen(0, 2, big"1e-40"; operator = :concentration,
        form = :log, precision = :quad, order = 2)
    @test tiny ≈ -5/big"1e-80" rtol=big"1e-25"
end

@testset "Explicit finite differences agree with analytic sensitivities" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T = precision === :quad ? BigFloat : Float64
        h = precision === :quad ? T(10)^(-10) : T(1e-4)
        tolerance = precision === :quad ? T(1e-17) : T(2e-7)
        for c in (T(1.25), complex(T(1.25), T(0.2))), metadata in (false, true)

            exact = jacobian_eigen(1, 2, c; precision, spheroid)
            fd = jacobian_eigen(1, 2, c; precision, spheroid, h,
                diagnostics = metadata, adaptive = metadata)
            if c isa Real
                @test (metadata ? fd.derivative : fd) ≈ exact rtol=tolerance
                if metadata
                    @test fd.metadata.method === :finite_difference
                    @test fd.metadata.suggested_action === :accept
                end
            else
                @test fd.d_dcreal ≈ exact.d_dcreal rtol=tolerance
                @test fd.d_dcimag ≈ im*exact.d_dcreal rtol=tolerance
                metadata && @test fd.metadata_dcimag.finite_flag
            end
            for (differentiate, points) in ((jacobian_smn, T[0.3]), (jacobian_rmn, T[2]))
                exact = differentiate(1, 2, c, points; precision, spheroid)
                fd = differentiate(1, 2, c, points; precision, spheroid, h,
                    diagnostics = metadata, adaptive = metadata)
                for field in keys(exact)
                    @test getproperty(fd, field) ≈ getproperty(exact, field) rtol=tolerance
                end
                if metadata
                    fields = c isa Real ? (:metadata_value, :metadata_derivative) :
                             (:metadata_value_dcreal, :metadata_value_dcimag,
                        :metadata_derivative_dcreal, :metadata_derivative_dcimag)
                    for field in fields
                        @test getproperty(fd, field).method === :finite_difference
                        @test getproperty(fd, field).finite_flag
                    end
                end
            end
        end
    end
end

@testset "Finite difference reliability and refinement" begin
    SW = SpheroidalWaves
    for precision in (:double, :quad)
        for (coarse, fine, flag, action) in (
            (1.0, 1.0, :good, :accept),
            (1.0001, 1.0, :warning, precision === :quad ? :accept : :retry_smaller_h),
            (2.0, 1.0, :poor, precision === :quad ? :retry_smaller_h : :use_quad),
            (Inf, 1.0, :poor, :retry_smaller_h))
            md = SW._jacobian_metadata(
                coarse, fine, 0.1; precision, rtol = 1e-6, atol = 1e-10)
            @test md.conditioning_flag === flag
            @test md.suggested_action === action
            @test md.finite_flag == isfinite(coarse)
        end
        # A cubic has a centered derivative error exactly h^2. Check which
        # stencil is returned and that poor consistency triggers refinement.
        steps = Float64[]
        calc(h) = (push!(steps, h); ((1+h)^3-(1-h)^3)/(2h))
        derivative, md = SW._finite_difference_diagnostics(calc, 0.5;
            precision, adaptive = true, rtol = 1e-6, atol = 1e-10)
        @test steps == [0.5, 0.25, 0.125]
        @test derivative == 3 + 0.125^2
        @test md.step_used == 0.125
        empty!(steps)
        derivative, md = SW._finite_difference_diagnostics(calc, 0.5;
            precision, adaptive = false, rtol = 1e-6, atol = 1e-10)
        @test steps == [0.5, 0.25]
        @test derivative == 3.25
        @test md.step_used == 0.5
    end
end

# Inline references from an independent 384-bit Legendre/Bessel expansion,
# differentiated with a five-point stencil. Refining 80 to 104 terms and
# h=1e-10 to 1e-12 changed every entry by less than 6e-39.
# Columns: lambda_c, S_c, S_xc, R1_c, R1_xc, R2_c, R2_xc.
const PARAMETER_DERIVATIVE_REFERENCES = [
    (:prolate,
        false,
        [
            complex(
                big"1.04134447575229814513447473541903867714463073707905555004405055453489178829171256160062531382864018350827428130329657",
                big"0.0"),
            complex(
                big"-0.0513121504515224585547864221694417917086229280508057027823224811374064888107610905943325949778400243738225354799840598",
                big"0.0"),
            complex(
                big"-0.0575658726725869537453051243420078914757267442108460241358803927302013767511080852806641535064428438274206372840179939",
                big"0.0"),
            complex(
                big"0.225202829586489086277711787707580587517145972444825371960898112643335747142790172008508019785219671357233778208491892",
                big"0.0"),
            complex(
                big"-0.0574677229022064974251895607226830070295496940382235334109509760074957737728000579377542721166413156684474449533536187",
                big"0.0"),
            complex(
                big"0.965053128286810610792950316720844271850598093237587568971049974021578540542189397819009981142133922674903951516199622",
                big"0.0"),
            complex(
                big"-0.7727583239143221462076209645045089919037824652722906050373778862105538244098701402935949841548020145434250599821084",
                big"0.0")
        ]),
    (:prolate,
        true,
        [
            complex(
                big"1.04359943711343590382026656914914633285292878432125839614374439401105335048409863991369895958598356366369041401414584",
                big"0.157196763124144530054729220357014166616443689611064350846565580089742736842297101109988609225967740057036881815855036"),
            complex(
                big"-0.051374059203743321904356214048130137981277741390027742699126361489045036234095830833906897939066272834228764719748657",
                big"-0.00799524485846281266095524472196471787825856593286405743766086518820353032723136677237558993734113419635123580492066403"),
            complex(
                big"-0.0580693451944157968154497623794223500835953416656005339320663897126113428588502412870488081130962122309442099063050749",
                big"-0.00717462601977087477138085669263023122665140255591845868801550894655210011112597929037657979857066289284343103694373339"),
            complex(
                big"0.23900280214911622726143094214644909571552630388645093722060846658162529876080650465204987303823627398382163838751721",
                big"-0.0520278158634895885114689462880911001412441947663646847428004533376577164336355677294141634905025031203767460046104753"),
            complex(
                big"-0.0429876290638710047490999532799681492509420548677983760394304544660943487804155595632931624673739365869930556041303746",
                big"-0.165011374802337686049997603692306233229708143457092189571295904563024646288089598176738290139298655035294447632325405"),
            complex(
                big"0.819962497886195281842152285492907315439296972071253781271346991850044070122819731517392945348271626120293125689743724",
                big"-0.30857960248008913823369185877541675612365616595439106831375116236939512416881469021841607680728518524015380210764287"),
            complex(
                big"-0.4812748249017850223258836225434896437469022544238978382437244003765408669049958344396349179048128554215000863305096",
                big"0.67749915645313249278248258872414407676427105582308173957659665609005128915372944953358814746635425248413523558215406")
        ]),
    (:oblate,
        false,
        [
            complex(
                big"-1.1020397556358890209529722275273234534563466579988070568933440598884194385224558813396771231459084873276106807760664",
                big"0.0"),
            complex(
                big"0.052301839645565119154344545437399061626321929092102364486384750597748237215699377040372020396493832702225498789690102",
                big"0.0"),
            complex(
                big"0.0700370318259530140216816676338459456231704747351403922131318318817051717888280995157390131816019453018906605229173932",
                big"0.0"),
            complex(
                big"0.182379140037948545780052563704044405297010874601869896919234972000163855440067799999524728463386786621079098517201509",
                big"0.0"),
            complex(
                big"-0.182442797866573472205938310741076874946741569944314980041948015114967834064847494334218297915573502243755645479860059",
                big"0.0"),
            complex(
                big"0.811576586279751348598399089244799771495851663469075871319187066638468408910525701962764747179058511137577617348336963",
                big"0.0"),
            complex(
                big"-0.197141629400369066928382650697059973589438758128878068892882433456093995499288841439810781663270333421355829723326374",
                big"0.0")
        ]),
    (:oblate,
        true,
        [
            complex(
                big"-1.09965225173240791668049824207011897460348401928725911686665322611707288072499856485639272603924280741133972303952984",
                big"-0.18606011922462646957021183952753546041287161510801812134856661709410493418489202874356579184477289755279328483481348"),
            complex(
                big"0.0522902619603366409825151556852014227558615798783338955845832903085602504939332158216181642896617190795626440494683444",
                big"0.00846333518879194908077391262275843538269047309886958997467128701671686730953622773448194616987483177628353870547571647"),
            complex(
                big"0.0695954602453029766953491313183064703234467015526888620144377381412059824414424332750657187600860979842433543036895731",
                big"0.0130950312549669321758912970709879375947817581484013455034884105139895650322403759517704416930246655750343697449082054"),
            complex(
                big"0.198083825356482368551115851954835874302798360395832758704604564748739373400730693959931790449557974628645509842033137",
                big"-0.0991192915483226248046917486379512451694178361045442680411116692257512997021991480918368377639434900960370418870557795"),
            complex(
                big"-0.180392665365300671786774644288565223820173370958903221191142012253615217372571293912159784683957615236786998003517322",
                big"-0.17957197901357865112801002911627628655878540444926039874312582578385572571092011801851181910611653806634564096707391"),
            complex(
                big"0.731486050642069795088602125676088066633082983378585866912168652889547200093225798632982893738534158351533196143848252",
                big"-0.174011317689039288211180312454906176729583992252123102227005584061484689126251442047706565855473778924414748034863727"),
            complex(
                big"-0.036775120702678747829352452010564976992972020323148230069992566459775728847203267181574972248181846703760915738721502",
                big"0.268326045993420373017460298023387114464461064401848927353798883795352876603021042312179887962900970648765282352897022")
        ])
]

@testset "Parameter derivatives and reusable expansions" begin
    for precision in (:double, :quad)
        lib = SpheroidalWaves.backend_library(; precision)
        (lib === nothing || !isfile(lib)) && continue
        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-28" : big"2e-12"
        radial_tolerance = precision === :quad ? big"1e-28" : big"5e-12"
        for (spheroid, complex_parameter, reference) in PARAMETER_DERIVATIVE_REFERENCES
            c = complex_parameter ? complex(T(5)/4, T(1)/5) : T(5)/4
            e = jacobian_eigen(1, 2, c; precision, spheroid, diagnostics = true)
            @test (complex_parameter ? e.d_dcreal : e.derivative) ≈ reference[1] rtol=tolerance
            metadata = complex_parameter ? e.metadata_dcreal : e.metadata
            @test metadata.method === :coefficients
            @test metadata.step_used === nothing
            s = jacobian_smn(1, 2, c, T[T(3) / 10]; precision, spheroid)
            sv = complex_parameter ? s.dvalue_dcreal : s.dvalue_dc
            sd = complex_parameter ? s.dderivative_dcreal : s.dderivative_dc
            @test only(sv) ≈ reference[2] rtol=tolerance
            @test only(sd) ≈ reference[3] rtol=tolerance
            unit = jacobian_smn(
                1, 2, c, T[T(3) / 10]; precision, spheroid, normalize = true)
            @test (complex_parameter ? unit.dvalue_dcreal : unit.dvalue_dc)*sqrt(T(12)/5) ≈
                  sv rtol=tolerance
            if complex_parameter
                @test e.d_dcimag ≈ im*reference[1] rtol=tolerance
                @test only(s.dvalue_dcimag) ≈ im*reference[2] rtol=tolerance
            end
            for kind in (1, 2)
                r = jacobian_rmn(
                    1, 2, c, T[2]; precision, spheroid, kind, diagnostics = true)
                rv = complex_parameter ? r.dvalue_dcreal : r.dvalue_dc
                rd = complex_parameter ? r.dderivative_dcreal : r.dderivative_dc
                @test only(rv) ≈ reference[2kind + 2] rtol=radial_tolerance
                @test only(rd) ≈ reference[2kind + 3] rtol=radial_tolerance
                rm = complex_parameter ? r.metadata_value_dcreal : r.metadata_value
                @test rm.method === :differentiated_expansion
                @test rm.step_used === nothing
                if complex_parameter
                    @test only(r.dvalue_dcimag) ≈ im*reference[2kind + 2] rtol=radial_tolerance
                    @test only(r.dderivative_dcimag) ≈ im*reference[2kind + 3] rtol=radial_tolerance
                end
            end
        end

        # Perturbation of P_0: lambda_c/c -> 2sigma/3 and
        # S_c/c -> -2sigma*P_2/9. This remains informative when finite
        # differences of the function round to zero in the requested precision.
        for spheroid in (:prolate, :oblate),
            c in (T(0), T(1e-20), complex(T(1e-20), T(2e-20)))

            sigma = spheroid === :prolate ? 1 : -1
            e = jacobian_eigen(0, 0, c; precision, spheroid)
            s = jacobian_smn(0, 0, c, T[T(3) / 10]; precision, spheroid)
            v = c isa Real ? e : e.d_dcreal
            sv = only(c isa Real ? s.dvalue_dc : s.dvalue_dcreal)
            sd = only(c isa Real ? s.dderivative_dc : s.dderivative_dcreal)
            if iszero(c)
                @test iszero(v) && iszero(sv) && iszero(sd)
            else
                # Keep the O(c^2) value and coordinate slope as well as the
                # O(c) sensitivities; treating a tiny c as zero loses these.
                lambda = eigenvalue(0, 0, c; precision, spheroid)
                wave = smn(0, 0, c, T[T(3) / 10]; precision, spheroid)
                @test lambda/c^2 ≈ sigma*one(T)/3 rtol=tolerance
                @test only(wave.derivative)/c^2 ≈ -sigma*(T(3)/10)/3 rtol=tolerance
                @test v/c ≈ 2sigma*one(T)/3 rtol=tolerance
                @test sv/c ≈ sigma*(1-3*(T(3)/10)^2)/9 rtol=tolerance
                @test sd/c ≈ -2sigma*(T(3)/10)/3 rtol=tolerance
            end
        end

        # Reconstruct using spherical-limit functions, and ensure callers cannot
        # corrupt a cached expansion by modifying returned coefficient arrays.
        for spheroid in (:prolate, :oblate)
            c = T(5)/4
            expansion = SpheroidalWaves._angular_coefficients(1, 2, c; precision, spheroid)
            @test expansion.converged
            @test expansion.quadrature_points == 0
            x = T[T(3) / 10, T(7) / 10]
            reconstructed = sum(d .* smn(1, l, zero(T), x; precision).value
            for (d, l) in zip(expansion.coefficients, expansion.degrees))
            @test reconstructed ≈ smn(1, 2, c, x; precision, spheroid).value rtol=tolerance
            reconstructed_dc = sum(d .* smn(1, l, zero(T), x; precision).value
            for (d, l) in zip(expansion.dcoefficients_dc, expansion.degrees))
            @test reconstructed_dc ≈ jacobian_smn(1, 2, c, x; precision, spheroid).dvalue_dc rtol=tolerance
            original = copy(expansion.coefficients)
            fill!(expansion.coefficients, zero(T))
            @test SpheroidalWaves._angular_coefficients(1, 2, c; precision, spheroid).coefficients ==
                  original
        end
        small_c = precision === :quad ? big"1e-12" : 1e-6
        r = rmn(0, 1, small_c, T[2]; precision, kind = 2)
        dr = jacobian_rmn(0, 1, small_c, T[2]; precision, kind = 2, diagnostics = true)
        # The leading second-kind n=1 term is proportional to c^-2.
        @test only(small_c .* dr.dvalue_dc ./ r.value) ≈ -2 rtol=(precision === :quad ?
                                                                  big"1e-20" : 1e-9)
        @test dr.metadata_value.step_used === nothing
        endpoint = jacobian_smn(1, 1, T(1.25), T[-1, 1]; precision, diagnostics = true)
        @test all(iszero, endpoint.dvalue_dc)
        @test all(isinf, endpoint.dderivative_dc)
        @test !endpoint.metadata_derivative.finite_flag
    end
    @test_throws ErrorException jacobian_eigen(0, 0, 1.0; h = Inf)
    @test_throws ArgumentError SpheroidalWaves._angular_coefficients(0, 0, 1.0; rtol = NaN)
end

@testset "Near-zero angular modes retain phase and scaling" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T = precision === :quad ? BigFloat : Float64
        tolerance = precision === :quad ? big"1e-28" : 2e-13
        c, x = complex(T(1e-20), T(2e-20)), T(3)/10
        for (m, n) in ((1, 1), (1, 3), (2, 3)), normalize in (false, true)

            spherical = smn(m, n, zero(c), x; precision, spheroid, normalize)
            wave = smn(m, n, c, x; precision, spheroid, normalize)
            @test wave.value ≈ spherical.value rtol=tolerance
            @test wave.derivative ≈ spherical.derivative rtol=tolerance
        end
        sigma = spheroid === :prolate ? 1 : -1
        scaled = smn(0, 0, c, x; precision, spheroid, scaled = true, logderivative = true)
        slope = only(scaled.derivative.mantissa .*
                     BigFloat(10) .^ scaled.derivative.exponent)
        @test slope/c^2 ≈ -sigma*x/3 rtol=tolerance
        @test only(scaled.logderivative)/c^2 ≈ -sigma*x/3 rtol=tolerance
        for parameter in (real(c), c)
            batch = smn(0, 0:2, parameter, x; precision, spheroid)
            @test batch.value[:, 1] ≈ smn(0, 0, parameter, x; precision, spheroid).value rtol=tolerance
            @test batch.derivative[1, 1]/parameter^2 ≈ -sigma*x/3 rtol=tolerance
        end
        @test accuracy(0, 0, c, [x]; precision, spheroid, target = :angular) == [-1]
    end
end
@testset "Complex coordinate parameter sensitivities" begin
    for (jac, z, options) in ((jacobian_smn, 1+0im, (;)),
        (jacobian_rmn, 1+0im, (;)),
        (jacobian_rmn, 1+0im, (; normalization = :static)),
        (jacobian_rmn, 0im, (; spheroid = :oblate)))
        c = big"1.25"+big"0.125"*im
        analytic = jac(0, 0, c, z; precision = :quad, options...)
        stencil = jac(0, 0, c, z; precision = :quad, h = big"1e-6", options...)
        @test analytic.dvalue_dcreal ≈ stencil.dvalue_dcreal rtol=big"1e-10"
        @test analytic.dderivative_dcreal≈stencil.dderivative_dcreal rtol=big"1e-10" atol=big"1e-40"
    end
    singular = jacobian_smn(0, 0, 1.25+0.125im, 1+0im; kind = 2, h = 1e-5)
    @test all(isnan, singular.dvalue_dcreal)
    for spheroid in (:prolate, :oblate),
        (f, jac, z) in ((smn, jacobian_smn, 0.3+0.2im), (rmn, jacobian_rmn, 1.6+0.3im)),
        kind in 1:2
        c = big"1.25"+big"0.125"*im
        result = jac(1, 2, c, z; spheroid, precision = :quad, kind, diagnostics = true)
        difference = jac(1, 2, c, z; spheroid, precision = :quad,
            kind, h = big"1e-6", diagnostics = true)
        @test result.dvalue_dcreal ≈ difference.dvalue_dcreal rtol=big"1e-10"
        @test result.dderivative_dcimag ≈ difference.dderivative_dcimag rtol=big"1e-10"
        @test result.dvalue_dcimag == im*result.dvalue_dcreal
        @test result.dderivative_dcimag == im*result.dderivative_dcreal
        @test result.metadata_value_dcreal.method === :differentiated_equation
        @test result.metadata_derivative_dcimag.conditioning_flag === :good
        second = jac(
            1, 2, c, z; spheroid, precision = :quad, kind, order = 2, diagnostics = true)
        h = big"1e-6"
        jp = jac(1, 2, c+h, z; spheroid, precision = :quad, kind)
        jm = jac(1, 2, c-h, z; spheroid, precision = :quad, kind)
        @test second.d2value_dcreal2 ≈ (jp.dvalue_dcreal-jm.dvalue_dcreal)/(2h) rtol=big"1e-10"
        @test second.d2value_dcimag2 == -second.d2value_dcreal2
        @test second.d2derivative_dcreal_dcimag == im*second.d2derivative_dcreal2
        real_result = jac(
            1, 2, real(c), z; spheroid, precision = :quad, kind, diagnostics = true)
        real_difference = jac(1, 2, real(c), z; spheroid, precision = :quad,
            kind, h, adaptive = false, diagnostics = true)
        @test real_result.dvalue_dc ≈ real_difference.dvalue_dc rtol=big"1e-10"
        @test real_result.dderivative_dc ≈ real_difference.dderivative_dc rtol=big"1e-10"
        # Differentiate the returned third derivative along an imaginary step.
        high = f(1, 2, c, z; spheroid, precision = :quad, kind, derivatives = 4)
        plus = f(1, 2, c, z+im*h; spheroid, precision = :quad, kind, derivatives = 3)
        minus = f(1, 2, c, z-im*h; spheroid, precision = :quad, kind, derivatives = 3)
        @test high.fourth_derivative ≈
              (plus.third_derivative-minus.third_derivative)/(2im*h) rtol=big"1e-9"
    end
    for spheroid in (:prolate, :oblate), kind in 1:2

        z, c, h = big"1.6"+big"0.3"*im, big"1.25", big"1e-6"
        result = jacobian_rmn(
            1, 2, c, z; spheroid, precision = :quad, kind, normalization = :static)
        plus = rmn(1, 2, c+h, z; spheroid, precision = :quad, kind, normalization = :static)
        minus = rmn(
            1, 2, c-h, z; spheroid, precision = :quad, kind, normalization = :static)
        @test result.dvalue_dc ≈ (plus.value-minus.value)/(2h) rtol=big"1e-10"
        @test all(iszero, jacobian_rmn(1, 2, 0, z; spheroid, kind, normalization = :static).dvalue_dc)
    end
    for f in (jacobian_smn, jacobian_rmn)
        @test_throws ArgumentError f(0, 0, 1, 0.3im; order = 3)
        @test_throws r"h must be finite and positive" f(0, 0, 1, 0.3im; h = 0)
        @test_throws r"rtol must be positive" f(0, 0, 1, 0.3im; rtol = -1)
    end
    @test_throws ArgumentError jacobian_smn(0, 0, 1, 0.3im; kind = 2, normalize = true)
end
