using SpheroidalWaves, Test

# Independent mpmath 1.3.0 references: scripts/generate_accuracy_references.py::legendre_mode, (dps, terms) = (110, 24), (140, 24), (140, 32).
const ARBITRARY_PRECISION_REFERENCES = (
    (spheroid = :prolate, c = 1//1, values = (
        ("0.319000055146892739783981985871826494103814376510741729114970059238053009763670379743178473870059809926688", "0.0"),
        ("0.998846198637306956174569062173684719853723656448361917202611635052765819775381160947447421157795223328458", "0.0"),
        ("-0.107375553229029149344775406440661507612183363536558525875293275723623535320412663687098871361155237685626", "0.0"),
        ("0.00182376400585983740216381880282305164314377100076956255832363512543755528097767133789908204621056456577148", "0.0"),
    )),
    (spheroid = :prolate, c = complex(1//1, 1//4), values = (
        ("0.303235577304105060785313344346188701012396500405900347760112187638888468163625765174587690518159352874124", "0.153352965506012687954704507348706414785976931529398947646079078147847292101453723983951842226194416353958"),
        ("0.999232288840524133096463116801833466442215717279772280521354160394495958507126487605553604601738495525586", "-0.00106084747003089812108916450597025832638806795320135074362656796214833690281975743416525958313126030121808"),
        ("-0.101904800169774246090681118395964930748183786113276206784680852033022815106964828417757753851495881969989", "-0.0519999887973533157285592857635032433861308591701146716982784477297271512894479656215853408868364628381157"),
        ("0.00119020250787946888063174683955329808638039614136031395971608309343656918293893382796707893353632662507858", "0.00168891240733499999005657134034105162084838396146200001977586302588731457001767060904451498028490218780736"),
    )),
    (spheroid = :oblate, c = 1//1, values = (
        ("-0.348602399470269089476372382847354292920767569504992053128206743341920042616829943858732111881131616819626", "0.0"),
        ("0.998690920312602560167567655815155285371540959667570014228333086152871840912665793146499415524856240728587", "0.0"),
        ("0.114368082844577740868710121249893286958441919006253975291841546598485625285558452594208271740939993729786", "0.0"),
        ("0.00197634822528126305380979457907817434263613938768717370029283881232848679882915599284358230355737440602016", "0.0"),
    )),
    (spheroid = :oblate, c = complex(1//1, 1//4), values = (
        ("-0.321883129754014593679374588303615476083721858296693323713886315177274424476460712582066203735686938450904", "-0.181098665410776759271204693665910000710095166563737449848924342434395006230367264607585795845087562112543"),
        ("0.999211790976822051132071318400150944955618245192655948251885647984138626834897164355932986079882011953378", "-0.0012466304160699010221888830556866920884166776267390894390276452601373686377073913097717176178869900706233"),
        ("0.106371191819943720049942512226922782525475164229014814432442529491013545059396559715937937306410845992993", "0.0585401058362556729545925995067605926520193604226642728190496301812162755055449908202610184906976350818728"),
        ("0.00121019315064169247997324166390006152131463108695380599579427965859450928561079289903373978839647900839429", "0.00187142947799564054500395173155620491928190628022222976801600097203505621223804469437970507254463970605194"),
    )),
)

@testset "Requested arbitrary precision" begin
    original = (precision(BigFloat), rounding(BigFloat))
    setprecision(BigFloat, 512) do
        @testset "Independent values beyond quad" begin
            for fixture in ARBITRARY_PRECISION_REFERENCES
                (; spheroid, c) = fixture
                reference = [complex(parse(BigFloat, re), parse(BigFloat, im))
                             for (re, im) in fixture.values]
                errors = BigFloat[]
                for bits in (160, 256)
                    result = dmn(0, 0, c; spheroid, precision_bits = bits)
                    lambda = eigenvalue(0, 0, c; spheroid, precision_bits = bits)
                    V = c isa Real ? BigFloat : Complex{BigFloat}
                    @test lambda isa V
                    @test eltype(result.coefficients) === V
                    @test result.converged
                    @test result.degrees[1:3] == [0, 2, 4]
                    @test result.terms == length(result.coefficients) == length(result.degrees)
                    @test result.terms == length(result.dcoefficients_dc)
                    for value in (lambda, result.eigenvalue, result.deigenvalue_dc,
                        result.coefficients..., result.dcoefficients_dc...,
                        result.tail_norm, result.relative_change)
                        @test precision(real(value)) == bits
                        value isa Complex && @test precision(imag(value)) == bits
                    end
                    @test isapprox(lambda, reference[1]; rtol = BigFloat(2)^(-bits+3))
                    @test isapprox(result.eigenvalue, reference[1]; rtol = BigFloat(2)^(-bits+3))
                    @test all(isapprox(a, b; rtol = BigFloat(2)^(-bits+4))
                              for (a, b) in zip(result.coefficients[1:3], reference[2:4]))
                    push!(errors, abs(lambda-reference[1]))
                end
                # Requested bits must improve measured accuracy beyond native precision.
                @test errors[1] < big"1e-45"
                @test errors[2] < errors[1]*BigFloat(2)^(-80)
                native = eigenvalue(0, 0, c; spheroid, precision = :quad)
                @test errors[1]*big"1e6" < abs(native-reference[1])
                # Either native precision selects the branch without capping requested accuracy.
                quad_seed = eigenvalue(0, 0, c; spheroid, precision = :quad,
                    precision_bits = 160)
                quad_coefficients = dmn(0, 0, c; spheroid, precision = :quad,
                    precision_bits = 160)
                @test precision(real(quad_seed)) == 160
                @test isapprox(quad_seed, reference[1]; rtol = BigFloat(2)^(-157))
                @test all(isapprox(a, b; rtol = BigFloat(2)^(-156)) for (a, b) in
                    zip(quad_coefficients.coefficients[1:3], reference[2:4]))
            end
        end

        @testset "Signed real parameters" begin
            for spheroid in (:prolate, :oblate)
                positive = dmn(1, 3, 2; spheroid, precision_bits = 160)
                negative = dmn(1, 3, -2; spheroid, precision_bits = 160)
                @test positive.degrees == negative.degrees
                @test positive.eigenvalue == negative.eigenvalue
                @test positive.coefficients == negative.coefficients
                @test positive.deigenvalue_dc == -negative.deigenvalue_dc
                @test positive.dcoefficients_dc == -negative.dcoefficients_dc
                @test eigenvalue(1, 3, -2; spheroid, precision_bits = 160) == positive.eigenvalue
                for value in (negative.eigenvalue, negative.deigenvalue_dc,
                    negative.coefficients..., negative.dcoefficients_dc...,
                    negative.tail_norm, negative.relative_change)
                    @test value isa BigFloat
                    @test precision(value) == 160
                end
            end
        end

        @testset "Normalization, bases, and parameter sensitivities" begin
            bits = 160
            tolerance = BigFloat(2)^(-bits+6)
            ferrers_norm(l) = BigFloat(2)*factorial(big(l+1))/((2l+1)*factorial(big(l-1)))
            for spheroid in (:prolate, :oblate), c in (5//4, complex(1//1, 1//4))
                plain = dmn(1, 2, c; spheroid, precision_bits = bits)
                unit = dmn(1, 2, c; spheroid, precision_bits = bits, normalize = true)
                orthogonal = dmn(1, 2, c; spheroid, precision_bits = bits,
                    normalize = true, basis = :orthonormal)
                norms = ferrers_norm.(unit.degrees)
                @test isapprox(sum(unit.coefficients.^2 .* norms), 1; rtol = tolerance)
                @test isapprox(sum(orthogonal.coefficients.^2), 1; rtol = tolerance)
                @test isapprox(unit.coefficients .* sqrt.(norms), orthogonal.coefficients;
                    rtol = tolerance)
                @test isapprox(unit.coefficients .* sqrt(ferrers_norm(2)), plain.coefficients;
                    rtol = tolerance)
                @test abs(sum(unit.coefficients .* unit.dcoefficients_dc .* norms)) < tolerance
                @test isapprox(unit.deigenvalue_dc, plain.deigenvalue_dc; rtol = tolerance)

                # Independent central differences have step and roundoff errors below 160 bits.
                h = BigFloat(2)^(-100)
                lower = dmn(1, 2, c-h; spheroid, precision_bits = 320, derivatives = false)
                upper = dmn(1, 2, c+h; spheroid, precision_bits = 320, derivatives = false)
                count = min(plain.terms, lower.terms, upper.terms)
                difference = (upper.coefficients[1:count]-lower.coefficients[1:count])/(2h)
                @test isapprox(plain.dcoefficients_dc[1:count], difference; rtol = tolerance)
                @test isapprox(plain.deigenvalue_dc,
                    (upper.eigenvalue-lower.eigenvalue)/(2h); rtol = tolerance)
                @test lower.dcoefficients_dc === nothing
                @test lower.deigenvalue_dc === nothing
            end
        end

        @testset "Large parameter refinement beyond 512 bits" begin
            result = dmn(0, 0, 50; precision_bits = 512)
            reference = dmn(0, 0, 50; precision_bits = 640)
            count = min(result.terms, reference.terms)
            tolerance = BigFloat(2)^(-508)
            @test result.converged
            @test precision(result.eigenvalue) == 512
            @test all(x -> precision(x) == 512, result.coefficients)
            @test result.tail_norm <= BigFloat(2)^(-512)
            @test result.relative_change <= BigFloat(2)^(-512)
            @test isapprox(result.eigenvalue, reference.eigenvalue; rtol = tolerance)
            @test isapprox(result.coefficients[1:count], reference.coefficients[1:count];
                rtol = tolerance)
            @test isapprox(result.dcoefficients_dc[1:count], reference.dcoefficients_dc[1:count];
                rtol = tolerance)
            @test isapprox(result.deigenvalue_dc, reference.deigenvalue_dc; rtol = tolerance)
        end

        @testset "Static limits and degree ranges" begin
            for spheroid in (:prolate, :oblate), c in (0, complex(0))
                result = dmn(2, 4, c; spheroid, precision_bits = 160)
                @test result.eigenvalue == eigenvalue(2, 4, c; spheroid, precision_bits = 160) == 20
                @test result.coefficients[findfirst(==(4), result.degrees)] == 1
                @test count(!iszero, result.coefficients) == 1
                @test all(iszero, result.dcoefficients_dc)
                @test iszero(result.deigenvalue_dc)
            end
            for spheroid in (:prolate, :oblate)
                results = dmn(1, 1:3, 5//4; spheroid, precision_bits = 160,
                    derivatives = false, basis = :orthonormal)
                @test length(results) == 3
                for (n, result) in zip(1:3, results)
                    single = dmn(1, n, 5//4; spheroid, precision_bits = 160,
                        derivatives = false, basis = :orthonormal)
                    @test result == single
                    @test result.coefficients !== single.coefficients
                    @test result.degrees !== single.degrees
                end
            end
        end

        @testset "Input fidelity, rounding, cache ownership, and tasks" begin
            # Preserve stored input bits above both native quad and ambient precision.
            delta = BigFloat(2)^(-180)
            c = 1+delta
            center = dmn(0, 0, 1; precision_bits = 256)
            perturbed = setprecision(BigFloat, 64) do
                before = (precision(BigFloat), rounding(BigFloat))
                result = dmn(0, 0, c; precision_bits = 256)
                @test (precision(BigFloat), rounding(BigFloat)) == before
                result
            end
            @test precision(perturbed.eigenvalue) == 256
            @test isapprox(perturbed.eigenvalue-center.eigenvalue,
                center.deigenvalue_dc*delta; rtol = BigFloat(2)^(-70))
            rational_value = setprecision(BigFloat, 64) do
                eigenvalue(0, 0, 1//3; precision_bits = 256)
            end
            @test isapprox(rational_value, eigenvalue(0, 0, BigFloat(1)/3; precision_bits = 256);
                rtol = BigFloat(2)^(-252))

            rounded = map((RoundDown, RoundUp, RoundDown)) do mode
                setrounding(BigFloat, mode) do
                    first = eigenvalue(0, 0, 1; precision_bits = 96)
                    second = eigenvalue(0, 0, 1; precision_bits = 96)
                    @test first == second
                    @test precision(first) == 96
                    @test rounding(BigFloat) === mode
                    first
                end
            end
            @test rounded[1] == rounded[3]
            @test rounded[1] < rounded[2]
            @test isapprox(rounded[1], center.eigenvalue; rtol = BigFloat(2)^(-92))
            @test isapprox(rounded[2], center.eigenvalue; rtol = BigFloat(2)^(-92))

            saved = copy(center.coefficients)
            fill!(center.coefficients, 0)
            @test dmn(0, 0, 1; precision_bits = 256).coefficients == saved
            settings = (precision(BigFloat), rounding(BigFloat))
            requested = (160, 256, 192, 320, 160, 256)
            reference = [eigenvalue(1, 2, 5//4; precision_bits = bits) for bits in requested]
            tasks = [Threads.@spawn begin
                         yield()
                         value = eigenvalue(1, 2, 5//4; precision_bits = bits)
                         (value, precision(value), precision(BigFloat), rounding(BigFloat))
                     end for bits in requested]
            for (task, bits, value) in zip(tasks, requested, reference)
                @test fetch(task) == (value, bits, settings...)
            end
            @test (precision(BigFloat), rounding(BigFloat)) == settings
        end

        @testset "Invalid requests and computation limits" begin
            for bits in (0, -1, 1, true, 160.0, :quad)
                @test_throws ArgumentError eigenvalue(0, 0, 1; precision_bits = bits)
                @test_throws ArgumentError dmn(0, 0, 1; precision_bits = bits)
            end
            for operator in (:concentration, :fourier)
                @test_throws ArgumentError eigenvalue(0, 0, 1; operator, precision_bits = 160)
            end
            @test_throws ArgumentError dmn(0, 0, 1; precision_bits = 160, basis = :invalid)
            @test_throws ArgumentError dmn(0, 0, 1; precision_bits = 160, rtol = 0)
            @test_throws ArgumentError dmn(0, 0, 1; precision_bits = 160, rtol = NaN)
            @test_throws r"max_terms" eigenvalue(0, 0, 1; precision_bits = 256, max_terms = 1)
            @test_throws r"max_terms" dmn(0, 0, 1; precision_bits = 256, max_terms = 1)
            # A loose coefficient tolerance must not weaken an explicit bit request.
            strict = dmn(0, 0, 1; precision_bits = 160)
            loose = dmn(0, 0, 1; precision_bits = 160, rtol = 0.1)
            @test isapprox(strict.coefficients, loose.coefficients; rtol = BigFloat(2)^(-155))
            # A tighter tolerance refines calculations while retaining the requested output width.
            tolerance = BigFloat(2)^(-250)
            tight = dmn(0, 0, 1; precision_bits = 160, rtol = tolerance)
            @test tight.tail_norm <= tolerance
            @test tight.relative_change <= tolerance
            @test tight.terms >= strict.terms
            @test precision(tight.eigenvalue) == 160
            @test precision(tight.tail_norm) == 160
            @test precision(tight.relative_change) == 160
            @test isapprox(tight.eigenvalue, strict.eigenvalue; rtol = BigFloat(2)^(-156))
            @test eigenvalue(0, 0, 1) isa Float64
            @test eigenvalue(0, 0, 1; precision = :quad) isa BigFloat
        end
    end
    @test (precision(BigFloat), rounding(BigFloat)) == original
end
