using SpheroidalWaves, Test

const PRECISION_WAVE_CASE = get(ENV, "SPHEROIDALWAVES_PRECISION_CASE", "all")
precision_wave_case(name) = PRECISION_WAVE_CASE=="all" || name in split(PRECISION_WAVE_CASE, ',')

function test_wave_bits(result, bits)
    for channel in values(result)
        channel isa NamedTuple && (channel = channel.mantissa)
        for value in channel
            @test precision(real(value))==bits
            value isa Complex && @test precision(imag(value))==bits
        end
    end
end

@testset "Arbitrary-precision waves" begin
    original = (precision(BigFloat), rounding(BigFloat))
    setprecision(BigFloat, 512) do
        if precision_wave_case("references") || startswith(PRECISION_WAVE_CASE, "reference-")
            @testset "Independent wave references" begin
                for line in eachline(joinpath(@__DIR__, "precision_wave_references.tsv"))
                    startswith(line, '#') && continue
                    family, geometry, ci, zi, yr, yi, dr, di = split(line, '\t')
                    name = "reference-$family-$geometry-$ci-$zi"
                    startswith(PRECISION_WAVE_CASE, "reference-") && PRECISION_WAVE_CASE!=name && continue
                    spheroid = Symbol(geometry)
                    angular = startswith(family, "angular")
                    kind = parse(Int, string(last(family)))
                    parameter = ci=="0" ? 1 : complex(1//1, 1//4)
                    coordinate = angular ? 1//4 : 2//1
                    zi=="0" || (coordinate = complex(coordinate, 1//8))
                    reference = complex.(parse.(BigFloat, [yr, dr]), parse.(BigFloat, [yi, di]))
                    evaluate = angular ? smn : rmn
                    errors = BigFloat[]
                    @testset "$name" begin
                        for bits in (160, 256)
                            result = evaluate(0, 0, parameter, [coordinate]; spheroid, kind,
                                precision_bits = bits)
                            test_wave_bits(result, bits)
                            actual = [only(result.value), only(result.derivative)]
                            @test all(isapprox(a, b; rtol = BigFloat(2)^(-bits+5))
                                for (a, b) in zip(actual, reference))
                            push!(errors, maximum(abs.((actual.-reference)./reference)))
                        end
                        @test errors[1]<BigFloat(2)^(-150)
                        @test errors[2]<errors[1]*BigFloat(2)^(-75)
                    end
                end
            end
        end

        if precision_wave_case("derivatives")
            @testset "Coordinate derivatives and normalization" begin
                # Independent mpmath Rodrigues expansion at 140 decimal digits and 40 terms.
                reference = parse.(BigFloat, [
                    "-0.9810303645284801495208013853358912324341863677607980647850523898460369053711354359338488325247887408",
                    "0.309649688659690203591415049178354606857683011389599290408852199466018232726208293765243320492869139",
                    "1.281042838578304673704606148612255930974980839330677629289488083999061630737295628790717893306030985",
                    "0.5562933351230701785927347286923792401927600691399871250449318963101738249606263584279011241969874035",
                    "3.228613384213223902151693627395703382660261172217180154289260380705731058780673546765210536720077844"])
                for bits in (160, 256)
                    plain = smn(1, 1, 1, [1//4]; precision_bits = bits, derivatives = 4)
                    test_wave_bits(plain, bits)
                    @test all(isapprox(only(channel), value; rtol = BigFloat(2)^(-bits+5))
                        for (channel, value) in zip(values(plain), reference))
                    unit = smn(1, 1, 1, [1//4]; precision_bits = bits, normalize = true,
                        derivatives = 4)
                    @test all(isapprox(only(channel), value*sqrt(BigFloat(3)/4);
                        rtol = BigFloat(2)^(-bits+5)) for (channel, value) in zip(values(unit), reference))
                    scaled = smn(1, 1, 1, [1//4]; precision_bits = bits, derivatives = 4,
                        scaled = true, logderivative = true)
                    test_wave_bits(scaled, bits)
                    for name in keys(plain)
                        channel = getproperty(scaled, name)
                        @test isapprox(only(channel.mantissa)*BigFloat(10)^only(channel.exponent),
                            only(getproperty(plain, name)); rtol = BigFloat(2)^(-bits+5))
                    end
                    @test isapprox(only(scaled.logderivative), reference[2]/reference[1];
                        rtol = BigFloat(2)^(-bits+5))
                end
            end
        end

        if precision_wave_case("static")
            @testset "Static solutions and endpoint limits" begin
                bits = 160
                tolerance = BigFloat(2)^(-bits+6)
                for spheroid in (:prolate, :oblate)
                    p = smn(0, 3, 0, [-1, 1]; spheroid, precision_bits = bits, derivatives = 4)
                    expected = ([-1, 1], [6, 6], [-15, 15], [15, 15], [0, 0])
                    @test all(isapprox(a, b; rtol = tolerance, atol = tolerance)
                        for (a, b) in zip(values(p), expected))
                    test_wave_bits(p, bits)
                    p1 = smn(1, 1, 0, [-1, 1]; spheroid, precision_bits = bits, derivatives = 4)
                    @test p1.value==[0, 0]
                    @test p1.derivative==[-Inf, Inf]
                    @test p1.second_derivative==[Inf, Inf]
                    @test p1.third_derivative==[-Inf, Inf]
                    @test p1.fourth_derivative==[Inf, Inf]
                    qend = smn(0, 0, 0, [-1, 1]; spheroid, kind = 2,
                        precision_bits = bits, derivatives = 4)
                    @test qend.value==[-Inf, Inf]
                    @test qend.derivative==[Inf, Inf]
                    @test qend.second_derivative==[-Inf, Inf]
                    @test qend.third_derivative==[Inf, Inf]
                    @test qend.fourth_derivative==[-Inf, Inf]
                    for z in (BigFloat(1)/4, complex(BigFloat(1)/4, BigFloat(1)/8))
                        q = smn(0, 0, 0, [z]; spheroid, kind = 2,
                            precision_bits = bits, derivatives = 4)
                        u = 1-z^2
                        exact = (atanh(z), inv(u), 2z/u^2, 2(1+3z^2)/u^3, 24z*(1+z^2)/u^4)
                        @test all(isapprox(only(a), b; rtol = tolerance) for (a, b) in zip(values(q), exact))
                        test_wave_bits(q, bits)
                    end
                    for z in (BigFloat(2), complex(BigFloat(2), BigFloat(1)/8))
                        regular = rmn(0, 0, 0, [z]; spheroid, normalization = :static,
                            precision_bits = bits, derivatives = 4)
                        @test only(regular.value)==1
                        @test all(iszero(only(channel)) for channel in values(regular)[2:end])
                        irregular = rmn(0, 0, 0, [z]; spheroid, normalization = :static,
                            precision_bits = bits, kind = 2)
                        expected_value = spheroid===:prolate ? -atanh(inv(z)) : -atan(inv(z))
                        @test isapprox(only(irregular.value), expected_value; rtol = tolerance)
                        @test isapprox(only(irregular.derivative), inv(z^2-(spheroid===:prolate ? 1 : -1)); rtol = tolerance)
                        test_wave_bits(irregular, bits)
                    end
                end
            end
        end

        if precision_wave_case("sensitivities")
            @testset "Analytic parameter sensitivities" begin
                bits = 160
                h = BigFloat(2)^(-90)
                for spheroid in (:prolate, :oblate), kind in (1, 2), angular in (true, false)
                    parameter = spheroid===:prolate ? complex(BigFloat(1), BigFloat(1)/4) : BigFloat(1)
                    coordinate = complex(angular ? BigFloat(1)/4 : BigFloat(2), BigFloat(1)/8)
                    evaluate, jacobian = angular ? (smn, jacobian_smn) : (rmn, jacobian_rmn)
                    derivative = jacobian(0, 0, parameter, [coordinate]; spheroid, kind,
                        precision_bits = bits, diagnostics = true)
                    lower = evaluate(0, 0, parameter-h, [coordinate]; spheroid, kind,
                        precision_bits = 320)
                    upper = evaluate(0, 0, parameter+h, [coordinate]; spheroid, kind,
                        precision_bits = 320)
                    names = parameter isa Real ? (:dvalue_dc, :dderivative_dc) :
                        (:dvalue_dcreal, :dderivative_dcreal)
                    for (name, channel) in zip(names, (:value, :derivative))
                        actual = only(getproperty(derivative, name))
                        expected = (only(getproperty(upper, channel))-only(getproperty(lower, channel)))/(2h)
                        @test isapprox(actual, expected; rtol = BigFloat(2)^(-bits+7))
                        @test precision(real(actual))==bits
                        @test precision(imag(actual))==bits
                    end
                    if parameter isa Complex
                        @test derivative.dvalue_dcimag==im .* derivative.dvalue_dcreal
                        @test derivative.dderivative_dcimag==im .* derivative.dderivative_dcreal
                        @test derivative.metadata_value_dcreal.finite_flag
                    else
                        @test derivative.metadata_value.finite_flag
                    end
                end
            end
        end

        if precision_wave_case("dispatch")
            @testset "Ranges, combined calls, and backend independence" begin
                bits = 160
                for spheroid in (:prolate, :oblate)
                    ranged = smn(0, 0:1, 7//8, [1//4]; spheroid, precision_bits = bits)
                    combined = spheroidal(0, 0:1, 7//8; eta = 1//4, x = 2,
                        coefficients = true, spheroid, precision_bits = bits)
                    @test ranged==combined.angular
                    @test size(combined.radial.value)==(1, 2)
                    @test length(combined.coefficients)==2
                    for n in 0:1
                        scalar = smn(0, n, 7//8, [1//4]; spheroid, precision_bits = bits)
                        @test ranged.value[:, n+1]==scalar.value
                    end
                end
                SW = SpheroidalWaves
                saved = copy(SW._backend_libraries)
                try
                    SW._backend_libraries[:double] = SW._backend_libraries[:quad] = nothing
                    @test_throws r"No backend library configured" SW._require_backend_library(:double)
                    @test_throws r"No backend library configured" SW._require_backend_library(:quad)
                    for spheroid in (:prolate, :oblate)
                        c = complex(9//8, 1//8)
                        angular = smn(0, 0, c, [complex(1//4, 1//8)]; spheroid, precision_bits = bits)
                        radial = rmn(0, 0, c, [2]; spheroid, precision_bits = bits)
                        @test all(isfinite, angular.value)
                        @test all(isfinite, radial.value)
                        @test isfinite(eigenvalue(0, 0, c; spheroid, precision_bits = bits))
                        @test dmn(0, 0, c; spheroid, precision_bits = bits).converged
                    end
                finally
                    merge!(SW._backend_libraries, saved)
                end
            end
        end

        if precision_wave_case("storage")
            @testset "Real coordinates with complex Legendre coefficients" begin
                bits = 160
                points = [-1//1, -1//2, 0//1, 1//2, 1//1]
                for spheroid in (:prolate, :oblate), mode in (RoundNearest, RoundDown, RoundUp)
                    setrounding(BigFloat, mode) do
                        parameter = complex(1//1, 1//8)
                        real_points = smn(2, 3, parameter, points;
                            spheroid, precision_bits = bits, derivatives = 4)
                        complex_points = smn(2, 3, parameter, complex.(points);
                            spheroid, precision_bits = bits, derivatives = 4)
                        @test all(channel -> eltype(channel)==Complex{BigFloat}, values(real_points))
                        @test all(isapprox(a, b; rtol = BigFloat(2)^(-bits+5), atol = BigFloat(2)^(-bits+5))
                            for (left, right) in zip(values(real_points), values(complex_points))
                            for (a, b) in zip(left, right))
                        @test real_points.value[[1, 3, 5]]==[0, 0, 0]
                        @test all(precision(real(value))==bits && precision(imag(value))==bits
                            for channel in values(real_points) for value in channel)
                        @test rounding(BigFloat)===mode
                    end
                end
            end

            @testset "Radial propagation anchor in coordinate batches" begin
                SW = SpheroidalWaves
                complex_parameter = complex(1//1, 1//8)
                for (spheroid, parameter, kinds) in ((:prolate, 1//1, (1, 2)),
                        (:prolate, complex_parameter, (1,)), (:oblate, 1//1, (1, 3)),
                        (:oblate, complex_parameter, (1, 2, 4)))
                    pair = SW._solve_coefficient_eigenpair_at_precision(0, 1, parameter,
                        spheroid, :double, 160, nothing, 512, true)
                    SW._with_swprecision(precision(real(pair.lambda))) do
                        for kind in kinds
                            points = SW._SWFloat[2, 3//2, 2]
                            batch = SW._radial_expansion_batch(pair, points, kind)
                            anchor = SW._radial_expansion_data(pair, SW._SWFloat(2), kind)
                            interior = SW._radial_expansion_batch(pair, SW._SWFloat[3//2], kind)
                            @test batch.states[1]==anchor.state==batch.states[3]
                            @test batch.states[2]==only(interior.states)
                            @test batch.tail==max(anchor.relative_tail, interior.tail)
                            @test all(isfinite(value) for state in batch.states for value in state)
                        end
                    end
                end
            end
        end

        if precision_wave_case("settings")
            @testset "Input bits, settings, and invalid requests" begin
                bits = 160
                x = 1-BigFloat(2)^(-180)
                value = setprecision(BigFloat, 64) do
                    before = (precision(BigFloat), rounding(BigFloat))
                    result = smn(1, 1, 0, [x]; precision_bits = bits)
                    @test (precision(BigFloat), rounding(BigFloat))==before
                    result
                end
                @test isapprox(only(value.value), -sqrt((1-x)*(1+x)); rtol = BigFloat(2)^(-bits+5))
                @test isapprox(only(value.derivative), x/sqrt((1-x)*(1+x)); rtol = BigFloat(2)^(-bits+5))
                test_wave_bits(value, bits)
                directed = map((RoundDown, RoundUp)) do mode
                    setrounding(BigFloat, mode) do
                        result = smn(0, 0, 1, [1//4]; precision_bits = bits)
                        @test rounding(BigFloat)===mode
                        only(result.value)
                    end
                end
                @test directed[1]<directed[2]
                for evaluate in (smn, rmn)
                    point = evaluate===smn ? complex(1//4, 1//8) : complex(2//1, 1//8)
                    @test_throws ErrorException evaluate(-1, 0, 1, [point]; precision_bits = bits)
                    @test_throws ArgumentError evaluate(0, 0, 1, [point]; precision_bits = bits, max_terms = 0)
                    @test_throws ArgumentError evaluate(0, 0, 1, [point]; precision_bits = 1)
                end
                @test_throws ArgumentError smn(0, 0, 1, [complex(1//4, 1//8)];
                    precision_bits = bits, kind = 2, normalize = true)
                @test_throws DomainError smn(0, 0, 1, [complex(2)]; precision_bits = bits)
                @test_throws r"max_terms" smn(0, 0, 1, [1//4]; precision_bits = 256, max_terms = 1)
            end
        end
    end
    @test (precision(BigFloat), rounding(BigFloat))==original
end
