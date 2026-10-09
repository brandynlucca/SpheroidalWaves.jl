using SpheroidalWaves, Test

@testset "Quad text and split-double transfer" begin
    SW = SpheroidalWaves
    setprecision(BigFloat, 256) do
        values = [big"1"+big"2"^(-80), -big"1.234567890123456789", big"0"]
        hi, lo = SW._split_real_vector_to_double_pairs(values)
        @test any(!iszero, lo)
        @test SW._combine_split_parts(hi, lo) ≈ values rtol=big"1e-31"
        @test SW._combine_split_complex_parts(hi, lo, -hi, -lo) ≈
              complex.(values, -values) rtol=big"1e-31"
        text = SW._format_fortran_input(values)
        @test SW._parse_fortran_output(text, length(values)) ≈ values rtol=big"1e-55"
        scalar = SW._format_fortran_input(values[1])
        @test SW._parse_fortran_output(scalar, 1) ≈ values[1:1] rtol=big"1e-55"
        @test_throws r"payload overflow" SW._format_fortran_input(big"1.25"; width = 1)
        @test_throws r"payload overflow" SW._format_fortran_input(values; width = 1)
        @test isempty(SW._parse_fortran_output(UInt8[], 0))
        # Conversion to Float64 can round a mantissa up to ten.
        scaled = SW._decimal_scaled([big"9.9999999999999999999", big"0", BigFloat(Inf)], Float64)
        @test scaled.mantissa == [1.0, 0.0, Inf]
        @test scaled.exponent == [1, 0, 0]
    end
end

@testset "Quad calculations accept ordinary numeric inputs" begin
    lib = SpheroidalWaves.backend_library(; precision = :quad)
    if lib !== nothing && isfile(lib)
        # Exactly representable inputs isolate calculation precision from input
        # rounding. An independent reference checks R1 for m=0, n=1, c=1, x=2.
        expected = big"0.45603690333332372100494242600202824193682868459324"
        r = rmn(0, 1, 1.0, 2.0; precision = :quad, kind = 1)
        @test eltype(r.value) === Complex{BigFloat}
        @test only(r.value) ≈ expected rtol=big"1e-28"
        @test abs(only(r.value)-expected) <
              abs(BigFloat(Float64(expected))-expected)/big"1e10"
        s = smn(1, 1, 0.0, 0.5; precision = :quad)
        @test only(s.value) ≈ -sqrt(big"0.75") rtol=big"1e-30"

        # Ordinary arrays must preserve the same quad values, derivatives, and
        # scaling exponents as explicitly widened versions of those inputs.
        for spheroid in (:prolate, :oblate)
            x = [1.5, 2.0]
            ordinary = rmn(
                1, 1:2, 1.25, x; precision = :quad, spheroid, kind = 3, scaled = true)
            widened = rmn(1, 1:2, big(1.25), BigFloat.(x); precision = :quad,
                spheroid, kind = 3, scaled = true)
            @test eltype(ordinary.value.mantissa) === Complex{BigFloat}
            @test ordinary == widened
        end
    else
        @info "Skipping ordinary-input quad checks: quad backend unavailable"
    end
end

@testset "Quad inputs accept higher BigFloat precision" begin
    lib = SpheroidalWaves.backend_library(; precision = :quad)
    if lib !== nothing && isfile(lib)
        for spheroid in (:prolate, :oblate), complex_parameter in (false, true)

            results = map((256, 768)) do bits
                setprecision(BigFloat, bits) do
                    c = big(5)/4
                    complex_parameter && (c = complex(c, big(1)/5))
                    (; angular = smn(0, 1, c, big(3)/10; precision = :quad, spheroid),
                        radial = rmn(
                            0, 1, c, big(4)/3; precision = :quad, spheroid, kind = 2),
                        lambda = eigenvalue(0, 1, c; precision = :quad, spheroid))
                end
            end
            for result in results[2:end]
                @test result.lambda ≈ results[1].lambda rtol=big"1e-31"
                for target in (:angular, :radial), field in (:value, :derivative)

                    @test getproperty(getproperty(result, target), field) ≈
                          getproperty(getproperty(results[1], target), field) rtol=big"1e-31"
                end
            end
        end
    else
        @info "Skipping higher-precision input checks: quad backend unavailable"
    end
end

@testset "Quad sweep preserves grid and branch decisions" begin
    # All three coordinates collapse to 1.0 in Float64.
    delta = big"1e-20"
    grid = [big"1", big"1"+delta, big"1"+2delta]
    evaluator(m, n, c) = c + n*delta
    for points in (grid, reverse(grid))
        s = eigenvalue_sweep(0, 0, points; precision = :quad, branch_lock = false,
            use_jacobian_predictor = false, evaluator)
        @test eltype(s.c) === BigFloat
        @test eltype(s.lambda) === BigFloat
        @test s.c == points
        @test s.lambda == points
    end
    # Candidate values collapse too; rounding would pick the first wrong degree.
    branch(m, n, c) = n == 1 ? big"1" : big"1" + (n+1)*delta
    s = eigenvalue_sweep(0, 1, grid; precision = :quad,
        use_jacobian_predictor = false, evaluator = branch)
    @test s.selected_n == [1, 1, 1]
    @test s.lambda == ones(BigFloat, 3)
    @test_throws ErrorException eigenvalue_sweep(0, 0, [Inf]; precision = :quad)
end

@testset "Complex quad precision survives native transfer" begin
    lib = SpheroidalWaves.backend_library(; precision = :quad)
    if lib === nothing || !isfile(lib)
        @info "Skipping complex quad transfer tests: backend unavailable."
    else
        delta = big"1e-20"
        for spheroid in (:prolate, :oblate)
            c = complex(big"1.25", big"0.2")
            eta = big"0.3"
            x = big"2"
            s = smn(0, 1, c, [eta, eta+delta]; precision = :quad, spheroid)
            r = rmn(0, 1, c, [x, x+delta]; precision = :quad, spheroid)
            @test eltype(s.value) === Complex{BigFloat}
            @test eltype(s.derivative) === Complex{BigFloat}
            @test eltype(r.value) === Complex{BigFloat}
            @test eltype(r.derivative) === Complex{BigFloat}
            # If either coordinates or outputs round to doubles, these slopes
            # are zero. Independently returned derivatives provide the check.
            @test (s.value[2]-s.value[1])/delta ≈ s.derivative[1] rtol=big"1e-10"
            @test (r.value[2]-r.value[1])/delta ≈ r.derivative[1] rtol=big"1e-10"
            for perturbation in (delta, im*delta)
                shifted = smn(0, 1, c+perturbation, [eta]; precision = :quad, spheroid)
                shifted_r = rmn(0, 1, c+perturbation, [x]; precision = :quad, spheroid)
                @test only(shifted.value) != s.value[1]
                @test only(shifted_r.value) != r.value[1]
                @test eigenvalue(0, 1, c+perturbation; precision = :quad, spheroid) !=
                      eigenvalue(0, 1, c; precision = :quad, spheroid)
            end
            eig = eigenvalue(0, 1, c; precision = :quad, spheroid)
            @test eig isa Complex{BigFloat}
            @test eig != Complex{BigFloat}(ComplexF64(eig))
            # Both native radial output channels must retain quad precision.
            # Direct Hankel waves have separate independent-reference tests.
            second = rmn(0, 1, c, [x]; precision = :quad, spheroid, kind = 2)
            expected = inv(c*(x^2+(spheroid === :prolate ? -1 : 1)))
            @test only(r.value[1] .* second.derivative-r.derivative[1] .* second.value) ≈
                  expected rtol=big"1e-27"
            # The complex route at real c agrees with the separately wrapped
            # real solver, at a tolerance that double precision cannot satisfy.
            for normalize in (false, true)
                real_s = smn(1, 1, big"1.25", [eta]; precision = :quad, spheroid, normalize)
                complex_s = smn(1, 1, complex(big"1.25", big"0"), [eta];
                    precision = :quad, spheroid, normalize)
                @test complex_s.value ≈ real_s.value rtol=big"1e-26"
                @test complex_s.derivative ≈ real_s.derivative rtol=big"1e-26"
            end
            @test eigenvalue(1, 2, complex(big"0", big"0"); precision = :quad, spheroid) isa
                  Complex{BigFloat}
            @test length(accuracy(
                0, 1, c, [eta]; precision = :quad, spheroid, target = :angular)) == 1
            radial_accuracy = accuracy(
                0, 1, c, [x]; precision = :quad, spheroid, target = :radial)
            @test length(radial_accuracy) == 1
            @test all(a -> -1 <= a <= 33, radial_accuracy)
            # Native mode=0 and the public sweep retain small grid spacings.
            grid = [big"1.25", big"1.25"+delta]
            sweep = eigenvalue_sweep(
                0, 1, grid; precision = :quad, spheroid, branch_lock = false)
            @test sweep.lambda ==
                  [eigenvalue(0, 1, z; precision = :quad, spheroid) for z in grid]
            @test sweep.lambda[1] != sweep.lambda[2]
        end
        # The rejected radial fallback must not return a double-valued result.
        @test_throws DomainError rmn(0, 1, big"0", [big(4)/3]; precision = :quad)
        old_library = SpheroidalWaves.backend_library(; precision = :double)
        if old_library !== nothing && isfile(old_library)
            @test !SpheroidalWaves._has_required_quad_abi(old_library)
            try
                SpheroidalWaves.set_backend_library!(old_library; precision = :quad)
                err = try
                    smn(0, 1, 1.0+0.2im, [0.3]; precision = :quad)
                catch e
                    e
                end
                @test err isa ErrorException
                @test occursin("Rebuild the native backend", sprint(showerror, err))
            finally
                SpheroidalWaves.set_backend_library!(lib; precision = :quad)
            end
        end
    end
end

@testset "Angular spherical limit and Condon–Shortley phase" begin
    # Closed-form Ferrers polynomials, independently of the implementation's
    # recurrence. Both normalization choices use the same phase.
    for precision in (:double, :quad), spheroid in (:prolate, :oblate),
        complex_c in (false, true), normalize in (false, true)
        T = precision === :quad ? BigFloat : Float64
        eta = T.([-0.5, 0.0, 0.25])
        u = 1 .- eta .^ 2
        c = complex_c ? complex(zero(T)) : zero(T)
        cases = (
            (0, 2, (3eta .^ 2 .- 1) ./ 2, 3eta),
            (1, 1, -sqrt.(u), eta ./ sqrt.(u)),
            (1, 2, -3eta .* sqrt.(u), -3 .* (1 .- 2eta .^ 2) ./ sqrt.(u)),
            (1, 3, -(3//2) .* (5eta .^ 2 .- 1) .* sqrt.(u),
                (3//2) .* eta .* (15eta .^ 2 .- 11) ./ sqrt.(u)),
            (1, 4, -(5//2) .* (7eta .^ 3 .- 3eta) .* sqrt.(u),
                (5//2) .* (28eta .^ 4 .- 27eta .^ 2 .+ 3) ./ sqrt.(u)),
            (2, 2, 3u, -6eta),
            (2, 3, 15eta .* u, 15 .* (1 .- 3eta .^ 2)),
            (3, 3, -15u .* sqrt.(u), 45eta .* sqrt.(u))
        )
        for (m, n, values, derivatives) in cases
            norm_squared = T(2) / (2n + 1) * T(factorial(n + m)) / factorial(n - m)
            scale = normalize ? inv(sqrt(norm_squared)) : one(T)
            result = smn(m, n, c, eta; precision, spheroid, normalize)
            expected_type = complex_c ? Complex{T} : T
            @test eltype(result.value) == expected_type
            @test eltype(result.derivative) == expected_type
            @test result.value≈scale .* values rtol=64eps(T) atol=64eps(T)
            @test result.derivative≈scale .* derivatives rtol=64eps(T) atol=64eps(T)
        end
    end

    # Exact endpoint limits: m=1 really has divergent derivatives; the other
    # listed cases have finite limits and must not produce 0/0 NaNs.
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        s = smn(0, 2, 0.0, [-1.0, 1.0]; precision, spheroid)
        @test s.value == [1, 1]
        @test s.derivative == [-3, 3]
        @test smn(1, 1, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [-Inf, Inf]
        @test smn(1, 2, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [Inf, Inf]
        @test smn(2, 2, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [6, -6]
        @test smn(3, 3, 0.0, [-1.0, 1.0]; precision, spheroid).derivative == [0, 0]
        @test smn(2, 2, 0.0, [-1.0, 1.0]; precision, spheroid, normalize = true).derivative ≈
              [6, -6] ./ sqrt(48 / 5)
    end

    setprecision(BigFloat, 256) do
        x = big"0.123456789012345678901234567890123456789"
        result = smn(1, 1, big"0", [x]; precision = :quad)
        @test abs(only(result.value) + sqrt(1 - x^2)) < big"1e-70"
        @test abs(only(result.derivative) - x / sqrt(1 - x^2)) < big"1e-70"
        # Unit normalization stays representable even when unnormalized P_m^m
        # would overflow double precision.
        m = 201
        expected = -sqrt(BigFloat(2m + 1) * binomial(big(2m), m) / big(2)^(2m + 1))
        normalized = smn(m, m, 0.0, 0.0; normalize = true)
        @test only(normalized.value) ≈ Float64(expected) rtol=1e-13
        @test only(normalized.derivative) == 0
    end

    @test_throws ErrorException smn(-1, 1, 0.0, 0.3)
    @test_throws ErrorException smn(2, 1, 0.0, 0.3)
    @test_throws ErrorException smn(1, 1, 0.0, 1.1)
    @test_throws ErrorException smn(1, 1, 0.0, NaN)
end

@testset "Regular endpoint second derivatives" begin
    SW = SpheroidalWaves
    setprecision(BigFloat, 256) do
        # P_3^2(x)=15x(1-x^2), so its endpoint slopes are -30 and its
        # second derivatives are +90 at -1 and -90 at +1.
        spherical = SW._evaluate_coefficient_vector((m = 2, degrees = [3]),
            [sqrt(big"240.0"/7)], [-1, 1]; second_derivative = true)
        @test spherical.value == [0, 0]
        @test spherical.derivative ≈ [-30, -30] rtol=big"1e-70"
        @test spherical.second_derivative ≈ [90, -90] rtol=big"1e-70"
        for precision in (:double, :quad), spheroid in (:prolate, :oblate)
            # At m=2, the endpoint limit of the angular equation gives
            # S''(x)=x*(lambda-sigma*c^2-3)*S'(x)/3 for x=+/-1.
            # c=5 selects coefficient evaluation through the public API.
            result = smn(2, 3, 5, [-1, 1]; precision, spheroid, second_derivative = true)
            lambda = eigenvalue(2, 3, 5; precision, spheroid)
            sigma = spheroid === :prolate ? 1 : -1
            expected = [-1, 1] .* (lambda-sigma*25-3) .* result.derivative ./ 3
            @test result.second_derivative ≈ expected rtol=(precision === :quad ?
                                                            big"1e-27" : 1e-12)
        end
        for precision in (:double, :quad), m in 0:5

            n, c = m+1, big"1.25"
            points = BigFloat[-1, -1 + big"1e-8", 1 - big"1e-8", 1]
            # A differentiated Legendre expansion provides an independent
            # reference for the native endpoint Taylor reconstruction.
            plan = SW._coefficient_plan(m, n, c; precision = :quad)
            scale = SW._coefficient_phase(plan, :quad)*sqrt(SW._ferrers_norm2(m, n, BigFloat))
            reference = SW._evaluate_coefficient_vector(plan, scale .* plan.v, points;
                second_derivative = true)
            result = smn(
                m, n, c, points; precision, second_derivative = true, scaled = true)
            for field in (:value, :derivative, :second_derivative)
                encoded = getproperty(result, field)
                actual = encoded.mantissa .* BigFloat(10) .^ encoded.exponent
                expected = getproperty(reference, field)
                @testset "$precision m=$m $field" begin
                    for i in eachindex(points)
                        @test actual[i] ≈ expected[i] rtol=(precision === :quad ?
                                                            big"1e-25" : 1e-11)
                    end
                end
            end
            # The regular radial solution has the same endpoint power.
            radial = rmn(
                m, n, c, BigFloat[1, 1 + big"1e-8"]; precision, second_derivative = true)
            if isodd(m) && m < 4
                @test isinf(real(radial.second_derivative[1]))
            elseif m > 4
                @test radial.second_derivative[1] == 0
            else
                @test isfinite(radial.second_derivative[1])
                @test radial.second_derivative[2] ≈ radial.second_derivative[1] rtol=1e-5
            end
        end
    end
end
