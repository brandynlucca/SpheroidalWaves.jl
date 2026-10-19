using SpheroidalWaves, Test, Libdl

@testset "Shared real degree expansions" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        T=precision===:quad ? BigFloat : Float64
        tol=precision===:quad ? big"1e-25" : big"3e-12"
        c=T(1.25)
        points=T.([-0.7, 0, 0.3])
        # Start above m to check native degree indexing, and use odd m to
        # check the Condon-Shortley phase in both output formats.
        for (normalize, scaled) in ((false, false), (true, true))
            block=smn(1, 2:3, c, points; spheroid, precision, normalize,
                scaled, logderivative = true, second_derivative = true)
            for (i, n) in enumerate(2:3)
                one=smn(1, n, c, points; spheroid, precision, normalize,
                    scaled, logderivative = true, second_derivative = true)
                for field in (:value, :derivative, :second_derivative)
                    a, b=getproperty(block, field), getproperty(one, field)
                    if scaled
                        @test a.mantissa[:, i]≈b.mantissa rtol=tol atol=tol
                        @test a.exponent[:, i]==b.exponent
                    else
                        @test a[:, i]≈b rtol=tol atol=tol
                    end
                end
                @test all(isapprox(a, b; rtol = tol, atol = tol, nans = true)
                for (a, b) in zip(block.logderivative[:, i], one.logderivative))
            end
        end
        x=T.(spheroid===:prolate ? [1.2, 2.0] : [0.0, 0.4, 2.0])
        for kind in 1:4, scaled in (false, true)

            block=rmn(1, 2:4, c, x; spheroid, precision, kind, scaled,
                logderivative = true, second_derivative = true)
            for (i, n) in enumerate(2:4)
                one=rmn(1, n, c, x; spheroid, precision, kind, scaled,
                    logderivative = true, second_derivative = true)
                for field in (:value, :derivative, :second_derivative)
                    a, b=getproperty(block, field), getproperty(one, field)
                    if scaled
                        @test a.mantissa[:, i]≈b.mantissa rtol=tol atol=tol
                        @test a.exponent[:, i]==b.exponent
                    else
                        @test a[:, i]≈b rtol=tol atol=tol
                    end
                end
                @test all(isapprox(a, b; rtol = tol, atol = tol, nans = true)
                for (a, b) in zip(block.logderivative[:, i], one.logderivative))
            end
        end
        # A mixed range keeps the refined low degrees and batches the rest.
        block=smn(0, 0:3, T(3), T(0.3); spheroid, precision, scaled = true)
        for (i, n) in enumerate(0:3)
            one=smn(0, n, T(3), T(0.3); spheroid, precision, scaled = true)
            @test block.value.mantissa[:, i]≈one.value.mantissa rtol=tol
            @test block.value.exponent[:, i]==one.value.exponent
        end
        # Independent pair consistency, both geometries.
        r1=rmn(1, 2:4, c, x; spheroid, precision)
        r2=rmn(1, 2:4, c, x; spheroid, precision, kind = 2)
        weight=c .* (x .^ 2 .+ (spheroid===:prolate ? -1 : 1))
        @test weight .* (r1.value .* r2.derivative-r1.derivative .* r2.value)≈ones(length(x), 3) rtol=tol atol=tol
        @test size(smn(1, 2:2, c, T(0.3); spheroid, precision).value)==(1, 1)
        symbol=precision===:quad ? :spheroidal_degrees_scaled_text :
               :spheroidal_degrees_scaled_double
        lib=SpheroidalWaves.backend_library(; precision)
        if Libdl.dlsym_e(SpheroidalWaves._require_backend_handle(lib), symbol)!=C_NULL
            @test SpheroidalWaves._shared_real_degree_range(1, 2:4, c, points, spheroid,
                precision, :angular, 0, true, false, false)!==nothing
        end
        # These values exceed Float64's exponent range. No reconstruction in
        # Float64 may occur before the scaled mantissas are returned.
        block=smn(200, 200:201, T(1), T(0.3); spheroid, precision, scaled = true)
        @test all(>(308), block.value.exponent)
        for (i, n) in enumerate(200:201)
            one=smn(200, n, T(1), T(0.3); spheroid, precision, scaled = true)
            @test block.value.mantissa[:, i]≈one.value.mantissa rtol=tol
            @test block.value.exponent[:, i]==one.value.exponent
        end
    end
end

@testset "Combined evaluations and native coefficient blocks" begin
    SW = SpheroidalWaves
    function coefficient_error(a, b)
        count = min(length(a), length(b))
        sqrt(sum(abs2, a[1:count]-b[1:count]) + sum(abs2, a[(count + 1):end]) +
             sum(abs2, b[(count + 1):end]))/sqrt(sum(abs2, b))
    end
    for precision in (:double, :quad)
        tolerance = precision === :quad ? big"1e-26" : 3e-12
        for derivatives in (false, true), kind in 1:4

            result = spheroidal(1, 2:4, 1.25; eta = [-0.4, 0.25], x = [1.2, 2.0],
                coefficients = true, precision, normalize = true, kind, derivatives,
                basis = :orthonormal)
            angular = smn(1, 2:4, 1.25, [-0.4, 0.25]; precision, normalize = true)
            radial = rmn(1, 2:4, 1.25, [1.2, 2.0]; precision, kind)
            @test isapprox(result.angular.value, angular.value; rtol = tolerance)
            @test isapprox(result.radial.value, radial.value; rtol = tolerance)
            @test hasproperty(result.angular, :derivative) == derivatives
            @test hasproperty(result.radial, :derivative) == derivatives
            if derivatives
                @test isapprox(result.angular.derivative, angular.derivative; rtol = tolerance)
                @test isapprox(result.radial.derivative, radial.derivative; rtol = tolerance)
            end
            @test length(result.coefficients) == 3
            for (i, degree) in enumerate(2:4)
                actual = result.coefficients[i]
                reference = dmn(
                    1, degree, 1.25; precision, normalize = true, basis = :orthonormal)
                @test coefficient_error(actual.coefficients, reference.coefficients) <
                      tolerance
                @test isapprox(sum(abs2, actual.coefficients), 1; rtol = tolerance)
                @test actual.dcoefficients_dc === nothing
                @test actual.deigenvalue_dc === nothing
            end
        end
        # Guarded low-degree angular functions retain their independent MPFR path.
        mixed = spheroidal(1, 1:5, 4.0; eta = [0, 0.4], x = 1.2, precision,
            normalize = true, derivatives = false)
        @test isapprox(mixed.angular.value,
            smn(1, 1:5, 4.0, [0, 0.4]; precision, normalize = true).value; rtol = tolerance)
        # Exact parity zeros must not shorten the next degree's native series.
        for (m, n, c) in ((0, 0:74, 21.0), (19, 19:100, 90.0))
            broadside = spheroidal(m, n, c; eta = 0.0, precision,
                normalize = true, derivatives = false)
            @test isapprox(broadside.angular.value,
                smn(m, n, c, 0.0; precision, normalize = true).value; rtol = tolerance)
        end
        # Scalar degree output shapes agree with smn/rmn, including endpoints.
        for eta in ([-1.0, 0.25, 1.0], [0.25+0.125im])
            combined = spheroidal(1, 2, 1.25; eta, x = 1.2, precision,
                coefficients = true, parameter_derivatives = true)
            @test isapprox(combined.angular.value, smn(1, 2, 1.25, eta; precision).value; rtol = tolerance)
            @test combined.radial.value isa Vector
            @test combined.coefficients.dcoefficients_dc !== nothing
        end
        for spheroid in (:prolate, :oblate), c in (0.0, 1.25+0.125im)

            result = spheroidal(1, 2, c; eta = 0.25, coefficients = true,
                spheroid, precision, normalize = true, basis = :orthonormal)
            @test isapprox(result.angular.value,
                smn(1, 2, c, 0.25; spheroid, precision, normalize = true).value; rtol = tolerance)
            reference = dmn(
                1, 2, c; spheroid, precision, normalize = true, basis = :orthonormal)
            @test isapprox(result.coefficients.coefficients, reference.coefficients; rtol = tolerance)
        end
        # The static coefficient-only request needs no native wave computation.
        static = spheroidal(
            1, 2, 0; coefficients = true, parameter_derivatives = true, precision)
        @test keys(static) == (:coefficients,)
        @test all(iszero, static.coefficients.dcoefficients_dc)
        # Reconstruct with the independently normalized static Ferrers basis.
        expansion = dmn(1, 2, 1.25; precision, normalize = true,
            derivatives = false, basis = :orthonormal)
        points = [-0.4, 0.25, 0.7]
        reconstructed = sum(d .* smn(1, l, 0, points; precision, normalize = true).value
        for (d, l) in zip(expansion.coefficients, expansion.degrees))
        @test isapprox(
            reconstructed, smn(1, 2, 1.25, points;
                precision, normalize = true).value;
            rtol = tolerance)
        # Native coefficients are not cached returned arrays.
        block = dmn(1, 2:4, 1.25; precision, derivatives = false, basis = :orthonormal)
        saved = copy(block[1].coefficients)
        fill!(block[1].coefficients, 0)
        @test dmn(1, 2:4, 1.25; precision, derivatives = false,
            basis = :orthonormal)[1].coefficients == saved
        pointer = Libdl.dlsym_e(
            SW._require_backend_handle(SW.backend_library(; precision)),
            :spheroidal_prolate_block_text)
        if pointer != C_NULL
            native = SW._native_prolate_coefficient_block(
                1, 2:4, 1.25, :prolate, precision,
                true, :orthonormal, nothing, 512)
            @test all(!isnothing, native)
            @test all(d -> isnan(d.relative_change), native)
        end
    end
    # The basis changes representation, including sensitivities, but not S.
    for spheroid in (:prolate, :oblate)
        raw = dmn(1, 2, 1.25; spheroid, precision = :quad)
        orthogonal = dmn(1, 2, 1.25; spheroid, precision = :quad, basis = :orthonormal)
        factors = [sqrt(BigFloat(2)*l*(l+1)/(2l+1)) for l in raw.degrees]
        @test isapprox(orthogonal.coefficients, raw.coefficients .* factors; rtol = big"1e-65")
        @test isapprox(orthogonal.dcoefficients_dc, raw.dcoefficients_dc .* factors; rtol = big"1e-65")
    end
    # A strict requested tolerance must reject fixed native precision.
    strict = dmn(1, 2, 1.25; precision = :quad, derivatives = false,
        basis = :orthonormal, rtol = big"1e-45")
    @test isfinite(strict.relative_change)
    @test strict.relative_change <= big"1e-45"
    @test strict.tail_norm <= big"1e-45"
    # Large-parameter coefficients retain separate guarded pointwise angular evaluation.
    for (m, n, c) in ((0, 0, 90.0), (19, 52, 250.0), (38, 120, 250.0))
        native = dmn(m, n, c; precision = :quad, derivatives = false,
            normalize = true, basis = :orthonormal)
        reference = dmn(m, n, c; precision = :quad, normalize = true, basis = :orthonormal)
        @test coefficient_error(native.coefficients, reference.coefficients) < big"1e-27"
    end
    high = spheroidal(200, 200:201, 1.0; eta = 0.3, precision = :quad)
    @test isapprox(high.angular.value, smn(200, 200:201, 1.0, 0.3; precision = :quad).value; rtol = big"1e-26")
    @test maximum(abs, high.angular.value) > big"1e308"
    @test_throws ArgumentError spheroidal(0, 0, 1)
    @test_throws ArgumentError spheroidal(0, 0, 1; eta = Float64[])
    @test_throws ArgumentError dmn(0, 0, 1; basis = :invalid)
    @test_throws ArgumentError dmn(0, 0:2, 1; rtol = 0)
    @test_throws ArgumentError dmn(0, 4:5, 1; max_terms = 1)
    # Calls share only thread-private native caches. Return data are call-local.
    evaluate(c) = spheroidal(1, 2:4, c; eta = [0.25, -0.4], x = 1.2,
        coefficients = true, precision = :quad, normalize = true, kind = 3,
        derivatives = false, basis = :orthonormal)
    parameters = [1.25, 2.25, 1.25, 2.25]
    reference = evaluate.(parameters)
    concurrent = fetch.([Threads.@spawn evaluate(c) for c in parameters])
    for (actual, expected) in zip(concurrent, reference)
        @test actual.angular == expected.angular
        @test actual.radial == expected.radial
        @test [d.coefficients for d in actual.coefficients] ==
              [d.coefficients for d in expected.coefficients]
    end
end

@testset "Coefficient degree ranges preserve sensitivities and fallback" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        tolerance = precision === :quad ? big"1e-26" : 3e-12
        # Default parameter sensitivities require the guarded solver for every degree.
        results = dmn(1, 2:4, 1.25; spheroid, precision, normalize = true)
        @test length(results) == 3
        for (degree, actual) in zip(2:4, results)
            reference = dmn(1, degree, 1.25; spheroid, precision, normalize = true)
            @test actual.degrees == reference.degrees
            for field in (:coefficients, :dcoefficients_dc, :eigenvalue, :deigenvalue_dc)
                @test isapprox(getproperty(actual, field), getproperty(reference, field);
                    rtol = tolerance)
            end
        end
        # Returned sensitivities must not alias the cached expansion or another degree.
        saved = deepcopy(results)
        fill!(results[1].coefficients, 0)
        fill!(results[1].dcoefficients_dc, 0)
        @test results[2] == saved[2]
        @test dmn(1, 2:4, 1.25; spheroid, precision, normalize = true) == saved
    end
    for precision in (:double, :quad)
        # Fixed native precision cannot meet this tolerance. Each degree must fall back.
        results = dmn(1, 2:4, 1.25; precision, derivatives = false,
            basis = :orthonormal, rtol = big"1e-45")
        for (degree, actual) in zip(2:4, results)
            reference = dmn(1, degree, 1.25; precision,
                basis = :orthonormal, rtol = big"1e-45")
            @test isfinite(actual.relative_change)
            @test actual.relative_change <= big"1e-45"
            @test actual.tail_norm <= big"1e-45"
            @test actual.dcoefficients_dc === nothing
            @test actual.deigenvalue_dc === nothing
            # Sensitivity convergence can retain more terms than values-only evaluation.
            count = min(actual.terms, reference.terms)
            @test actual.degrees[1:count] == reference.degrees[1:count]
            difference = sum(abs2, actual.coefficients[1:count]-reference.coefficients[1:count]) +
                         sum(abs2, actual.coefficients[(count + 1):end]) +
                         sum(abs2, reference.coefficients[(count + 1):end])
            @test sqrt(difference/sum(abs2, reference.coefficients)) <= big"1e-45"
            @test actual.eigenvalue == reference.eigenvalue
        end
    end
end

@testset "Combined oblate geometry dispatch" begin
    for precision in (:double, :quad)
        tolerance = precision === :quad ? big"1e-26" : 3e-12
        for derivatives in (false, true), kind in 1:4
            # Oblate coordinates include the origin and values below the prolate boundary.
            points = [0.0, 0.4]
            actual = spheroidal(1, 2:4, 1.25; spheroid = :oblate,
                eta = [-0.4, 0.25], x = points, coefficients = true,
                precision, kind, derivatives, normalize = true, basis = :orthonormal)
            angular = smn(1, 2:4, 1.25, [-0.4, 0.25];
                spheroid = :oblate, precision, normalize = true)
            radial = rmn(1, 2:4, 1.25, points; spheroid = :oblate, precision, kind)
            @test size(actual.angular.value) == (2, 3)
            @test size(actual.radial.value) == (2, 3)
            @test all(isfinite, actual.radial.value)
            @test isapprox(actual.angular.value, angular.value; rtol = tolerance)
            @test isapprox(actual.radial.value, radial.value; rtol = tolerance)
            @test hasproperty(actual.angular, :derivative) == derivatives
            @test hasproperty(actual.radial, :derivative) == derivatives
            if derivatives
                @test isapprox(actual.angular.derivative, angular.derivative; rtol = tolerance)
                @test isapprox(actual.radial.derivative, radial.derivative; rtol = tolerance)
            end
            # Coefficients must retain the oblate eigenproblem and requested basis.
            for (i, degree) in enumerate(2:4)
                reference = dmn(1, degree, 1.25; spheroid = :oblate, precision,
                    normalize = true, basis = :orthonormal, derivatives = false)
                @test isapprox(actual.coefficients[i].coefficients,
                    reference.coefficients; rtol = tolerance)
                @test isapprox(actual.coefficients[i].eigenvalue,
                    reference.eigenvalue; rtol = tolerance)
                @test actual.coefficients[i].dcoefficients_dc === nothing
            end
        end
        for c in (1.25, 1.25+0.125im)
            # Scalar-degree fallback preserves complex-parameter and sensitivity conventions.
            actual = spheroidal(1, 2, c; spheroid = :oblate, eta = 0.25, x = 0.4,
                coefficients = true, precision, kind = 3, parameter_derivatives = true)
            angular = smn(1, 2, c, 0.25; spheroid = :oblate, precision)
            radial = rmn(1, 2, c, 0.4; spheroid = :oblate, precision, kind = 3)
            coefficients = dmn(1, 2, c; spheroid = :oblate, precision)
            @test actual.angular.value isa Vector
            @test actual.radial.value isa Vector
            for field in (:value, :derivative)
                @test isapprox(getproperty(actual.angular, field),
                    getproperty(angular, field); rtol = tolerance)
                @test isapprox(getproperty(actual.radial, field),
                    getproperty(radial, field); rtol = tolerance)
            end
            @test isapprox(actual.coefficients.coefficients,
                coefficients.coefficients; rtol = tolerance)
            @test isapprox(actual.coefficients.dcoefficients_dc,
                coefficients.dcoefficients_dc; rtol = tolerance)
            @test isapprox(actual.coefficients.deigenvalue_dc,
                coefficients.deigenvalue_dc; rtol = tolerance)
            # Omitted output families remain absent on the oblate fallback.
            radial_only = spheroidal(1, 2, c; spheroid = :oblate, x = 0.4,
                precision, kind = 3, derivatives = false)
            @test keys(radial_only) == (:radial,)
            @test keys(radial_only.radial) == (:value,)
            @test isapprox(radial_only.radial.value, radial.value; rtol = tolerance)
        end
    end
    # Each geometry retains its own coordinate validation rather than sharing prolate bounds.
    @test_throws ArgumentError spheroidal(1, 2, 1.25; spheroid = :invalid, eta = 0.25)
    @test_throws ErrorException spheroidal(1, 2, 1.25; spheroid = :prolate, x = 0.4)
    @test_throws ErrorException spheroidal(1, 2, 1.25; spheroid = :oblate, x = -0.4)
end
