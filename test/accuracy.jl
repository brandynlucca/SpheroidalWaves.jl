using SpheroidalWaves, Test

@testset "Complex-coordinate convergence diagnostics" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate),
        target in (:angular, :radial)
        R = precision===:quad ? BigFloat : Float64
        points = target===:angular ? Complex{R}[0.25 + 0.125im, 0.5 + 0.25im] :
                 Complex{R}[2 + 0.25im, 1.5 + 0.25im]
        options = target===:angular ? (;) : (; normalization = :static)
        report = accuracy(0, 0, 0, points; spheroid, precision, target,
            kind = 2, diagnostics = true, options...)
        actual = (target===:angular ? smn : rmn)(
            0, 0, 0, points; spheroid, precision, kind = 2, options...)
        reference = target===:angular ? atanh.(points) :
                    spheroid===:prolate ? -atanh.(inv.(points)) : atan.(points) .- R(pi)/2
        @test report.digits ==
              accuracy(0, 0, 0, points; spheroid, precision, target, kind = 2, options...)
        @test all(d -> d>=(precision===:quad ? 28 : 12), report.digits)
        for (i, diagnostic) in enumerate(report.diagnostics)
            @test diagnostic.method === :continuation_refinement
            @test diagnostic.converged && diagnostic.finite_flag
            @test diagnostic.conditioning_flag === :good
            @test abs(actual.value[i]-reference[i]) <=
                  R(10)^(-report.digits[i])*abs(reference[i])
            @test diagnostic.coordinate_condition ≈
                  abs(points[i]*actual.derivative[i]/actual.value[i])
            @test iszero(diagnostic.parameter_condition)
        end
    end
    for target in (:angular, :radial)
        z = target===:angular ? 0.25+0.125im : 2+0.25im
        report = accuracy(1, 2, 1.25+0.125im, [z]; target, diagnostics = true)
        @test only(report.digits)>=12
        @test only(report.diagnostics).converged
    end
    zero_report = accuracy(0, 1, 0, [0im]; target = :angular, diagnostics = true)
    @test zero_report.digits == [-1]
    @test only(zero_report.diagnostics).finite_flag
    @test isinf(only(zero_report.diagnostics).coordinate_condition)
    endpoints = accuracy(
        1, 1, 1.25, [1+0im, 0.25+0.125im]; target = :angular, diagnostics = true)
    @test endpoints.digits[1] == -1 && endpoints.digits[2]>=12
    @test endpoints.diagnostics[1].conditioning_flag === :unavailable
    real_report = accuracy(0, 0, 0, [0.25]; target = :angular, diagnostics = true)
    @test real_report.digits == [-1] && real_report.diagnostics === nothing
    @test_throws ArgumentError accuracy(
        0, 0, 1, [0.25im]; target = :angular, kind = 2, normalize = true)
    @test accuracy(0, 0, 1, [2+0.25im]; normalize = true) == accuracy(0, 0, 1, [2+0.25im])

    SW = SpheroidalWaves
    F = SW._SWFloat
    SW._with_swprecision(384) do
        plan = SW._sensitivity_plan(0, 0, F(0), :prolate, :double)
        state = Complex{F}.((1, 1, 0, 0))
        unresolved = SW._coordinate_accuracy_record(
            state, Complex{F}.((2, 1, 0, 0)), plan, 1, :double, false)
        @test unresolved.digits == -1 && unresolved.diagnostics.conditioning_flag === :poor
        underflow = Complex{F}.((F("1e-400"), F("1e-400"), 0, 0))
        @test SW._coordinate_accuracy_record(
            underflow, underflow, plan, 1, :double, false).digits == 0
        overflow = Complex{F}.((F("1e400"), 1, 0, 0))
        @test SW._coordinate_accuracy_record(
            overflow, overflow, plan, 1, :double, false).diagnostics.conditioning_flag ===
              :singular
        sensitive = Complex{F}.((1, 10000, 0, 0))
        warning = SW._coordinate_accuracy_record(
            sensitive, sensitive, plan, 1, :double, false)
        @test warning.diagnostics.conditioning_flag === :warning
        @test warning.digits < 12
    end
end

@testset "Endpoint series estimates remain unavailable" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate)

        lib=SpheroidalWaves.backend_library(; precision)
        if lib === nothing || !isfile(lib)
            continue
        end
        T=precision === :quad ? BigFloat : Float64
        for c in (T(5)/4, complex(T(5)/4, T(1)/5)), m in (0, 2)

            estimates=accuracy(
                m, m+1, c, [-big"1", -1+big"1e-12", big"0.3", 1-big"1e-12", big"1"];
                spheroid, precision, target = :angular)
            @test estimates[[1, 2, 4, 5]] == [-1, -1, -1, -1]
            @test estimates[3] == only(accuracy(
                m, m+1, c, [big"0.3"]; spheroid, precision, target = :angular))
        end
    end
end

@testset "Accuracy reports unavailable estimates honestly" begin
    for precision in (:double, :quad), spheroid in (:prolate, :oblate),
        c in (0.0, 0.0+0.0im, big"0.0")
        @test accuracy(0, 2, c, [-1.0, 0.3, 1.0]; precision, spheroid, target = :angular) ==
              [-1, -1, -1]
        # Singular endpoint derivatives are not assigned fictitious digit counts.
        @test accuracy(1, 1, c, [-1.0, 1.0]; precision, spheroid, target = :angular) ==
              [-1, -1]
        for (m, n, eta) in ((-1, 0, [0.3]), (2, 1, [0.3]), (0, 1, [1.1]),
            (0, 1, [NaN]), (0, 1, [Inf]), (0, 1, Float64[]))
            @test_throws ErrorException smn(m, n, c, eta; precision, spheroid)
            @test_throws ErrorException accuracy(
                m, n, c, eta; precision, spheroid, target = :angular)
        end
    end
    for c in (NaN, Inf, complex(1.0, Inf)), target in (:angular, :radial)

        @test_throws ErrorException accuracy(0, 1, c, [0.3]; target)
    end
    for (m, n, c, x, kind) in ((-1, 0, 1.0, [2.0], 2), (2, 1, 1.0, [2.0], 2),
        (0, 1, 1.0, [NaN], 2), (0, 1, 1.0, Float64[], 2),
        (0, 1, 1.0, [0.5], 2), (0, 1, 1.0, [2.0], 5))
        @test_throws ErrorException rmn(m, n, c, x; kind)
        @test_throws ErrorException accuracy(m, n, c, x; kind)
    end
    @test_throws DomainError accuracy(0, 1, big"0", [big"2"]; precision = :quad)
    @test SpheroidalWaves._reported_accuracy([14, -1, 999, 12], [1.0, 1.0, 1.0, Inf], :double) ==
          [14, -1, -1, -1]
end

@testset "Accuracy tracks the requested native evaluation" begin
    for precision in (:double, :quad)
        lib = SpheroidalWaves.backend_library(; precision)
        if lib === nothing || !isfile(lib)
            @info "Skipping numerical accuracy tests: backend unavailable." precision
            continue
        end
        for spheroid in (:prolate, :oblate), c in (big"1.25", complex(big"1.25", big"0.2"))

            for kind in (1, 3, 4)
                @test accuracy(0, 1, c, [big"2"]; precision, spheroid, kind) == [-1]
            end
            angular = accuracy(0, 1, c, [big"0.3"]; precision, spheroid, target = :angular)
            radial = accuracy(0, 1, c, [big"2"]; precision, spheroid, kind = 2)
            @test length(angular) == length(radial) == 1
            @test all(a -> -1 <= a <= (precision === :quad ? 33 : 16), [angular; radial])
        end
    end
    # A quad radial coordinate just outside x=1 must not round onto the boundary
    # in the accuracy wrapper. Compare the estimate with the native evaluation
    # that also returns the function, ensuring both use the same inputs.
    lib = SpheroidalWaves.backend_library(; precision = :quad)
    if lib !== nothing && isfile(lib)
        x = [big"1"+big"1e-20"]
        r = SpheroidalWaves._call_real_rmn(:psms, 0, 1, big"1.25", x;
            precision = :quad, kind = 2, with_accuracy = true)
        @test all(isfinite, r.value)
        @test accuracy(0, 1, big"1.25", x; precision = :quad, kind = 2) == r.accuracy
        @test r.accuracy != [-1]
        eta = [big"1"-big"1e-20"]
        a = SpheroidalWaves._call_real_smn(:psms, 1, 1, big"1.25", eta;
            precision = :quad, with_accuracy = true)
        @test only(a.value) != 0
        @test accuracy(1, 1, big"1.25", eta; precision = :quad, target = :angular) ==
              a.accuracy
        old_library = SpheroidalWaves.backend_library(; precision = :double)
        if old_library !== nothing && isfile(old_library)
            try
                SpheroidalWaves.set_backend_library!(old_library; precision = :quad)
                for target in (:angular, :radial)
                    err = try
                        accuracy(0, 1, big"1.25", target === :angular ? eta : x;
                            precision = :quad, target, kind = (target === :angular ? 1 : 2))
                    catch e
                        e
                    end
                    @test err isa ErrorException
                    @test occursin("Rebuild the native backend", sprint(showerror, err))
                end
            finally
                SpheroidalWaves.set_backend_library!(lib; precision = :quad)
            end
        end
    end
end
