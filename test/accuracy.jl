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
        wide = Complex{BigFloat}.(points)
        reference = target===:angular ? atanh.(wide) :
                    spheroid===:prolate ? -atanh.(inv.(wide)) :
                    atan.(wide) .- BigFloat(pi)/2
        reference_derivative = target===:angular ? inv.(1 .- wide .^ 2) :
                               inv.(wide .^ 2 .- (spheroid===:prolate ? 1 : -1))
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
            @test diagnostic.value.status === diagnostic.derivative.status === :estimated
            @test abs(actual.value[i]-reference[i]) <= diagnostic.value.absolute_error
            @test abs(actual.derivative[i]-reference_derivative[i]) <=
                  diagnostic.derivative.absolute_error
            @test diagnostic.working_precision[2] > diagnostic.working_precision[1]
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
    @test only(zero_report.diagnostics).value.status === :estimated
    @test isnan(only(zero_report.diagnostics).value.relative_error)
    endpoints = accuracy(
        1, 1, 1.25, [1+0im, 0.25+0.125im]; target = :angular, diagnostics = true)
    @test endpoints.digits[1] == -1 && endpoints.digits[2]>=12
    @test endpoints.diagnostics[1].conditioning_flag === :unavailable
    @test endpoints.diagnostics[1].value.status === :unavailable
    @test endpoints.diagnostics[1].derivative.status === :singularity
    real_report = accuracy(0, 0, 0, [0.25]; target = :angular, diagnostics = true)
    @test real_report.digits == [-1] &&
          only(real_report.diagnostics).value.status === :estimated
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
              :nonfinite
        sensitive = Complex{F}.((1, 10000, 0, 0))
        warning = SW._coordinate_accuracy_record(
            sensitive, sensitive, plan, 1, :double, false)
        @test warning.diagnostics.conditioning_flag === :warning
        @test warning.digits < 12
    end
end

@testset "Separate value and derivative accuracy estimates" begin
    function check_estimates(actual, diagnostic, value, derivative)
        for (field, reference) in ((:value, value), (:derivative, derivative))
            estimate = getproperty(diagnostic, field)
            @test estimate.status === :estimated
            @test abs(only(getproperty(actual, field))-reference) <= estimate.absolute_error
            @test estimate.absolute_error isa BigFloat
            @test !iszero(reference) && estimate.relative_error >= 0
        end
        @test diagnostic.working_precision[2] > diagnostic.working_precision[1]
        @test diagnostic.terms[2] > diagnostic.terms[1]
    end
    # Independent c=1 spectral references bracketing a zero and a stationary
    # point. These exact dyadic inputs are also used by numerical_regressions.jl.
    zeros = (
        (:prolate, 2, 2402525102,
            big"-3.646498962489905788454236914910610385883e-10",
            big"1.703323695364498173128344421120509893166775201302"),
        (:oblate, 3, 1952290389,
            big"-0.45130974952979024172775506093481452827491145456698",
            big"-0.00000000088576585426959176119068231968610374136907004850296"))
    for precision in (:double, :quad), (spheroid, n, k, value, derivative) in zeros

        normalize = spheroid === :oblate
        point = k/2.0^32
        report = accuracy(0, n, 1, [point]; spheroid, precision, normalize,
            target = :angular, diagnostics = true)
        actual = smn(0, n, 1, point; spheroid, precision, normalize)
        factor = normalize ? sqrt(BigFloat(2n+1)/2) : BigFloat(1)
        check_estimates(actual, only(report.diagnostics), factor*value, factor*derivative)
        @test report.digits ==
              accuracy(0, n, 1, [point]; spheroid, precision, normalize, target = :angular)
    end
    # Large-bandwidth cancellation uses the independently converged fixture
    # from oblate_precision.jl, with the same 256-bit input construction.
    point = BigFloat(3)/10
    for precision in (:double, :quad)
        report = accuracy(0, 0, 50, [point]; spheroid = :oblate, precision,
            target = :angular, diagnostics = true)
        actual = smn(0, 0, 50, point; spheroid = :oblate, precision)
        check_estimates(actual, only(report.diagnostics),
            big"9.67772408518644447309773518088656173136537178202329502e-15",
            big"4.763828624750672125456540294407534653167498239290307627e-13")
        @test report.digits == [-1]
    end
    # Read only four interior pairs from the independent matrix. The broad
    # cut survey and its Python generator are never executed by package tests.
    rows = [split(line, '\t')
            for line in readlines(joinpath(@__DIR__, "accuracy_references.tsv"))
            if startswith(line, "radial")]
    function radial_reference(spheroid, ci, kind)
        pair = map(("radial1", "radial2")) do family
            row = only(filter(
                r -> r[1] == family && r[2] == string(spheroid) &&
                     r[3:9] == ["0", "0", "1", ci, "2", "0", "interior"],
                rows))
            (complex(parse(BigFloat, row[10]), parse(BigFloat, row[11])),
                complex(parse(BigFloat, row[12]), parse(BigFloat, row[13])))
        end
        kind <= 2 ? pair[kind] : pair[1] .+ (kind == 3 ? im : -im) .* pair[2]
    end
    for precision in (:double, :quad), kind in 1:4

        spheroid = isodd(kind) ? :prolate : :oblate
        c = kind <= 2 ? 1 : 1+0.25im
        reference = radial_reference(spheroid, kind <= 2 ? "0" : "0.25", kind)
        report = accuracy(0, 0, c, [2]; spheroid, precision, kind, diagnostics = true)
        actual = rmn(0, 0, c, 2; spheroid, precision, kind)
        check_estimates(actual, only(report.diagnostics), reference...)
        @test only(report.diagnostics).method === :radial_refinement
    end
    # The ordinary double output underflows. Its absolute error must remain
    # nonzero, using the independent decaying-wave fixture in radial_hankel.jl.
    report = accuracy(2, 4, 20im, [40]; spheroid = :oblate, kind = 3, diagnostics = true)
    diagnostic = only(report.diagnostics)
    @test diagnostic.value.status === diagnostic.derivative.status === :estimated
    @test diagnostic.value.digits == diagnostic.derivative.digits == 0
    @test diagnostic.value.absolute_error >=
          big"3.802038016847395155618029612498744199455358714591274911e-351"
    @test diagnostic.derivative.absolute_error >=
          big"7.611800331549399620260330117854328694240686039155645507e-350"
end

@testset "Accuracy diagnostic failure states" begin
    unavailable = accuracy(0, 0, 0, [0.25]; target = :angular, kind = 2, diagnostics = true)
    @test only(unavailable.diagnostics).value.status === :unavailable
    @test isnan(only(unavailable.diagnostics).value.absolute_error)
    static = accuracy(0, 1, 0, [2]; normalization = :static, diagnostics = true)
    @test only(static.diagnostics).derivative.status === :unavailable
    regular = accuracy(0, 1, 0, [0, 1]; target = :angular, diagnostics = true)
    @test regular.diagnostics[1].value.status === :estimated
    @test regular.diagnostics[1].value.digits == -1
    @test isnan(regular.diagnostics[1].value.relative_error)
    @test isfinite(regular.diagnostics[1].value.absolute_error)
    @test regular.diagnostics[2].value.status === :unavailable
    singular = accuracy(1, 1, 0, [1]; target = :angular, diagnostics = true)
    @test only(singular.diagnostics).value.status === :unavailable
    @test only(singular.diagnostics).derivative.status === :singularity
    for (target, point) in ((:angular, 1), (:radial, 1))
        report = accuracy(0, 0, 1, [point]; target, kind = 2, diagnostics = true)
        @test only(report.diagnostics).value.status === :singularity
        @test only(report.diagnostics).derivative.status === :singularity
    end
    overflow = accuracy(200, 200, 0, [0.3]; target = :angular, diagnostics = true)
    @test only(overflow.diagnostics).value.status === :nonfinite
    SW = SpheroidalWaves
    failed = SW._real_accuracy_diagnostics(0, 0, 0, [0.25], :prolate, :double,
        :angular, 1, false, :standard;
        reference = (args...; kwargs...) -> error("Unresolved reference"))
    @test only(failed).value.status === only(failed).derivative.status ===
          :refinement_failed
    @test only(failed).value.digits == -1
    @test_throws ArgumentError SW._accuracy_refinement(() -> throw(ArgumentError("Invalid option")))
    SW._with_swprecision(384) do
        F = SW._SWFloat
        failed = SW._accuracy_component(F(1), F(1), F(2), F(2), :double)
        @test failed.status === :refinement_failed && failed.digits == -1

        # Failed continuation refinement reuses the original state. Identical
        # fallback states must not imply convergence or a zero error estimate.
        point = Complex{F}(0.25, 0.125)
        plan = SW._sensitivity_plan(0, 1, F(0), :prolate, :double)
        state = Complex{F}.((point, 1, 0, 0))
        record = SW._coordinate_accuracy_record(
            state, state, plan, point, :double, false; refinement_failed = true)
        @test record.digits == -1
        @test !record.diagnostics.converged
        @test record.diagnostics.finite_flag
        @test record.diagnostics.conditioning_flag === :poor
        for component in (record.diagnostics.value, record.diagnostics.derivative)
            @test component.status === :refinement_failed
            @test component.digits == -1
            @test isnan(component.absolute_error)
            @test isnan(component.relative_error)
        end
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
