# Run explicitly with: julia --project scripts/validate_accuracy.jl
# Optional first argument selects a different generated TSV reference file.
# Broad parameter/cut validation stays outside the package's CI test suite.
using SpheroidalWaves, Test

function read_references(path)
    lines = filter(line -> !startswith(line, "#") && !isempty(line), readlines(path))
    fields = Tuple(Symbol.(split(first(lines), '\t')))
    return [NamedTuple{fields}(Tuple(split(line, '\t'))) for line in lines[2:end]]
end

function reference_pair(row)
    return (complex(parse(BigFloat, row.value_re), parse(BigFloat, row.value_im)),
        complex(parse(BigFloat, row.slope_re), parse(BigFloat, row.slope_im)))
end

function reference_key(row, family = row.family)
    (family, row.geometry, row.m, row.n, row.c_re, row.c_im, row.x_re, row.x_im)
end

function validate_references(rows)
    lookup = Dict(reference_key(row) => row for row in rows)
    groups = Dict{Tuple, Vector{eltype(rows)}}()
    for row in rows
        row.family == "radial2" && continue
        key = (row.family, row.geometry, row.m, row.n, row.c_re, row.c_im,
            iszero(parse(BigFloat, row.x_im)))
        push!(get!(groups, key, eltype(rows)[]), row)
    end
    @testset "Independent accuracy matrix" begin
        for key in sort!(collect(keys(groups))), precision in (:double, :quad)

            family, geometry, ms, ns, cr, ci, real_coordinate = key
            m, n = parse(Int, ms), parse(Int, ns)
            spheroid = Symbol(geometry)
            T = precision === :quad ? BigFloat : Float64
            tolerance = precision === :quad ? big"2e-26" : big"2e-12"
            c = iszero(parse(T, ci)) ? parse(T, cr) : complex(parse(T, cr), parse(T, ci))
            samples = groups[key]
            points = real_coordinate ? [parse(T, row.x_re) for row in samples] :
                     [complex(parse(T, row.x_re), parse(T, row.x_im)) for row in samples]
            radial = family == "radial1"
            options = radial ?
                      [(kind, normalization) for kind in 1:4
                       for normalization in (kind <= 2 ? (:standard, :static) :
                                             (:standard,))] :
                      family == "angular1" ? [(1, false), (1, true)] : [(2, false)]
            for (kind, normalization) in options
                @testset "$precision $geometry $family m=$m n=$n c=$c real_coordinate=$real_coordinate kind=$kind normalization=$normalization" begin
                    expected = map(samples) do row
                        pair = reference_pair(row)
                        if radial
                            second = reference_pair(lookup[reference_key(row, "radial2")])
                            pair = kind == 1 ? pair :
                                   kind == 2 ? second :
                                   pair .+ (kind == 3 ? im : -im) .* second
                            normalization == :static &&
                                (pair = pair .* Complex{BigFloat}(c)^(kind == 1 ? -n : n+1))
                        elseif normalization
                            norm2 = 2BigFloat(factorial(big(n+m))) /
                                    ((2n+1)*BigFloat(factorial(big(n-m))))
                            pair = pair ./ sqrt(norm2)
                        end
                        pair
                    end
                    result = try
                        radial ?
                        rmn(m, n, c, points; spheroid, precision, kind,
                            normalization, scaled = true) :
                        smn(m, n, c, points; spheroid, precision, kind,
                            normalize = normalization, scaled = true)
                    catch error
                        println("Evaluation error for ", key, ": ", sprint(showerror, error))
                        nothing
                    end
                    @test result !== nothing
                    if result !== nothing
                        errors = zeros(BigFloat, 2)
                        for (i, row) in enumerate(samples),
                            (j, field) in enumerate((:value, :derivative))

                            component = getproperty(result, field)
                            actual = component.mantissa[i]*BigFloat(10)^component.exponent[i]
                            # Relative comparisons, with an absolute floor for exact zeros.
                            error = abs(actual-expected[i][j])/max(abs(expected[i][j]), big"1e-50")
                            errors[j] = max(errors[j], error)
                            if !(error <= tolerance)
                                println((; precision, family, geometry, m, n, c, kind,
                                    normalization, x_re = row.x_re, x_im = row.x_im,
                                    field, actual, expected = expected[i][j], error))
                            end
                            @test error <= tolerance
                        end
                        println((;
                            precision, family, geometry, m, n, c, kind, normalization,
                            real_coordinate, value_error = errors[1], derivative_error = errors[2]))
                    end
                end
            end
            flush(stdout)
        end
    end
end

path = isempty(ARGS) ? joinpath(@__DIR__, "..", "test", "accuracy_references.tsv") :
       first(ARGS)
println((;
    julia = VERSION, package = pkgversion(SpheroidalWaves), references = abspath(path)))
setprecision(BigFloat, 256) do
    validate_references(read_references(path))
end
