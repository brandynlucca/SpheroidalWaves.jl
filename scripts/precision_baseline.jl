#!/usr/bin/env julia
# Capture native outputs at fixed inputs; invoke with an output TSV path.
using Libdl, SHA, SpheroidalWaves

ENV["SPHEROIDALWAVES_REQUIRE_BACKENDS"] = "true"
include(joinpath(@__DIR__, "..", "test", "require_backends.jl"))

const BASELINE_WIDTH = 256
const BASELINE_BITS = 113

function baseline_text(values; terminated = false)
    result = fill(terminated ? UInt8(0) : UInt8(' '), BASELINE_WIDTH, length(values))
    for (column, value) in enumerate(values)
        text = string(value)
        ncodeunits(text) < BASELINE_WIDTH || error("Baseline text exceeds its buffer")
        result[1:ncodeunits(text), column] .= codeunits(text)
    end
    result
end

function baseline_component(buffer, channel; terminated = false)
    bytes = buffer[:, channel]
    if terminated
        stop = findfirst(iszero, bytes)
        stop === nothing && error("Unterminated native MPFR output")
        stop == 1 && return BigFloat(0)
        bytes = bytes[1:stop-1]
    end
    parse(BigFloat, strip(String(bytes)))
end

function baseline_pointer(backend, symbol)
    SpheroidalWaves._symbol_pointer(SpheroidalWaves.backend_library(; precision = backend), symbol)
end

function baseline_mpfr(prolate, complex_input, m, n, c, coordinate, mode)
    parameter = baseline_text((c, complex_input ? "0.25" : "0"); terminated = true)
    points = baseline_text((coordinate,); terminated = true)
    output = zeros(UInt8, BASELINE_WIDTH, 10)
    accuracy = zeros(Cint, 1)
    status = Ref{Cint}(0)
    pointer = baseline_pointer(:mpfr, :spheroidal_mpfr_batch)
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
            Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        prolate, complex_input, m, n, mode, 0, BASELINE_BITS, 0, 1, BASELINE_WIDTH,
        parameter, points, output, accuracy, status)
    channels = mode == 0 ? (1:1) : mode == 3 ? (2:5) : (2:3)
    values = status[] == 0 ? [complex(baseline_component(output, 2i-1; terminated = true),
        baseline_component(output, 2i; terminated = true)) for i in channels] : Complex{BigFloat}[]
    (; status = Int(status[]), accuracy = Int(accuracy[1]), values)
end

function baseline_native_eigenvalue(backend, prolate, complex_input, m, n, c)
    prefix = complex_input ? (prolate ? :cprolate : :coblate) : (prolate ? :psms : :oblate)
    status = Ref{Cint}(0)
    re = Ref{Cdouble}(0)
    im = Ref{Cdouble}(0)
    if backend === :quad && complex_input
        parameter = baseline_text((c, "0.25"))
        points = baseline_text(("0",))
        output = fill(UInt8(' '), BASELINE_WIDTH, 8)
        exponents = zeros(Cint, 8)
        accuracy = zeros(Cint, 1)
        pointer = baseline_pointer(backend, Symbol(prefix, "_batch_quad_text"))
        ccall(pointer, Cvoid,
            (Cint, Cint, Cint, Cint, Cint, Ptr{UInt8}, Ptr{UInt8}, Cint,
                Ptr{UInt8}, Ptr{Cint}, Ptr{Cint}, Ref{Cint}),
            m, n, 0, 0, 1, parameter, points, BASELINE_WIDTH, output, exponents, accuracy, status)
        value = status[] == 0 ? complex(baseline_component(output, 1), baseline_component(output, 2)) : complex(BigFloat(NaN))
    elseif backend === :quad
        pointer = baseline_pointer(backend, Symbol(prefix, "_eigenvalue_quad_fullsplit"))
        ccall(pointer, Cvoid,
            (Cint, Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}),
            m, n, parse(Float64, c), 0, re, im, status)
        value = complex(BigFloat(re[]) + BigFloat(im[]))
    elseif complex_input
        pointer = baseline_pointer(backend, Symbol(prefix, "_eigenvalue_c8"))
        ccall(pointer, Cvoid,
            (Cint, Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}),
            m, n, parse(Float64, c), 0.25, re, im, status)
        value = complex(BigFloat(re[]), BigFloat(im[]))
    else
        pointer = baseline_pointer(backend, Symbol(prefix, "_eigenvalue_double"))
        ccall(pointer, Cvoid, (Cint, Cint, Cdouble, Ref{Cdouble}, Ref{Cint}),
            m, n, parse(Float64, c), re, status)
        value = complex(BigFloat(re[]))
    end
    (; status = Int(status[]), accuracy = -1,
        values = status[] == 0 ? [value] : Complex{BigFloat}[])
end

function baseline_native(backend, prolate, complex_input, m, n, c, coordinate, mode)
    mode == 0 && return baseline_native_eigenvalue(backend, prolate, complex_input, m, n, c)
    parameter = baseline_text((c, complex_input ? "0.25" : "0"))
    x = parse(BigFloat, coordinate)
    mode >= 2 && prolate && (x -= 1)
    points = baseline_text((x,))
    output = fill(UInt8(' '), BASELINE_WIDTH, 8)
    exponents = zeros(Cint, 8)
    status = Ref{Cint}(0)
    pointer = baseline_pointer(backend, :spheroidal_scaled_text)
    ccall(pointer, Cvoid,
        (Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
            Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ref{Cint}),
        !prolate, complex_input, mode == 1 ? 1 : 2, mode == 1 ? 0 : mode-1,
        m, n, 1, BASELINE_WIDTH, parameter, points, output, exponents, status)
    channels = mode == 3 ? (1:4) : (1:2)
    values = status[] == 0 ? [complex(
        baseline_component(output, 2i-1)*BigFloat(10)^exponents[2i-1],
        baseline_component(output, 2i)*BigFloat(10)^exponents[2i]) for i in channels] : Complex{BigFloat}[]
    (; status = Int(status[]), accuracy = -1, values)
end

function baseline_snapshot(io)
    fixtures = ((m = 0, n = 0, c = "1", x = "2"),
        (m = 1, n = 2, c = "5", x = "1.25"),
        (m = 3, n = 6, c = "20", x = "1.25"))
    successes = Dict(backend => 0 for backend in (:double, :quad, :mpfr))
    rejections = 0
    differences = Dict{Tuple{Symbol, Int}, BigFloat}()
    println(io, "family\tfixture\tmode\tbackend\tstatus\taccuracy\tchannel\treal\timag")
    for prolate in (true, false), complex_input in (false, true), (index, fixture) in enumerate(fixtures), mode in 0:3
        family = (complex_input ? "complex_" : "real_") * (prolate ? "prolate" : "oblate")
        coordinate = mode == 0 ? "0" : mode == 1 ? "0.25" : fixture.x
        results = Dict{Symbol, Any}()
        for backend in (:double, :quad, :mpfr)
            evaluate() = backend === :mpfr ? baseline_mpfr(prolate, complex_input,
                fixture.m, fixture.n, fixture.c, coordinate, mode) :
                baseline_native(backend, prolate, complex_input,
                    fixture.m, fixture.n, fixture.c, coordinate, mode)
            result = evaluate()
            repeated = evaluate()
            isequal(result, repeated) || error("Repeated native call changed: $family/$index/$mode/$backend")
            results[backend] = result
            result.status == 0 || (backend === :mpfr && result.status == 4) ||
                error("Unexpected status $(result.status): $family/$index/$mode/$backend")
            index == 1 && result.status != 0 && error("Required small native case rejected: $family/$mode/$backend")
            if result.status == 0
                successes[backend] += 1
                all(isfinite, result.values) || error("Nonfinite native result: $family/$index/$mode/$backend")
                for (channel, value) in enumerate(result.values)
                    println(io, join((family, index, mode, backend, result.status, result.accuracy,
                        channel, real(value), imag(value)), '\t'))
                end
            else
                rejections += 1
                println(io, join((family, index, mode, backend, result.status, result.accuracy, 0, "", ""), '\t'))
            end
            flush(io)
        end
        reference = results[:quad].values
        for backend in (:double, :mpfr)
            result = results[backend]
            result.status == 0 || continue
            difference = maximum(abs(a-b)/max(abs(b), BigFloat("1e-100")) for (a, b) in zip(result.values, reference))
            key = (backend, mode)
            differences[key] = max(get(differences, key, BigFloat(0)), difference)
        end
    end
    println(stderr, "Native cases per backend: ", successes, "; MPFR workspace rejections: ", rejections,
        "; each case executed twice; MPFR bits: ", BASELINE_BITS)
    for key in sort!(collect(keys(differences)); by = x -> (string(x[1]), x[2]))
        println(stderr, key[1], " mode ", key[2], " max relative difference from quad: ", differences[key])
    end
end

function main()
    length(ARGS) == 1 || error("Usage: julia --project=. scripts/precision_baseline.jl OUTPUT.tsv")
    for backend in (:double, :quad, :mpfr)
        path = SpheroidalWaves.backend_library(; precision = backend)
        println(stderr, backend, " ", realpath(path), " sha256=", bytes2hex(open(sha256, path)))
    end
    mkpath(dirname(abspath(ARGS[1])))
    setprecision(BigFloat, 256) do
        open(ARGS[1], "w") do io
            baseline_snapshot(io)
        end
    end
end

main()
