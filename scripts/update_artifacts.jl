#!/usr/bin/env julia

using SHA
using Tar
using TOML

const ARTIFACT_PREFIX = "spheroidal_backend_"
const SUPPORTED_TRIPLETS = Dict(
    "x86_64-linux-gnu" => Dict("arch" => "x86_64", "os" => "linux", "libc" => "glibc"),
    "x86_64-w64-mingw32" => Dict("arch" => "x86_64", "os" => "windows"),
    "x86_64-apple-darwin" => Dict("arch" => "x86_64", "os" => "macos"),
    "aarch64-linux-gnu" => Dict("arch" => "aarch64", "os" => "linux", "libc" => "glibc"),
    "aarch64-apple-darwin" => Dict("arch" => "aarch64", "os" => "macos"),
)

function usage(io::IO=stdout)
    println(io, "Usage:")
    println(io, "  julia scripts/update_artifacts.jl <Artifacts.toml> <repository> <release-tag> <tarball>...")
    println(io)
    println(io, "Tarballs must be named:")
    println(io, "  spheroidal_backend_<double|quad>-<platform-triplet>.tar.gz")
end

sha256_file(path::AbstractString) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function tree_hash_file(path::AbstractString)
    gzip = Sys.which(Sys.iswindows() ? "gzip.exe" : "gzip")
    if gzip === nothing && Sys.iswindows()
        candidates = String[
            joinpath(get(ENV, "ProgramFiles", ""), "Git", "usr", "bin", "gzip.exe"),
            joinpath(get(ENV, "ProgramFiles(x86)", ""), "Git", "usr", "bin", "gzip.exe"),
            raw"C:\msys64\usr\bin\gzip.exe",
        ]
        gzip = findfirst(isfile, candidates)
        gzip = gzip === nothing ? nothing : candidates[gzip]
    end
    gzip === nothing && error(
        "gzip is required to compute the artifact tree hash for compressed tarballs. " *
        "Install gzip or add it to PATH.",
    )
    return string(Tar.tree_hash(Cmd([gzip, "-dc", abspath(path)])))
end

function platform_from_tarball(path::AbstractString)
    filename = basename(path)
    matched = match(r"^spheroidal_backend_(double|quad)-(.+)\.tar\.gz$", filename)
    matched === nothing && error("Unexpected artifact filename: $filename")

    precision, triplet = matched.captures
    platform = get(SUPPORTED_TRIPLETS, triplet, nothing)
    platform === nothing && error("Unsupported platform triplet in $filename: $triplet")
    return precision, triplet, platform
end

function release_url(repository::AbstractString, tag::AbstractString, path::AbstractString)
    return "https://github.com/$repository/releases/download/$tag/$(basename(path))"
end

function artifact_entry(
    path::AbstractString,
    repository::AbstractString,
    tag::AbstractString,
    platform::Dict{String,String},
)
    entry = Dict{String,Any}(platform)
    entry["git-tree-sha1"] = tree_hash_file(path)
    entry["download"] = [Dict(
        "url" => release_url(repository, tag, path),
        "sha256" => sha256_file(path),
    )]
    return entry
end

function main(args)
    if length(args) < 4
        usage(stderr)
        error("Expected at least 4 arguments, got $(length(args)).")
    end

    artifacts_path, repository, tag = args[1:3]
    tarballs = args[4:end]
    all(isfile, tarballs) || error("One or more artifact tarballs do not exist.")

    entries = Dict(
        "double" => Dict{String,Dict{String,Any}}(),
        "quad" => Dict{String,Dict{String,Any}}(),
    )

    for tarball in tarballs
        precision, triplet, platform = platform_from_tarball(tarball)
        haskey(entries[precision], triplet) && error("Duplicate $precision artifact for $triplet")
        entries[precision][triplet] = artifact_entry(tarball, repository, tag, platform)
    end

    double_platforms = Set(keys(entries["double"]))
    quad_platforms = Set(keys(entries["quad"]))
    double_platforms == quad_platforms || error(
        "Artifact platform mismatch: double=$(sort!(collect(double_platforms))), " *
        "quad=$(sort!(collect(quad_platforms)))",
    )
    isempty(double_platforms) && error("No complete double/quad artifact pairs were supplied.")

    data = if isfile(artifacts_path) && !isempty(strip(read(artifacts_path, String)))
        TOML.parsefile(artifacts_path)
    else
        Dict{String,Any}()
    end

    for precision in ("double", "quad")
        name = ARTIFACT_PREFIX * precision
        data[name] = [entries[precision][triplet] for triplet in sort!(collect(keys(entries[precision])))]
    end

    open(artifacts_path, "w") do io
        TOML.print(io, data; sorted=true)
    end

    println("Updated $artifacts_path for release $tag:")
    for triplet in sort!(collect(double_platforms))
        println("  $triplet (double and quad)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
