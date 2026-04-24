#!/usr/bin/env julia

using SHA
using Tar
using TOML

function usage()
    println("Usage:")
    println("  julia scripts/update_artifacts.jl <Artifacts.toml> <r8_url> <r8_tarball> <r16_url> <r16_tarball>")
    println("  julia scripts/update_artifacts.jl <Artifacts.toml> <r8_url> <r8_sha256> <r8_git_tree_sha1> <r16_url> <r16_sha256> <r16_git_tree_sha1>")
    println("")
    println("Example:")
    println("  julia scripts/update_artifacts.jl Artifacts.toml https://.../r8.tar.gz r8.tar.gz https://.../r16.tar.gz r16.tar.gz")
    println("  julia scripts/update_artifacts.jl Artifacts.toml https://.../r8.tar.gz <sha256> <tree> https://.../r16.tar.gz <sha256> <tree>")
end

function set_artifact!(data::Dict{String,Any}, name::String, url::String, sha256::String, tree::String)
    data[name] = Dict(
        "git-tree-sha1" => tree,
        "download" => [Dict("url" => url, "sha256" => sha256)],
    )
end

function sha256_file(path::String)
    open(path, "r") do io
        return bytes2hex(sha256(io))
    end
end

function tree_hash_tarball(path::String)
    hash = Tar.tree_hash(path)
    return string(hash)
end

function main(args)
    if !(length(args) in (5, 7))
        usage()
        error("Expected 5 or 7 arguments, got $(length(args)).")
    end

    artifacts_path = args[1]
    r8_url = args[2]
    r16_url = args[length(args) == 5 ? 4 : 5]

    if length(args) == 5
        r8_tarball = args[3]
        r16_tarball = args[5]
        if !isfile(r8_tarball)
            error("r8 tarball not found: $r8_tarball")
        end
        if !isfile(r16_tarball)
            error("r16 tarball not found: $r16_tarball")
        end
        r8_sha = sha256_file(r8_tarball)
        r8_tree = tree_hash_tarball(r8_tarball)
        r16_sha = sha256_file(r16_tarball)
        r16_tree = tree_hash_tarball(r16_tarball)
    else
        r8_sha = args[3]
        r8_tree = args[4]
        r16_sha = args[6]
        r16_tree = args[7]
    end

    data = if isfile(artifacts_path) && !isempty(strip(read(artifacts_path, String)))
        TOML.parsefile(artifacts_path)
    else
        Dict{String,Any}()
    end

    set_artifact!(data, "spheroidal_backend_r8", r8_url, r8_sha, r8_tree)
    set_artifact!(data, "spheroidal_backend_r16", r16_url, r16_sha, r16_tree)

    open(artifacts_path, "w") do io
        TOML.print(io, data)
    end

    println("Updated $(artifacts_path) with artifact bindings for r8 and r16.")
    println("r8: sha256=$r8_sha tree=$r8_tree")
    println("r16: sha256=$r16_sha tree=$r16_tree")
end

main(ARGS)
