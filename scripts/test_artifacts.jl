#!/usr/bin/env julia

# Run in the disposable artifact-build checkout: this replaces Artifacts.toml
# with bindings to the freshly packaged libraries before starting Julia tests.
using Pkg
using Pkg.Artifacts
using Pkg.PlatformEngines: unpack
using Tar

include("update_artifacts.jl")

function test_artifacts(dist::AbstractString, triplet::AbstractString)
    project = dirname(@__DIR__)
    artifacts_toml = joinpath(project, "Artifacts.toml")
    for precision in ("double", "quad")
        tarball = abspath(dist, "spheroidal_backend_$precision-$triplet.tar.gz")
        # Use Pkg's unpacking and read-only artifact installation, so Windows
        # permissions are exercised exactly where the published DLL failed.
        hash = create_artifact() do root
            unpack(tarball, root)
        end
        @assert string(hash) == tree_hash_file(tarball)
        bind_artifact!(artifacts_toml, "spheroidal_backend_$precision", hash; force=true)
    end

    # A fresh process must discover the artifacts rather than retain the local
    # build paths from an earlier module initialization. Overrides must not mask
    # a packaging error.
    code = """
        using Pkg, Artifacts, SpheroidalWaves
        for precision in (:double, :quad)
            name = "spheroidal_backend_\$precision"
            root = artifact_path(artifact_hash(name, $(repr(artifacts_toml))))
            filename = SpheroidalWaves._backend_filename("spheroidal_batch_\$precision")
            path = SpheroidalWaves.backend_library(; precision)
            @assert path == joinpath(root, filename)
            @info "Testing packaged backend" precision path
            result = SpheroidalWaves.smn(0, 0, 1.0, [0.25]; precision)
            @assert isfinite(only(result.value))
        end
        Pkg.test()
        """
    env = copy(ENV)
    pop!(env, "SPHEROIDALWAVES_LIBRARY_DOUBLE", nothing)
    pop!(env, "SPHEROIDALWAVES_LIBRARY_QUAD", nothing)
    run(setenv(`$(Base.julia_cmd()) --startup-file=no --project=$project -e $code`, env))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("Usage: julia scripts/test_artifacts.jl <artifact-dist> <triplet>")
    test_artifacts(ARGS...)
end
