#!/usr/bin/env julia
# deps/build.jl
# Julia package build script: compiles Fortran batch modules into shared library
# Runs automatically during package installation

using Artifacts
using Libdl

# ============================================================================
# Configuration
# ============================================================================
const SCRIPT_DIR = dirname(@__FILE__)
const PROJECT_DIR = dirname(SCRIPT_DIR)
const BUILD_DIR = joinpath(PROJECT_DIR, "build")
const BUILD_LOG = joinpath(SCRIPT_DIR, "build_output.txt")
const ARTIFACTS_TOML = joinpath(PROJECT_DIR, "Artifacts.toml")
const SELECTED_FORTRAN = Ref{Union{Nothing,String}}(nothing)

function detect_library_dir()
    candidates = [joinpath(BUILD_DIR, "lib"), joinpath(BUILD_DIR, "bin"), BUILD_DIR]
    for dir in candidates
        if isdir(dir)
            return dir
        end
    end
    return joinpath(BUILD_DIR, "lib")
end

# ============================================================================
# Helper Functions
# ============================================================================
function info_msg(msg)
    println("[SpheroidalWaves.jl] $msg")
end

function error_msg(msg)
    println("[SpheroidalWaves.jl ERROR] $msg")
end

function warn_msg(msg)
    println("[SpheroidalWaves.jl WARNING] $msg")
end

function _compiler_target(path::String)
    try
        return lowercase(strip(readchomp(Cmd([path, "-dumpmachine"]))))
    catch
        return ""
    end
end

function _target_matches_julia_arch(target::String)
    if isempty(target)
        return false
    end
    if Sys.ARCH == :x86_64
        return occursin("x86_64", target) || occursin("mingw64", target)
    elseif Sys.ARCH == :i686
        return occursin("i686", target) || occursin("mingw32", target)
    end
    return true
end

function _gfortran_candidates()
    candidates = String[]
    exe = Sys.iswindows() ? "gfortran.exe" : "gfortran"
    p = Sys.which(exe)
    if p !== nothing
        push!(candidates, p)
    end

    if Sys.iswindows()
        try
            lines = split(readchomp(`where gfortran`), '\n')
            for line in lines
                s = strip(replace(line, '\r' => ""))
                if !isempty(s)
                    push!(candidates, s)
                end
            end
        catch
            # ignore when where.exe is unavailable
        end
    end

    return unique(candidates)
end

function choose_fortran_compiler()
    gfortran_paths = _gfortran_candidates()
    for path in gfortran_paths
        target = _compiler_target(path)
        if _target_matches_julia_arch(target)
            return (path, "gfortran", target)
        end
    end
    if !isempty(gfortran_paths)
        path = first(gfortran_paths)
        return (path, "gfortran", _compiler_target(path))
    end

    for compiler in ["ifort", "ifx"]
        exe = Sys.iswindows() ? "$compiler.exe" : compiler
        path = Sys.which(exe)
        if path !== nothing
            return (path, compiler, "")
        end
    end

    return nothing
end

function choose_cmake_generator()
    if Sys.iswindows()
        return "MinGW Makefiles"
    end
    return nothing
end

# ============================================================================
# Step 1: Check for CMake
# ============================================================================
function check_cmake()
    cmake_exe = Sys.iswindows() ? "cmake.exe" : "cmake"
    cmake_path = Sys.which(cmake_exe)
    
    if cmake_path === nothing
        error_msg("CMake not found. Please install CMake 3.15 or later.")
        error_msg("  Windows: https://cmake.org/download/")
        error_msg("  macOS: brew install cmake")
        error_msg("  Linux: apt-get install cmake (or equivalent)")
        return false
    end
    
    info_msg("Found CMake: $cmake_path")
    return true
end

# ============================================================================
# Step 2: Check for Fortran Compiler
# ============================================================================
function check_fortran_compiler()
    chosen = choose_fortran_compiler()
    if chosen !== nothing
        path, compiler, target = chosen
        SELECTED_FORTRAN[] = path
        if !isempty(target)
            info_msg("Found Fortran compiler: $path (using $compiler, target=$target)")
        else
            info_msg("Found Fortran compiler: $path (using $compiler)")
        end
        return true
    end
    
    error_msg("No Fortran compiler found (gfortran, ifort, or ifx).")
    error_msg("  Windows: Install Intel Fortran or MinGW-w64 (gfortran)")
    error_msg("  macOS: brew install gcc")
    error_msg("  Linux: apt-get install gfortran (or equivalent)")
    return false
end

# ============================================================================
# Step 3: Verify Source Files Exist
# ============================================================================
function verify_sources()
    base_solvers_double = [
        "prolate_swf.f90",
        "oblate_swf.f90",
        "complex_prolate_swf.f90",
        "complex_oblate_swf.f90",
    ]
    
    batch_wrappers_double = [
        "psms_batch_fortran.f90",
        "oblate_batch_fortran.f90",
        "complex_prolate_batch_fortran.f90",
        "complex_oblate_batch_fortran.f90",
    ]

    base_solvers_quad = [
        "prolate_swf_quad.f90",
        "oblate_swf_quad.f90",
        "complex_prolate_swf_quad.f90",
        "complex_oblate_swf_quad.f90",
    ]

    batch_wrappers_quad = [
        "psms_batch_fortran_quad.f90",
        "oblate_batch_fortran_quad.f90",
        "complex_prolate_batch_fortran_quad.f90",
        "complex_oblate_batch_fortran_quad.f90",
    ]

    all_sources = [base_solvers_double; batch_wrappers_double; base_solvers_quad; batch_wrappers_quad]
    deps_dir = joinpath(PROJECT_DIR, "deps")

    for src in all_sources
        src_path = joinpath(deps_dir, src)
        if !isfile(src_path)
            error_msg("Source file not found: $src_path")
            return false
        end
    end
    
    info_msg("All $(length(all_sources)) source files verified.")
    return true
end

# ============================================================================
# Step 4: Run CMake Build
# ============================================================================
function run_build()
    # Create build directory
    mkpath(BUILD_DIR)
    
    # Run CMake configure
    info_msg("Running CMake configure...")
    configure_cmd = [
        "cmake",
        "-S", PROJECT_DIR,
        "-B", BUILD_DIR,
    ]

    generator = choose_cmake_generator()
    if generator !== nothing
        push!(configure_cmd, "-G", generator)
    end

    if SELECTED_FORTRAN[] !== nothing
        push!(configure_cmd, "-DCMAKE_Fortran_COMPILER=$(SELECTED_FORTRAN[])")
    end
    
    info_msg("Command: $(join(configure_cmd, " "))")
    
    try
        run(Cmd(configure_cmd))
    catch e
        error_msg("CMake configure failed: $e")
        return false
    end
    
    # Run CMake build
    info_msg("Running CMake build...")
    build_cmd = ["cmake", "--build", BUILD_DIR, "--config", "Release"]
    
    info_msg("Command: $(join(build_cmd, " "))")
    
    try
        open(BUILD_LOG, "w") do io
            run(pipeline(Cmd(build_cmd), stdout=io, stderr=io))
        end
    catch e
        error_msg("CMake build failed: $e")
        if isfile(BUILD_LOG)
            error_msg("Last build log lines:")
            for line in split(read(BUILD_LOG, String), '\n')[max(1, end - 39):end]
                if !isempty(line)
                    println(line)
                end
            end
            error_msg("Full build log: $BUILD_LOG")
        end
        return false
    end
    
    info_msg("Build completed successfully.")
    return true
end

# ============================================================================
# Step 5: Verify Library Exists
# ============================================================================
function _backend_filename(stem::String)
    if Sys.iswindows()
        return "$stem.dll"
    elseif Sys.isapple()
        return "lib$stem.dylib"
    else
        return "lib$stem.so"
    end
end

function _backend_filename_candidates(stem::String)
    names = String[]
    push!(names, _backend_filename(stem))
    if !Sys.iswindows()
        unprefixed = Sys.isapple() ? "$stem.dylib" : "$stem.so"
        if unprefixed ∉ names
            push!(names, unprefixed)
        end
    end
    return names
end

function _find_built_library(lib_dir::String, stem::String)
    for name in _backend_filename_candidates(stem)
        path = joinpath(lib_dir, name)
        if isfile(path)
            return path
        end
    end
    return nothing
end

function verify_libraries_built()
    lib_dir = detect_library_dir()
    stems = ["spheroidal_batch_double", "spheroidal_batch_quad"]
    ok = true
    for stem in stems
        path = _find_built_library(lib_dir, stem)
        if path === nothing
            error_msg("Library not found after build for $stem in $lib_dir")
            ok = false
        else
            info_msg("Library verified: $path")
        end
    end
    return ok
end

# ============================================================================
# Step 6: Report Library Paths
# ============================================================================
function report_library_paths()
    lib_dir = detect_library_dir()
    lib_double_path = _find_built_library(lib_dir, "spheroidal_batch_double")
    lib_quad_path = _find_built_library(lib_dir, "spheroidal_batch_quad")
    if lib_double_path === nothing || lib_quad_path === nothing
        error_msg("Cannot report library paths because one or more built libraries were not found.")
        return false
    end
    info_msg("double library path: $lib_double_path")
    info_msg("quad library path: $lib_quad_path")
    return true
end

function _find_library_in_root(root::String, stem::String)
    for dir in (root, joinpath(root, "lib"), joinpath(root, "bin"))
        path = _find_built_library(dir, stem)
        path === nothing || return path
    end
    return nothing
end

function _artifact_library(name::String, stem::String)
    isfile(ARTIFACTS_TOML) || return nothing
    hash = Artifacts.artifact_hash(name, ARTIFACTS_TOML)
    hash === nothing && return nothing
    Artifacts.artifact_exists(hash) || return nothing
    return _find_library_in_root(Artifacts.artifact_path(hash), stem)
end

function use_prebuilt_artifacts()
    double = _artifact_library("spheroidal_backend_double", "spheroidal_batch_double")
    quad = _artifact_library("spheroidal_backend_quad", "spheroidal_batch_quad")
    if double !== nothing && quad !== nothing
        info_msg("Using prebuilt backend artifacts for this platform.")
        info_msg("double artifact library: $double")
        info_msg("quad artifact library: $quad")
        info_msg("Local CMake/Fortran compilation is not required.")
        return true
    end
    return false
end

# ============================================================================
# Main Build Workflow
# ============================================================================
function main()
    info_msg("SpheroidalWaves Fortran batch build starting...")

    # Pkg installs non-lazy artifacts before running this build script. Avoid
    # compiling locally when both precision backends are already available.
    if use_prebuilt_artifacts()
        return true
    end

    info_msg("Building dual precision backends (double and quad)")
    info_msg("Build directory: $BUILD_DIR")
    
    # Check prerequisites
    if !check_cmake()
        return false
    end
    
    if !check_fortran_compiler()
        return false
    end
    
    if !verify_sources()
        return false
    end
    
    # Run build
    if !run_build()
        return false
    end
    
    # Verify output
    if !verify_libraries_built()
        return false
    end
    
    # Report the library locations. The module probes these known build paths.
    if !report_library_paths()
        return false
    end
    
    info_msg("Build successful!")
    info_msg("")
    info_msg("Next steps:")
    info_msg("  1. Library is ready at: $(detect_library_dir())")
    info_msg("  2. Use SpheroidalWaves module in Julia")
    info_msg("  3. Call smn/rmn with precision=:double or precision=:quad")
    info_msg("")
    
    return true
end

# ============================================================================
# Entry Point
# ============================================================================
if !main()
    is_ci = get(ENV, "CI", "false") == "true"
    
    if is_ci
        error("Build failed in CI. This indicates a real installation problem. See errors above.")
    else
        warn_msg("Local backend build failed. This is expected if:")
        warn_msg("  - You don't have a Fortran compiler installed")
        warn_msg("  - CMake is not available")
        warn_msg("")
        warn_msg("You can still use SpheroidalWaves if:")
        warn_msg("  1. Pre-built artifacts are available (automatic download)")
        warn_msg("  2. You set environment variables with backend paths")
        warn_msg("  3. You install a Fortran compiler and rebuild")
        warn_msg("")
        warn_msg("To install a Fortran compiler:")
        warn_msg("  - Ubuntu/Debian: sudo apt-get install gfortran cmake")
        warn_msg("  - macOS: brew install gcc cmake")
        warn_msg("  - Windows: install MinGW-w64 (gfortran) and CMake")
        warn_msg("")
        warn_msg("Then rebuild with: julia> import Pkg; Pkg.build(\"SpheroidalWaves\")")
    end
end

