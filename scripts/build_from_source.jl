#!/usr/bin/env julia
# scripts/build_from_source.jl
# Optional source build; normal installations use prebuilt SpheroidalWaves_jll libraries.

using Libdl

# ============================================================================
# Configuration
# ============================================================================
const SCRIPT_DIR = dirname(@__FILE__)
const PROJECT_DIR = dirname(SCRIPT_DIR)
const BUILD_DIR = abspath(get(ENV, "SWF_BUILD_DIR", joinpath(PROJECT_DIR, "build")))
const BUILD_TYPE = get(ENV, "SWF_BUILD_TYPE", "Release")
const BUILD_LOG = joinpath(BUILD_DIR, "build_output.txt")
const SELECTED_FORTRAN = Ref{Union{Nothing, String}}(nothing)
const BUILD_QUAD = lowercase(get(ENV, "SWF_BUILD_QUAD", "on")) ∉ ("off", "false", "0")
const BUILD_MPFR = lowercase(get(ENV, "SWF_BUILD_MPFR", "on")) ∉ ("off", "false", "0")
const BUILD_TESTS = lowercase(get(ENV, "SWF_BUILD_TESTS", "off")) ∉ ("off", "false", "0")
const RUN_TESTS = BUILD_TESTS && lowercase(get(ENV, "SWF_RUN_TESTS", "on")) ∉ ("off", "false", "0")
const BACKENDS = ["double"; BUILD_QUAD ? ["quad"] : []; BUILD_MPFR ? ["mpfr"] : []]

# ============================================================================
# Helper Functions
# ============================================================================
function info_msg(msg)
    println("[SpheroidalWaves.jl] $msg")
end

function error_msg(msg)
    println("[SpheroidalWaves.jl ERROR] $msg")
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
    if haskey(ENV, "FC")
        path = Sys.which(ENV["FC"])
        path === nothing && error("Fortran compiler specified by FC was not found: $(ENV["FC"])")
        return (path, basename(path), _compiler_target(path))
    end
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
    haskey(ENV, "CMAKE_GENERATOR") && return ENV["CMAKE_GENERATOR"]
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
    all_sources = ["precision.inc", "precision.F90", "scaled_batch.f90"]
    for geometry in ("prolate", "oblate", "complex_prolate", "complex_oblate")
        push!(all_sources, "$(geometry)_kernel.F90")
    end
    for wrapper in ("psms", "oblate", "complex_prolate", "complex_oblate")
        push!(all_sources, "$(wrapper)_batch_fortran.f90")
        BUILD_QUAD && push!(all_sources, "$(wrapper)_batch_fortran_quad.f90")
    end
    BUILD_MPFR && append!(all_sources,
        ["mpfr_real.f90", "mpfr_complex.f90", "mpfr_batch.f90",
         "mpfr_arithmetic.h", "mpfr_arithmetic.c", "mpc_arithmetic.c"])
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
        "-DCMAKE_BUILD_TYPE=$BUILD_TYPE",
        "-DSWF_BUILD_QUAD=$(BUILD_QUAD ? "ON" : "OFF")",
        "-DSWF_BUILD_MPFR=$(BUILD_MPFR ? "ON" : "OFF")",
        "-DSWF_BUILD_TESTS=$(BUILD_TESTS ? "ON" : "OFF")"
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
    build_cmd = ["cmake", "--build", BUILD_DIR, "--config", BUILD_TYPE, "--parallel", "1"]

    info_msg("Command: $(join(build_cmd, " "))")

    try
        open(BUILD_LOG, "w") do io
            output = Pipe()
            process = run(pipeline(Cmd(build_cmd), stdout = output, stderr = output); wait = false)
            close(output.in)
            try
                for line in eachline(output)
                    println(io, line)
                    flush(io)
                    println(line)
                end
            finally
                close(output)
            end
            wait(process)
            success(process) || error("CMake build exited with code $(process.exitcode)")
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
# Step 5: Verify Libraries and Entry Points
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

function _find_built_library(stem::String)
    directories = [joinpath(BUILD_DIR, subdir, config)
                   for subdir in ("lib", "bin", "") for config in (BUILD_TYPE, "")]
    for dir in directories, name in _backend_filename_candidates(stem)
        path = joinpath(dir, name)
        if isfile(path)
            return path
        end
    end
    return nothing
end

function verify_libraries_built()
    symbols = Dict(
        "double" => (:psms_eigenvalue_double, :oblate_eigenvalue_double,
                     :cprolate_eigenvalue_c8, :coblate_eigenvalue_c8, :spheroidal_scaled_text),
        "quad" => (:psms_eigenvalue_quad_fullsplit, :oblate_eigenvalue_quad_fullsplit,
                   :cprolate_batch_quad_text, :coblate_batch_quad_text, :spheroidal_scaled_text),
        "mpfr" => (:spheroidal_mpfr_batch,))
    for backend in BACKENDS
        path = _find_built_library("spheroidal_batch_$backend")
        if path === nothing
            error_msg("Library not found after build for $backend in $BUILD_DIR")
            return false
        end
        try
            Libdl.dlopen(path) do handle
                for symbol in symbols[backend]
                    Libdl.dlsym(handle, symbol)
                end
            end
        catch e
            error_msg("Cannot load $backend library or its required symbols: $e")
            return false
        end
        info_msg("$backend library verified: $path")
        info_msg("SPHEROIDALWAVES_LIBRARY_$(uppercase(backend))=$path")
    end
    return true
end

# ============================================================================
# Main Build Workflow
# ============================================================================
function main()
    info_msg("SpheroidalWaves native backend build starting...")
    info_msg("Building $(join(BACKENDS, ", ")) backends ($BUILD_TYPE)")
    BUILD_MPFR && info_msg("MPFR/MPC requires C/Fortran compilers and MPFR/MPC/GMP development files.")
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

    if RUN_TESTS
        try
            run(Cmd(`ctest --build-config $BUILD_TYPE -j1 --output-on-failure`; dir=BUILD_DIR))
        catch e
            error_msg("Native backend tests failed: $e")
            return false
        end
    end

    info_msg("Build successful!")

    return true
end

# ============================================================================
# Entry Point
# ============================================================================
if !main()
    error_msg("Source build failed; see errors above and $BUILD_LOG.")
    info_msg("Normal installations use prebuilt SpheroidalWaves_jll libraries and do not run this script.")
    exit(1)
end
