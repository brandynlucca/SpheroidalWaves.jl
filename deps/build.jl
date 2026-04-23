#!/usr/bin/env julia
# deps/build.jl
# Julia package build script: compiles Fortran batch modules into shared library
# Runs automatically during package installation

using Pkg
using Libdl

# ============================================================================
# Configuration
# ============================================================================
const SCRIPT_DIR = dirname(@__FILE__)
const PROJECT_DIR = dirname(SCRIPT_DIR)
const BUILD_DIR = joinpath(PROJECT_DIR, "build")
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
    println("[SpheroidalWaveFunctions.jl] $msg")
end

function error_msg(msg)
    println("[SpheroidalWaveFunctions.jl ERROR] $msg")
end

function warn_msg(msg)
    println("[SpheroidalWaveFunctions.jl WARNING] $msg")
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
    base_solvers_r8 = [
        "prolate_swf.f90",
        "oblate_swf.f90",
        "complex_prolate_swf.f90",
        "complex_oblate_swf.f90",
    ]
    
    batch_wrappers_r8 = [
        "psms_batch_fortran.f90",
        "oblate_batch_fortran.f90",
        "complex_prolate_batch_fortran.f90",
        "complex_oblate_batch_fortran.f90",
    ]

    base_solvers_r16 = [
        "prolate_swf_r16.f90",
        "oblate_swf_r16.f90",
        "complex_prolate_swf_r16.f90",
        "complex_oblate_swf_r16.f90",
    ]

    batch_wrappers_r16 = [
        "psms_batch_fortran_r16.f90",
        "oblate_batch_fortran_r16.f90",
        "complex_prolate_batch_fortran_r16.f90",
        "complex_oblate_batch_fortran_r16.f90",
    ]

    all_sources = [base_solvers_r8; batch_wrappers_r8; base_solvers_r16; batch_wrappers_r16]
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

    if SELECTED_FORTRAN[] !== nothing
        push!(configure_cmd, "-DCMAKE_Fortran_COMPILER=$(SELECTED_FORTRAN[])")
    end
    
    # Add platform-specific flags
    if Sys.iswindows()
        push!(configure_cmd, "-G", "MinGW Makefiles")  # or "Unix Makefiles" if using Unix-like on Windows
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
        run(Cmd(build_cmd))
    catch e
        error_msg("CMake build failed: $e")
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

function verify_libraries_built()
    lib_dir = detect_library_dir()
    expected = [_backend_filename("spheroidal_batch_r8"), _backend_filename("spheroidal_batch_r16")]
    ok = true
    for name in expected
        path = joinpath(lib_dir, name)
        if !isfile(path)
            error_msg("Library not found after build: $path")
            ok = false
        else
            info_msg("Library verified: $path")
        end
    end
    return ok
end

# ============================================================================
# Step 6: Create Library Path Configuration
# ============================================================================
function configure_library_paths()
    lib_dir = detect_library_dir()
    lib_r8 = replace(joinpath(lib_dir, _backend_filename("spheroidal_batch_r8")), '\\' => '/')
    lib_r16 = replace(joinpath(lib_dir, _backend_filename("spheroidal_batch_r16")), '\\' => '/')

    config_script = """
    # Auto-generated: library path configuration
    # Set on package load in src/SpheroidalWaveFunctions.jl

    const SPHEROIDAL_BATCH_LIBRARY_R8 = "$lib_r8"
    const SPHEROIDAL_BATCH_LIBRARY_R16 = "$lib_r16"
    """

    config_file = joinpath(SCRIPT_DIR, "library_config.jl")
    write(config_file, config_script)

    info_msg("Library path configuration written to: $config_file")
    info_msg("r8 library path: $lib_r8")
    info_msg("r16 library path: $lib_r16")
    return true
end

# ============================================================================
# Main Build Workflow
# ============================================================================
function main()
    info_msg("SpheroidalWaveFunctions Fortran batch build starting...")
    info_msg("Building dual precision backends (r8 and r16)")
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
    
    # Configure Julia to find library
    if !configure_library_paths()
        return false
    end
    
    info_msg("Build successful!")
    info_msg("")
    info_msg("Next steps:")
    info_msg("  1. Library is ready at: $(detect_library_dir())")
    info_msg("  2. Use SpheroidalWaveFunctions module in Julia")
    info_msg("  3. Call smn/rmn with precision=:double or precision=:quad")
    info_msg("")
    
    return true
end

# ============================================================================
# Entry Point
# ============================================================================
if !main()
    error("Build failed. See errors above.")
end
