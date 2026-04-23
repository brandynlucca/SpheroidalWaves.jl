# SpheroidalWaveFunctions Build System

This document describes the build process for compiling the Fortran batch kernels into a shared library for Julia use.

## Architecture Overview

```
deps/
  ├── prolate_swf.f90              # Base prolate solver (Fortran)
  ├── oblate_swf.f90               # Base oblate solver (Fortran)
  ├── complex_prolate_swf.f90       # Base complex prolate solver (Fortran)
  ├── complex_oblate_swf.f90        # Base complex oblate solver (Fortran)
  ├── psms_batch_fortran.f90        # Batch wrapper for prolate (iso_c_binding)
  ├── oblate_batch_fortran.f90      # Batch wrapper for oblate (iso_c_binding)
  ├── complex_prolate_batch_fortran.f90    # Batch wrapper for complex prolate
  ├── complex_oblate_batch_fortran.f90     # Batch wrapper for complex oblate
  ├── build.jl                      # Julia build script (runs during pkg install)
  └── library_config.jl             # Generated at build time (library path config)

src/
  └── SpheroidalWaveFunctions.jl     # Julia FFI layer (ccall wrappers)

CMakeLists.txt                      # Cross-platform build configuration

build/                              # Build artifacts (generated)
  ├── lib/spheroidal_batch.so       # Shared library (Unix/Linux)
  ├── lib/spheroidal_batch.dll      # Shared library (Windows)
  └── lib/libspheroidal_batch.dylib # Shared library (macOS)
```

## Public API

The package exposes **two entry points** that intelligently dispatch to the appropriate Fortran solver:

```julia
using SpheroidalWaveFunctions

# Angular functions
vals = smn(m, n, c, η; spheroid=:prolate, precision=:double, normalize=false)

# Radial functions
vals = rmn(m, n, c, x; spheroid=:prolate, precision=:double, kind=1)
```

**Dispatch Logic:**
- If `c` is `Real`: uses double-precision real solver
- If `c` is `Complex`: uses complex-valued solver
- `spheroid=:prolate` (default) or `spheroid=:oblate` selects geometry
- `precision=:double` (default) or `precision=:quad` selects precision at runtime
- `kind` parameter for radial functions selects Bessel-like (1–2) or Hankel-like (3–4) types

**Example Usage:**
```julia
# Real prolate angular function, double precision
vals = smn(1, 2, 1.5, [0.5, 0.6, 0.7])
println(vals.value)        # Array of function values
println(vals.derivative)   # Array of derivatives

# Complex prolate radial function, quad precision
vals = rmn(1, 2, 1.5 + 0.5im, [1.1, 1.2, 1.3]; precision=:quad)
println(typeof(vals.value[1]))  # ComplexF64

# Oblate geometry
vals = smn(1, 2, 1.5, [0.5, 0.6]; spheroid=:oblate)
```

## Precision (knd) Selection

### Runtime Precision Selection

Precision is selected at runtime via the `precision` keyword argument to `smn()` and `rmn()`:

```julia
# Double precision (knd=8) — default
vals = smn(1, 2, 1.5, [0.5, 0.6])

# Quad precision (knd=16) — if supported by library
vals = smn(1, 2, 1.5, [0.5, 0.6]; precision=:quad)
```

### Compile-Time Precision Setup

The current build system compiles with **double precision (knd=8)** by default. To support quad precision at runtime, the library must be built to include both precision variants.

**Future enhancement** (in development): Build system will compile separate batches of base solvers for knd=8 and knd=16, making both available simultaneously in the shared library. Then Julia's runtime dispatch via `precision=:double` or `precision=:quad` will select the appropriate entry points.

## Build Prerequisites

### 1. CMake
- **Version**: 3.15 or later
- **Installation**:
  - **Windows**: Download from [cmake.org](https://cmake.org/download/)
  - **macOS**: `brew install cmake`
  - **Linux**: `apt-get install cmake` (or equivalent)

### 2. Fortran Compiler
- **Supported**: gfortran, Intel Fortran (ifort/ifx)
- **Installation**:
  - **Windows (gfortran)**: Install MinGW-w64 with gfortran component
  - **Windows (Intel)**: Intel oneAPI Base Toolkit + Fortran Compiler
  - **macOS**: `brew install gcc`
  - **Linux**: `apt-get install gfortran` (or equivalent)

### 3. Julia 1.10+
- SpheroidalWaveFunctions requires Julia 1.10 or later

## Build Methods

### Method 1: Automatic Build (Recommended)

When you install the package using Julia's package manager, the build script runs automatically:

```julia
julia> import Pkg; Pkg.add("SpheroidalWaveFunctions")
```

Or if you have a local development copy:

```julia
julia> import Pkg; Pkg.develop("path/to/SpheroidalWaveFunctions")
```

The build process will:
1. Detect your system (Windows/Linux/macOS)
2. Find an available Fortran compiler
3. Run CMake to configure the build
4. Compile all 8 Fortran source files
5. Create a shared library (`.so`, `.dll`, or `.dylib`)
6. Register the library path for Julia

**Expected output**:
```
[SpheroidalWaveFunctions.jl] SpheroidalWaveFunctions Fortran batch build starting...
[SpheroidalWaveFunctions.jl] Found CMake: /usr/bin/cmake
[SpheroidalWaveFunctions.jl] Found Fortran compiler: /usr/bin/gfortran (using gfortran)
[SpheroidalWaveFunctions.jl] All 8 source files verified.
[SpheroidalWaveFunctions.jl] Running CMake configure...
[SpheroidalWaveFunctions.jl] Running CMake build...
[SpheroidalWaveFunctions.jl] Library verified: .../build/lib/libspheroidal_batch.so
[SpheroidalWaveFunctions.jl] Build successful!
```

### Method 2: Manual CMake Build

If you want to build manually or customize compiler flags:

```bash
# Navigate to package root
cd /path/to/SpheroidalWaveFunctions

# Create build directory
mkdir -p build && cd build

# Configure (default: double precision, knd=8)
cmake .. -DFORTRAN_PRECISION=8

# Build
cmake --build . --config Release

# Library output: build/lib/spheroidal_batch.so (or .dll/.dylib)
```

### Method 3: Rebuild After Modification

If you modify Fortran source files and need to rebuild:

```julia
julia> import Pkg; Pkg.build("SpheroidalWaveFunctions")
```

Or:

```bash
cd /path/to/SpheroidalWaveFunctions
rm -rf build
# Then re-run Method 2
```

## Precision Selection

The build system supports two precision backends:

- **Double precision (knd=8)**: Default, 64-bit floating point
- **Quad precision (knd=16)**: Extended precision, 128-bit floating point

### Building with Different Precision

**Using Julia package build**:
```julia
julia> ENV["FORTRAN_PRECISION"] = "16"
julia> import Pkg; Pkg.build("SpheroidalWaveFunctions")
```

**Using CMake directly**:
```bash
cmake .. -DFORTRAN_PRECISION=16
```

**Note**: Quad precision requires Fortran compiler support (gfortran or Intel Fortran with `-freal=16` or `-real_size 64` flags). Not all systems support quad precision.

## Troubleshooting

### "CMake not found"
- Install CMake from [cmake.org](https://cmake.org/download/)
- Or use package manager: `brew install cmake` (macOS), `apt-get install cmake` (Linux)

### "No Fortran compiler found"
- Install gfortran: 
  - Windows: MinGW-w64 (https://www.mingw-w64.org/)
  - macOS: `brew install gcc`
  - Linux: `apt-get install gfortran`
- Or install Intel Fortran from [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)

### "Source file not found"
- Ensure all `.f90` files exist in `deps/` directory
- Do not move or rename source files
- Run from the package root directory

### "Build failed with compilation error"
- Check error messages above; common issues:
  - Missing `-fPIC` flag (fixed in CMakeLists.txt)
  - Fortran standard compatibility (set to F2008)
  - Module load order (handled by CMake dependency order)
- Try with verbose output:
  ```bash
  cmake --build . --config Release --verbose
  ```

## Troubleshooting

### "CMake not found"
- Install CMake from [cmake.org](https://cmake.org/download/)
- Or use package manager: `brew install cmake` (macOS), `apt-get install cmake` (Linux)

### "No Fortran compiler found"
- Install gfortran: 
  - Windows: MinGW-w64 (https://www.mingw-w64.org/)
  - macOS: `brew install gcc`
  - Linux: `apt-get install gfortran`
- Or install Intel Fortran from [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)

### "Source file not found"
- Ensure all `.f90` files exist in `deps/` directory
- Do not move or rename source files
- Run from the package root directory

### "Build failed with compilation error"
- Check error messages above; common issues:
  - Missing `-fPIC` flag (fixed in CMakeLists.txt)
  - Fortran standard compatibility (set to F2008)
  - Module load order (handled by CMake dependency order)
- Try with verbose output:
  ```bash
  cmake --build . --config Release --verbose
  ```

### "Library loads but function call fails"
- Verify library path:
  ```julia
  julia> using SpheroidalWaveFunctions
  julia> SpheroidalWaveFunctions.backend_library()
  ```
- Should show path like `/path/to/build/lib/libspheroidal_batch.so`
- If empty, rebuild: `Pkg.build("SpheroidalWaveFunctions")`
- Or manually set:
  ```julia
  julia> SpheroidalWaveFunctions.set_backend_library!("/path/to/build/lib/libspheroidal_batch.so")
  ```

### "Error calling Fortran function"
- Verify status code in error message
- Test with known good inputs: `smn(1, 2, 1.0, [1.0])`
- If still failing, library may be for wrong precision (e.g., calling double-precision data with quad library)

## Library Path Discovery

After build, Julia discovers the library via:

1. **Automatic (via `library_config.jl`)**:
   - File created by build script at `deps/library_config.jl`
   - Loaded at module initialization (`__init__()` in `src/SpheroidalWaveFunctions.jl`)
   - Path automatically set via `set_backend_library!(path)`

2. **Manual override**:
   ```julia
   julia> using SpheroidalWaveFunctions
   julia> set_backend_library!("/custom/path/to/libspheroidal_batch.so")
   julia> prolate_smn_batch(1, 2, 1.0, [1.0])  # Now uses custom library
   ```

3. **Environment variable** (not yet implemented, future):
   ```bash
   export SPHEROIDAL_BATCH_LIBRARY="/path/to/libspheroidal_batch.so"
   julia
   ```

## Build Artifacts

### Output Files

After a successful build:

```
build/
├── lib/
│   ├── spheroidal_batch.so       (Linux/Unix)
│   ├── spheroidal_batch.dll      (Windows)
│   └── libspheroidal_batch.dylib (macOS)
├── CMakeFiles/
├── CMakeCache.txt
├── SpheroidalWaveFunctions.cmake
└── cmake_install.cmake
```

### Library Size

- Typical library size: 1–3 MB (depending on optimization, compiler)
- Contains 8 public entry points (batch functions for 4 solver types)
- Statically linked Fortran runtime (no additional dependencies)

## Performance Notes

### Compilation Time
- First build: 10–30 seconds (depends on system, compiler)
- Rebuild (cache): 2–5 seconds
- Quad precision: typically slightly slower due to extended-precision math

### Runtime Performance
- Batch operations: O(n) where n is number of evaluation points
- No per-call overhead (Fortran directly called via ccall)
- Characteristic-exponent reconstruction adds minimal cost

## Development Workflow

### For Package Maintainers

1. **Modify Fortran source**:
   - Edit files in `deps/`
   - Ensure iso_c_binding declarations stay consistent

2. **Test build**:
   ```bash
   cd /path/to/SpheroidalWaveFunctions
   rm -rf build
   mkdir build && cd build
   cmake ..
   cmake --build .
   ```

3. **Update version**:
   - Edit `Project.toml` version
   - Update `CHANGELOG` with build system changes

4. **Release**:
   - Ensure build passes on all platforms
   - Tag release: `git tag v1.0.0`
   - Consider pre-compiled binary distribution (future: BinaryBuilder.jl)

### For Users

1. **Install**:
   ```julia
   julia> import Pkg; Pkg.add("SpheroidalWaveFunctions")
   ```

2. **Use**:
   ```julia
   julia> using SpheroidalWaveFunctions
   julia> vals = prolate_smn_batch(1, 2, 1.5, [0.5, 0.6, 0.7])
   julia> vals.value  # Batch results
   ```

3. **Rebuild if needed**:
   ```julia
   julia> Pkg.build("SpheroidalWaveFunctions")
   ```

## Cross-Platform Compatibility

The build system is designed to work on:

- **Windows** (10, 11): Visual Studio, MinGW, Intel Fortran
- **macOS** (12+): Apple Clang + gfortran/Intel Fortran
- **Linux** (glibc 2.29+): GCC/gfortran, Intel Fortran

All platforms produce compatible `.so`/`.dll`/`.dylib` files with same ABI.

## Future Improvements

- [ ] Pre-compiled binaries via BinaryBuilder.jl
- [ ] Support for quad precision on all platforms
- [ ] Automatic precision selection based on system capabilities
- [ ] Optional GPU acceleration (CUDA/OpenCL)
- [ ] Environment variable for library path: `SPHEROIDAL_BATCH_LIBRARY`

## Questions or Issues?

See [README.md](README.md) for contact information and issue tracking.
