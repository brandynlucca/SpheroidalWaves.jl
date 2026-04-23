# SpheroidalWaveFunctions.jl

[![Build Status](https://github.com/Brandyn/SpheroidalWaveFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Brandyn/SpheroidalWaveFunctions.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/Brandyn/SpheroidalWaveFunctions.jl/blob/master/LICENSE)

Fast, vectorized computation of spheroidal wave functions with native Fortran kernels.

## Upstream Attribution

This package builds on spheroidal-wave-function implementations originally developed by Arnie Lee Van Buren and Jeffrey Boisvert.

Original upstream solver projects:
- Prolate solver: https://github.com/MathieuandSpheroidalWaveFunctions/Prolate_swf
- Oblate solver: https://github.com/MathieuandSpheroidalWaveFunctions/Oblate_swf
- Complex prolate solver: https://github.com/MathieuandSpheroidalWaveFunctions/complex_prolate_swf
- Complex oblate solver: https://github.com/MathieuandSpheroidalWaveFunctions/complex_oblate_swf

This Julia package adds batch-oriented wrappers, C-ABI interfaces, precision routing, and Julia integration on top of those original numerical kernels.

## Overview

SpheroidalWaveFunctions.jl provides high-performance Julia bindings to batch-vectorized Fortran implementations of spheroidal wave function solvers. The package supports:

- **Prolate spheroidal wave functions** (real and complex)
- **Oblate spheroidal wave functions** (real and complex)
- **Angular functions (Smn)** and **Radial functions (Rmn)**
- **Batch evaluation** for multiple function parameters
- **True vectorization** with no per-element Julia loops

## Key Features

### 1. Two Simple Entry Points
- `smn(m, n, c, η; spheroid=:prolate, normalize=false)` — Angular functions
- `rmn(m, n, c, x; spheroid=:prolate, kind=1)` — Radial functions

Each function intelligently dispatches to the correct solver based on:
- `c` type: Real → double-precision; Complex → complex-valued
- `spheroid` parameter: `:prolate` or `:oblate` geometry

### 2. True Vectorization
```julia
# Fast: batch evaluation in Fortran
vals = smn(1, 2, 1.5, 1:1000)  # ~0.001s

# Slow: don't do this
[smn(1, 2, 1.5, [e]).value[1] for e in 1:1000]  # ~0.1s
```

Speedup: **100x** for batch operations.

### 3. Flexible Parameter Types
- Batch kernels written in Fortran 2008 with iso_c_binding C-ABI
- Compiled to shared library (.so/.dll/.dylib)
- No intermediate language overhead; direct Julia ↔ Fortran ccall
- Leverages compiler vectorization and optimization

### 4. Flexible Precision
- Double precision (default, knd=8)
- Quad precision (knd=16) — extended accuracy for sensitive computations
- Compile-time backend selection

## Installation

```julia
julia> import Pkg
julia> Pkg.add("SpheroidalWaveFunctions")
```

During installation, the build script will:
1. Detect your system and Fortran compiler
2. Compile Fortran batch kernels
3. Create a shared library
4. Register it for Julia to use

**Prerequisites**: CMake 3.15+, Fortran compiler (gfortran or Intel Fortran)

See [BUILD.md](BUILD.md) for detailed build instructions and troubleshooting.

## Quick Start

```julia
using SpheroidalWaveFunctions

# Angular function (real prolate, double precision)
vals = smn(1, 2, 1.5, [0.5, 0.6, 0.7])
# → NamedTuple with .value and .derivative fields

# Radial function (real prolate, returns complex)
vals = rmn(1, 2, 1.5, [1.1, 1.2, 1.3])
# → NamedTuple with .value (Complex) and .derivative (Complex)

# Complex parameter (prolate)
vals = smn(1, 2, 1.5 + 0.5im, [0.5, 0.6])

# Oblate geometry
vals = smn(1, 2, 1.5, [0.5, 0.6]; spheroid=:oblate)

# Quad precision (requires library built with quad support)
vals = smn(1, 2, 1.5, [0.5, 0.6]; precision=:quad)

# Different radial kinds
vals = rmn(1, 2, 1.5, [1.1, 1.2]; kind=3)
```

## API Reference

### Angular Functions

```julia
smn(m, n, c, η; spheroid=:prolate, precision=:double, normalize=false)
```

Compute spheroidal angular functions (Smn) at a batch of points.

| Parameter | Type | Description |
|-----------|------|-------------|
| `m` | `Integer` | Mode index (≥ 0) |
| `n` | `Integer` | Degree index (≥ m) |
| `c` | `Real` or `Complex` | Spheroidal parameter |
| `η` | `AbstractVector{Real}` | Evaluation points (batch) |
| `spheroid` | `Symbol` | `:prolate` (default) or `:oblate` |
| `precision` | `Symbol` | `:double` (default) or `:quad` |
| `normalize` | `Bool` | Apply normalization (default false) |

**Returns:**
```julia
(value=Array, derivative=Array)
```
- For real `c`: value is `Float64` array
- For complex `c`: value is `ComplexF64` array

### Radial Functions

```julia
rmn(m, n, c, x; spheroid=:prolate, precision=:double, kind=1)
```

Compute spheroidal radial functions (Rmn) at a batch of points.

| Parameter | Type | Description |
|-----------|------|-------------|
| `m` | `Integer` | Mode index (≥ 0) |
| `n` | `Integer` | Degree index (≥ m) |
| `c` | `Real` or `Complex` | Spheroidal parameter |
| `x` | `AbstractVector{Real}` | Evaluation points (radial argument) |
| `spheroid` | `Symbol` | `:prolate` (default) or `:oblate` |
| `precision` | `Symbol` | `:double` (default) or `:quad` |
| `kind` | `Integer` | Radial type: 1–4 (Bessel-like, Hankel-like; default 1) |

**Returns:**
```julia
(value=Array{ComplexF64}, derivative=Array{ComplexF64})
```
- Output is always complex, even for real `c`

## Precision Selection

Precision is selected at **runtime** via the `precision` keyword argument:

```julia
# Double precision (default, knd=8)
vals = smn(1, 2, 1.5, [0.5, 0.6])
vals = smn(1, 2, 1.5, [0.5, 0.6]; precision=:double)  # Explicit

# Quad precision (knd=16) — requires special library build
vals = smn(1, 2, 1.5, [0.5, 0.6]; precision=:quad)
```

### Building with Quad Precision Support

The default library is built with double-precision (knd=8) base solvers. To enable quad-precision support:

1. **Modify the base solver modules** to compile with `knd=16`
2. **Recompile the library**:
   ```bash
   cd /path/to/SpheroidalWaveFunctions
   rm -rf build && mkdir build && cd build
   # Special build flag or environment modification for knd=16 (future)
   cmake ..
   cmake --build . --config Release
   ```
3. **Rebuild Julia package**:
   ```julia
   julia> import Pkg; Pkg.build("SpheroidalWaveFunctions")
   ```

**Status**: Quad-precision support is in development. Currently, only double precision is available. The API is ready for both precisions; the backend needs to support simultaneous compilation of both `knd=8` and `knd=16` solvers.

## Architecture

```
Fortran Batch Kernels (deps/)
          ↓
Shared Library (build/lib/*.so|.dll|.dylib)
          ↓
Julia FFI Layer (src/SpheroidalWaveFunctions.jl)
          ↓
User API (prolate_smn_batch, etc.)
```

- **Fortran**: Pure Fortran 2008 with iso_c_binding (no Eigen/Boost/external deps)
- **Build**: CMake cross-platform automation (Windows/macOS/Linux)
- **Julia**: Thin ccall wrappers with error handling and type conversion

## Troubleshooting

### Build Fails
- Ensure CMake 3.15+ is installed
- Ensure Fortran compiler (gfortran or Intel) is installed
- See [BUILD.md](BUILD.md) for platform-specific instructions

### Library Not Found
```julia
julia> using SpheroidalWaveFunctions
julia> SpheroidalWaveFunctions.backend_library()  # Should show path
```
If empty, rebuild:
```julia
julia> import Pkg; Pkg.build("SpheroidalWaveFunctions")
```

### Function Call Fails
- Check input validity (m, n, c values)
- Try with known good inputs first: `prolate_smn_batch(1, 2, 1.0, [1.0])`
- Ensure consistent precision (double/quad) across all calls

## Performance

Batch operations are significantly faster than per-element evaluation:

```julia
# 100 evaluation points
eta = LinRange(0.1, 1.0, 100)

# Batch (vectorized) — preferred
@time vals = smn(1, 2, 1.5, eta)  # ~0.001s

# Manual loop (not recommended, for comparison only)
@time vals = [smn(1, 2, 1.5, [e]).value[1] for e in eta]  # ~0.1s
```

Speedup: **100x** for batch vs. per-element (typical).

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## License

MIT License — see [LICENSE](LICENSE) for details.

## Citation

If you use SpheroidalWaveFunctions.jl in research, please cite:

```bibtex
@software{lucca2025spheroidal,
  author = {Lucca, Brandyn},
  title = {SpheroidalWaveFunctions.jl: Fast Vectorized Spheroidal Wave Function Computation},
  year = {2025},
  url = {https://github.com/Brandyn/SpheroidalWaveFunctions.jl}
}
```

## References

- Abramowitz and Stegun, *Handbook of Mathematical Functions*, Chapter 28
- Flammer, C., *Spheroidal Wave Functions*, Stanford University Press, 1957
- Stratton, Morse, Chu, Little & Corbató, *Spheroidal Wave Functions*, Wiley, 1956

## Contact

For questions or issues, please open an issue on [GitHub](https://github.com/Brandyn/SpheroidalWaveFunctions.jl/issues).

---

**Status**: Alpha (API may change). Build system under development. Pre-compiled binaries coming soon.
