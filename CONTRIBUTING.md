# Contributing to SpheroidalWaves.jl

Contributions can address numerical correctness, performance, documentation, or new spheroidal-function capabilities. Use [GitHub issues](https://github.com/brandynlucca/SpheroidalWaves.jl/issues) for bug reports, feature proposals, and usage questions. Discuss changes to public APIs, mathematical conventions, or dependencies in an issue before substantial implementation.

## Reporting bugs and requesting features

Search existing issues before opening a report. For numerical bugs, include:

- A minimal Julia example with exact inputs, keyword arguments, and input types.
- Expected and observed results, including error messages or numerical discrepancies.
- Geometry, mode indices, coordinates, normalization, and requested precision.
- Julia and package versions, operating system, and backend selection. Identify any source-built libraries or environment overrides.
- A reference calculation, mathematical identity, or citation supporting the expected result, when available.

For feature requests, describe the scientific use case, proposed behavior, and relevant mathematical conventions. For usage questions, include a minimal example and the intended calculation.

## Development setup

Use Julia 1.10 or later, consistent with [Project.toml](Project.toml). Start from a fork or an existing checkout and create a branch for the contribution. Run commands from the repository root.

Start Julia with the package environment:

```sh
julia --project=. --threads=2
```

Install dependencies and run the Julia tests:

```julia
using Pkg
Pkg.instantiate()
ENV["SPHEROIDALWAVES_TEST_GROUP"] = "all"
Pkg.test(julia_args=["--threads=2", "-O1"])
```

The package uses prebuilt libraries from `SpheroidalWaves_jll` on supported platforms. Full numerical validation requires double, quadruple, and arbitrary precision backends. An installed JLL release without the arbitrary precision backend requires a compatible JLL release or a source build.

## Building and testing compiled backends

For Fortran or C changes, build the backends from the working tree and run both compiled and Julia tests. Install CMake 3.15 or later, a Fortran compiler with OpenMP and `REAL(16)` support, a C compiler, and development headers and libraries for MPFR, MPC, and GMP. The [CI workflow](.github/workflows/CI.yml) contains platform setup examples.

Set `FC` to select a Fortran compiler and `CMAKE_PREFIX_PATH` when dependency discovery needs an installation prefix. On Windows, provide a compatible MinGW toolchain, `mingw32-make`, and required runtime libraries on `PATH`. The build script selects `MinGW Makefiles` on Windows and accepts a `CMAKE_GENERATOR` override.

Start a fresh Julia session from the repository root with `julia --project=. --threads=2`. Build all three backends and run the compiled regression tests:

```julia
ENV["SWF_BUILD_TESTS"] = "ON"
ENV["SWF_RUN_TESTS"] = "ON"
include("scripts/build_from_source.jl")
```

The default output directory is `build/`. Keep `SWF_BUILD_QUAD` and `SWF_BUILD_MPFR` enabled for full validation. Package initialization prefers JLL libraries over local build discovery, so select the compiled libraries before running Julia tests. In the same Julia session, using the default build location:

```julia
using Pkg, Libdl

library_dir = abspath("build", Sys.iswindows() ? "bin" : "lib")
library_prefix = Sys.iswindows() ? "" : "lib"
for backend in ("double", "quad", "mpfr")
    filename = "$(library_prefix)spheroidal_batch_$(backend).$(Libdl.dlext)"
    ENV["SPHEROIDALWAVES_LIBRARY_" * uppercase(backend)] = joinpath(library_dir, filename)
end
ENV["SPHEROIDALWAVES_REQUIRE_BACKENDS"] = "true"
ENV["SPHEROIDALWAVES_TEST_GROUP"] = "all"
Pkg.test(julia_args=["--threads=2", "-O1"])
```

For a custom build directory or configuration subdirectory, use the library paths printed by the build script. The backend requirement check verifies the selected library paths and required symbols before numerical tests.

## Numerical and code conventions

Preserve public signatures, return shapes, defaults, and mathematical conventions unless an agreed API change requires a revision. Follow the surrounding Julia, Fortran, and C style. Keep changes focused and avoid unrelated formatting.

- Add regression coverage for numerical fixes using an independent reference or mathematical identity. Record reference provenance and justify tolerances. Use absolute error near zeros and relative error away from zeros.
- Cover affected geometries, real or complex inputs, precision modes, normalization choices, and endpoint or branch behavior. Include precision-refinement checks for accuracy claims.
- Preserve thread safety and caller arithmetic settings. Check concurrent results against serial references for changes affecting caches, precision handling, or shared state.
- Keep test fixtures with the package tests and make tests self-contained. Avoid dependencies on external research directories or generated local outputs.

Document public arguments, return values, mathematical definitions, units or dimensionless conventions, and practical limitations. Preserve attribution for adapted algorithms and numerical references. Keep generated binaries, documentation output, and scratch files out of commits.

## Building the documentation

Start a separate Julia session from the repository root:

```sh
julia --project=docs
```

Use the working package checkout and build the documentation:

```julia
using Pkg
Pkg.develop(path=".")
Pkg.instantiate()
include("docs/make.jl")
```

Open `docs/build/index.html` for local review. Documenter checks documentation coverage for exported functions. Keep public examples aligned with API behavior and mathematical conventions.

## Submitting a pull request

Open a pull request against `main` and link any related issue. Describe the problem, resulting behavior, and numerical or API implications. Include validation commands and outcomes, plus any unavailable backends or untested platforms.

Run the full Julia suite for code changes, the compiled tests for backend changes, and the documentation build for documentation or public API changes. Review the diff for unrelated edits and generated files. Address failing continuous integration checks before requesting final review.

Contribute code under the repository [GPL-3.0 license](LICENSE) and retain applicable copyright and attribution notices. Logo artwork has a [separate license](docs/src/assets/LICENSE).
