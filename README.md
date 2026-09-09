# SpheroidalWaves.jl

[![Documentation](https://img.shields.io/badge/docs-latest-blue)](https://brandynlucca.github.io/SpheroidalWaves.jl)
[![CI](https://img.shields.io/github/actions/workflow/status/brandynlucca/SpheroidalWaves.jl/CI.yml?branch=main&label=CI)](https://github.com/brandynlucca/SpheroidalWaves.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19728040.svg)](https://doi.org/10.5281/zenodo.19728040)

Fast, vectorized computation of prolate and oblate spheroidal wave functions using native Fortran kernels.

The package supports real and complex spheroidal parameters, double and quad precision, angular and radial functions, eigenvalues, accuracy estimates, derivatives, and inverse eigenvalue calculations.

## Installation

```julia
import Pkg
Pkg.add("SpheroidalWaves")
```

Prebuilt backends are downloaded automatically on supported platforms. See [BUILD.md](BUILD.md) if a local backend build is required.

## Quick start

```julia
using SpheroidalWaves

# Angular function and derivative at several points
angular = smn(1, 2, 1.5, [-0.5, 0.0, 0.5])
angular.value
angular.derivative

# Radial function of the first kind
radial = rmn(1, 2, 1.5, [1.1, 1.2, 1.3]; kind=1)

# Oblate geometry
oblate = smn(1, 2, 1.5, [-0.5, 0.0, 0.5]; spheroid=:oblate)

# Complex spheroidal parameter
complex_result = smn(1, 2, 1.5 + 0.5im, [-0.5, 0.0, 0.5])

# Quad-precision backend
quad_result = rmn(1, 2, 1.5, [1.1, 1.2]; precision=:quad)

# Characteristic value
lambda = eigenvalue(1, 2, 1.5)
```

Pass arrays to `smn` and `rmn` to use the batch-oriented backend. Independent calls may also be scheduled concurrently with Julia threads.

## Main API

- `smn`: angular function and derivative
- `rmn`: radial function and derivative
- `eigenvalue` and `eigenvalue_sweep`: characteristic values
- `accuracy`: estimated numerical accuracy
- `radial_wronskian`: radial Wronskian
- `jacobian_eigen`, `jacobian_smn`, and `jacobian_rmn`: derivatives with respect to parameters
- `find_c_for_eigenvalue`: inverse characteristic-value calculation

The main options are `spheroid=:prolate` or `:oblate`, `precision=:double` or `:quad`, and radial `kind=1:4`.

See the [documentation](https://brandynlucca.github.io/SpheroidalWaves.jl) for definitions, argument restrictions, normalization, and complete examples.

## Attribution

The numerical kernels are based on the spheroidal-wave-function implementations by Arnie Lee Van Buren and Jeffrey Boisvert:

- [Prolate](https://github.com/MathieuandSpheroidalWaveFunctions/Prolate_swf)
- [Oblate](https://github.com/MathieuandSpheroidalWaveFunctions/Oblate_swf)
- [Complex prolate](https://github.com/MathieuandSpheroidalWaveFunctions/complex_prolate_swf)
- [Complex oblate](https://github.com/MathieuandSpheroidalWaveFunctions/complex_oblate_swf)

## Citation and license

Use the [Zenodo record](https://doi.org/10.5281/zenodo.19728040) to cite SpheroidalWaves.jl. The package is available under the [MIT License](LICENSE).
