<h1 align="center">
  <img src="docs/src/assets/logo.svg" alt="SpheroidalWaves.jl" width="300">
</h1>

[![Julia Registry](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fraw.githubusercontent.com%2FJuliaRegistries%2FGeneral%2Fmaster%2FS%2FSpheroidalWaves%2FVersions.toml&search=.*%5C%5B%22%28%5B%5E%22%5D%2B%29%22%5C%5D&flags=s&replace=v%241&label=Julia%20Registry&color=blue)](https://platform.juliahub.com/ui/Packages/General/SpheroidalWaves)
[![GitHub Version](https://img.shields.io/github/v/release/brandynlucca/SpheroidalWaves.jl?label=GitHub)](https://github.com/brandynlucca/SpheroidalWaves.jl)
[![Yggdrasil](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fraw.githubusercontent.com%2FJuliaPackaging%2FYggdrasil%2Fmaster%2FS%2FSpheroidalWaves%2Fbuild_tarballs.jl&search=version%5Cs%2A%3D%5Cs%2Av%22%28%5B%5E%22%5D%2B%29%22&replace=v%241&label=Yggdrasil&color=forestgreen)](https://platform.juliahub.com/ui/Packages/General/SpheroidalWaves_jll)
[![Julia Compatibility](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fbrandynlucca%2FSpheroidalWaves.jl%2Fmaster%2FProject.toml&query=%24.compat.julia&suffix=%2B&label=Julia&color=purple)](https://julialang.org)
[![GitHub last commit](https://img.shields.io/github/last-commit/brandynlucca/SpheroidalWaves.jl?label=Last%20commit)](https://github.com/brandynlucca/SpheroidalWaves.jl/commits/main)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19728040.svg)](https://doi.org/10.5281/zenodo.19728040)


[![Documentation (stable)](https://img.shields.io/badge/docs-stable-blue?label=Package%20documentation%20(stable))](https://brandynlucca.github.io/SpheroidalWaves.jl/stable/)
[![Documentation (latest)](https://img.shields.io/badge/docs-latest-blue?label=Package%20documentation%20(latest))](https://brandynlucca.github.io/SpheroidalWaves.jl/dev/)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL3-green.svg?label=License)](LICENSE)

[![CI](https://img.shields.io/github/actions/workflow/status/brandynlucca/SpheroidalWaves.jl/CI.yml?branch=main&label=Build%20status)](https://github.com/brandynlucca/SpheroidalWaves.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/brandynlucca/SpheroidalWaves.jl/graph/badge.svg?token=ZH28KZ4DTQ)](https://codecov.io/gh/brandynlucca/SpheroidalWaves.jl)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FSpheroidalWaves&query=total_requests&label=SpheroidalWaves%20Downloads)](https://juliapkgstats.com/pkg/SpheroidalWaves)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FSpheroidalWaves_jll&query=total_requests&label=SpheroidalWaves_jll%20Downloads)](https://juliapkgstats.com/pkg/SpheroidalWaves_jll)

Fast, vectorized computation of prolate and oblate spheroidal wave functions using native Fortran kernels.

The package supports real and complex spheroidal parameters, double and quad precision, angular and radial functions, eigenvalues, accuracy estimates, derivatives, and inverse eigenvalue calculations.

## Installation

```julia
import Pkg
Pkg.add("SpheroidalWaves")
```

Prebuilt backends are downloaded automatically on supported platforms.

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

Pass arrays to `smn` and `rmn` for batch evaluation. Calls are [thread-safe](https://brandynlucca.github.io/SpheroidalWaves.jl/dev/math-and-usage/).

## Main API

- `smn`: angular function and derivative
- `rmn`: radial function and derivative
- `eigenvalue` and `eigenvalue_sweep`: characteristic values
- `accuracy`: estimated numerical accuracy
- `radial_wronskian`: radial Wronskian
- `jacobian_eigen`, `jacobian_smn`, and `jacobian_rmn`: derivatives with respect to parameters
- `find_c_for_eigenvalue`: inverse characteristic-value calculation

The main options are `spheroid=:prolate` or `:oblate`, `precision=:double` or `:quad`, and radial `kind=1:4`.

For basis conversion and connection formulas, `dmn`, `kmn`, and `amn` provide [expansion and normalization quantities](https://brandynlucca.github.io/SpheroidalWaves.jl/dev/mathematical-tools/).

See the [documentation](https://brandynlucca.github.io/SpheroidalWaves.jl/dev/) for definitions, argument restrictions, normalization, and complete examples.

## Attribution

The numerical kernels are based on the spheroidal-wave-function implementations by Arnie Lee Van Buren and Jeffrey Boisvert:

- [Prolate](https://github.com/MathieuandSpheroidalWaveFunctions/Prolate_swf)
- [Oblate](https://github.com/MathieuandSpheroidalWaveFunctions/Oblate_swf)
- [Complex prolate](https://github.com/MathieuandSpheroidalWaveFunctions/complex_prolate_swf)
- [Complex oblate](https://github.com/MathieuandSpheroidalWaveFunctions/complex_oblate_swf)

## Citation and license

Use the [Zenodo record](https://doi.org/10.5281/zenodo.19728040) to cite SpheroidalWaves.jl. The package is available under the [GPL-3.0 License](https://www.gnu.org/licenses/gpl-3.0.html). The logo incorporates Julia's dots and is licensed separately under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
