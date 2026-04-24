# SpheroidalWaveFunctions.jl

SpheroidalWaveFunctions.jl provides fast Julia bindings to native Fortran kernels for spheroidal wave computations.

## Public API

- `smn`
- `rmn`
- `eigenvalue`
- `accuracy`
- `radial_wronskian`
- `jacobian_eigen`
- `jacobian_smn`
- `jacobian_rmn`
- `find_c_for_eigenvalue`

## Quick Start

```julia
using SpheroidalWaveFunctions

eta = [-0.5, 0.0, 0.5]
x = [1.1, 1.2]

s = smn(0, 2, 20.0, eta; spheroid=:prolate, precision=:double)
r = rmn(0, 1, 20.0, x; spheroid=:prolate, precision=:double, kind=1)

lambda = eigenvalue(0, 2, 20.0; spheroid=:prolate, precision=:double)
acc = accuracy(0, 2, 20.0, eta; target=:angular)
W = radial_wronskian(0, 1, 20.0, x; spheroid=:prolate, precision=:double)

j_lambda = jacobian_eigen(0, 1, 20.0; with_metadata=true)
root = find_c_for_eigenvalue(0, 1, lambda; bracket=(1.0, 40.0))
```

See the `API` page for signatures and docstrings, and `Math and Usage` for equations and conventions.
