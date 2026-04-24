# Mathematical Basis and Usage Guide

## Scope

This document describes the mathematical model behind every public API operation in `SpheroidalWaveFunctions.jl` and shows practical usage patterns.

Public operations covered:

- `smn`
- `rmn`
- `radial_wronskian`
- `eigenvalue`
- `accuracy`
- `jacobian_eigen`
- `jacobian_smn`
- `jacobian_rmn`
- `find_c_for_eigenvalue`

## Notation and Problem Setup

The spheroidal wave equation separates in spheroidal coordinates into angular and radial ordinary differential equations that share the same separation constant (eigenvalue) `lambda`.

- `m`: order, integer mode index
- `n`: degree, integer with `n >= m >= 0`
- `c`: spheroidal size parameter (real or complex)
- `eta`: angular coordinate in `[-1, 1]`
- `x`: radial coordinate (`x > 1` for the real prolate radial regime)
- `lambda_mn(c)`: separation constant for `(m, n, c)`

For real `c`, the library computes real-family spheroidal functions. For complex `c`, it computes complex-family spheroidal functions.

## Angular Function: `smn`

### Mathematical Basis

`smn` computes angular spheroidal functions of the first kind and their first derivatives with respect to `eta`:

$$
\frac{d}{d\eta}\left[(1-\eta^2)\frac{dS_{mn}}{d\eta}\right] + \left(\lambda_{mn}(c) - c^2(1-\eta^2) - \frac{m^2}{1-\eta^2}\right)S_{mn} = 0.
$$

Returned values:

- `value`: `S_{mn}(eta, c)`
- `derivative`: `dS_{mn}/deta`

### Usage

```julia
using SpheroidalWaveFunctions

eta = [-0.5, 0.0, 0.5]
s = smn(0, 2, 20.0, eta; spheroid=:prolate, precision=:double)
```

## Radial Function: `rmn`

### Mathematical Basis

`rmn` computes radial spheroidal functions and first derivatives with respect to `x`.

A radial equation form is:

$$
\frac{d}{dx}\left[(x^2-1)\frac{dR_{mn}}{dx}\right] + \left(\lambda_{mn}(c) - c^2x^2 - \frac{m^2}{x^2-1}\right)R_{mn} = 0.
$$

The `kind` selector chooses the radial family:

- `kind=1`: first kind
- `kind=2`: second kind
- `kind=3`, `kind=4`: linear combinations analogous to outgoing/incoming constructions

Returned values are complex-valued vectors for consistent representation across kinds.

### Usage

```julia
x = [1.1, 1.3, 1.6]
r = rmn(0, 1, 40.0, x; spheroid=:prolate, precision=:double, kind=1)
```

## Wronskian Diagnostic: `radial_wronskian`

### Mathematical Basis

For two linearly independent radial solutions `R1` and `R2`, the Wronskian is

$$
\mathscr{W}(x) = R_1(x)R_2'(x) - R_1'(x)R_2(x).
$$

For a correctly computed solution pair under fixed normalization, $\mathscr{W}(x)$ should be approximately constant in $x$. This is a consistency diagnostic for numerical stability.

### Usage

```julia
x = [1.1, 1.2, 1.4]
W = radial_wronskian(0, 1, 40.0, x; spheroid=:prolate, precision=:double)
```

## Separation Constant: `eigenvalue`

### Mathematical Basis

`eigenvalue` returns the separation constant `lambda_mn(c)` used by both angular and radial equations for the selected family.

Spherical limit (`c = 0`) is

$$
\lambda_{mn}(0) = n(n+1)
$$

for the conventions used by this package.

### Usage

```julia
lambda_real = eigenvalue(0, 2, 20.0; spheroid=:prolate, precision=:double)
lambda_complex = eigenvalue(0, 2, 20.0 + 0.2im; spheroid=:oblate, precision=:double)
```

## Accuracy Estimate: `accuracy`

### Mathematical Basis

`accuracy` returns backend-estimated decimal digits of reliability at each evaluation point.

- `target=:angular` reports angular accuracy estimates
- `target=:radial` reports radial accuracy estimates

This is not a strict interval bound; it is a solver quality estimate useful for gating downstream workflows.

### Usage

```julia
acc_s = accuracy(0, 2, 20.0, [-0.3, 0.0, 0.3]; target=:angular)
acc_r = accuracy(0, 2, 20.0, [1.1, 1.3]; target=:radial, kind=1)
```

## Jacobian APIs

All Jacobians are numerical derivatives with respect to parameter `c`.

- Real `c`: derivative with respect to scalar `c`
- Complex `c = a + ib`: partials with respect to `a` and `b`

Centered finite differences are used. With reliability mode enabled, derivatives are checked against a halved step and optionally refined.

### `jacobian_eigen`

#### Mathematical Basis

For real `c`:

$$
\frac{d\lambda}{dc} \approx \frac{\lambda(c+h)-\lambda(c-h)}{2h}
$$

For complex `c = a + ib`:

$$
\frac{\partial\lambda}{\partial a} \approx \frac{\lambda(a+h+ib)-\lambda(a-h+ib)}{2h},
\quad
\frac{\partial\lambda}{\partial b} \approx \frac{\lambda(a+i(b+h))-\lambda(a+i(b-h))}{2h}.
$$

#### Usage

```julia
j = jacobian_eigen(0, 1, 40.0; h=1e-6)
jm = jacobian_eigen(0, 1, 40.0 + 0.1im; h=1e-6, with_metadata=true)
```

### `jacobian_smn`

#### Mathematical Basis

Computes derivatives of both outputs from `smn` with respect to `c`.

For real `c`:

$$
\frac{\partial S}{\partial c}, \quad \frac{\partial}{\partial c}\left(\frac{dS}{d\eta}\right)
$$

approximated with centered finite differences.

#### Usage

```julia
eta = [0.0, 0.2]
js = jacobian_smn(0, 1, 40.0, eta; h=1e-6)
```

### `jacobian_rmn`

#### Mathematical Basis

Computes derivatives of both outputs from `rmn` with respect to `c`:

$$
\frac{\partial R}{\partial c}, \quad \frac{\partial}{\partial c}\left(\frac{dR}{dx}\right)
$$

using centered finite differences.

#### Usage

```julia
x = [1.1, 1.2]
jr = jacobian_rmn(0, 1, 40.0, x; kind=1, h=1e-6)
```

## Root Finding: `find_c_for_eigenvalue`

### Mathematical Basis

Solves for real `c` in

```math
\lambda_{mn}(c) - \lambda_{\text{target}} = 0
```

using a bracketed hybrid method with guaranteed bisection fallback and optional
Jacobian-guided acceleration.

### Usage

```julia
lambda_target = eigenvalue(0, 1, 30.0; spheroid=:prolate, precision=:double)
root = find_c_for_eigenvalue(0, 1, lambda_target;
                             bracket=(10.0, 50.0),
                             spheroid=:prolate,
                             precision=:double)
```

## Jacobian Reliability Metadata

Set `with_metadata=true` on Jacobian APIs to receive derivative quality diagnostics.

Metadata fields:

- `step_used`: step size used for reported derivative
- `relative_change_when_halving_step`: consistency metric between step levels
- `finite_flag`: `true` when derivative estimate is finite
- `conditioning_flag`: one of `:good`, `:warning`, `:poor`
- `suggested_action`: one of `:accept`, `:retry_smaller_h`, `:use_quad`

Example:

```julia
jm = jacobian_eigen(0, 1, 40.0; with_metadata=true)
if jm.metadata.suggested_action != :accept
    # tighten step, or switch precision=:quad
end
```

## Parameter Regimes and Practical Guidance

- Use `precision=:double` for standard workloads.
- Use `precision=:quad` for sensitive regimes, high `|c|`, or when Jacobian metadata indicates poor conditioning.
- Validate radial domains physically before solve calls.
- For inverse or root workflows, use Jacobian metadata and `accuracy` together to decide acceptance.

## Spherical Limit

When `c = 0`, several quantities reduce to classical spherical values:

- `lambda_mn(0) = n(n+1)`
- angular functions reduce to associated Legendre-family behavior under package conventions
- radial behavior follows corresponding spherical-function limits

These identities are used in package tests for baseline correctness checks.

## End-to-End Example

```julia
using SpheroidalWaveFunctions

m, n = 0, 1
c = 60.0
eta = [-0.4, 0.0, 0.4]
x = [1.1, 1.3]

s = smn(m, n, c, eta; spheroid=:prolate, precision=:double)
r = rmn(m, n, c, x; spheroid=:prolate, precision=:double, kind=1)
lambda = eigenvalue(m, n, c; spheroid=:prolate, precision=:double)
acc_s = accuracy(m, n, c, eta; target=:angular)
acc_r = accuracy(m, n, c, x; target=:radial, kind=1)
W = radial_wronskian(m, n, c, x; spheroid=:prolate, precision=:double)

j_lambda = jacobian_eigen(m, n, c; with_metadata=true)
j_s = jacobian_smn(m, n, c, eta; with_metadata=true)
j_r = jacobian_rmn(m, n, c, x; kind=1, with_metadata=true)
```

## Related Documents

- `README.md` for package overview and quick start
- `docs/fortran-api-user-facing.md` for low-level Fortran argument mapping
