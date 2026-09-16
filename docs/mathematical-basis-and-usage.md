# Mathematical Basis and Usage Guide

## Scope

This document describes the mathematical model behind every public API operation in `SpheroidalWaves.jl` and shows practical usage patterns.

See the [Mathematical Tools guide](src/mathematical-tools.md) for normalized
Wronskian errors, scaled outputs, and logarithmic and second derivatives.

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
\frac{d}{d\eta}\left[(1-\eta^2)\frac{dS_{mn}}{d\eta}\right] + \left(\lambda_{mn}(c) - \sigma c^2\eta^2 - \frac{m^2}{1-\eta^2}\right)S_{mn} = 0.
$$

Here `sigma=1` for prolate and `sigma=-1` for oblate functions. Angular values,
coordinate derivatives, and parameter Jacobians include the Condon–Shortley
phase `(-1)^m`. With `normalize=false`, the functions use Meixner–Schäfke
normalization. For complex parameters, the eigenmode and normalization continue
vertically from the real axis. Parameter Jacobians preserve the local branch.

At exactly `c=0`, both geometries use the associated Legendre limit for all
`0 <= m <= n`, with analytic endpoint derivatives. Unity normalization multiplies
the result by `sqrt((2n+1)(n-m)! / (2(n+m)!))`. See the
[angular convention notes](src/math-and-usage.md#angular-phase-and-spherical-limit)
for examples and endpoint behavior.

Returned values:

- `value`: `S_{mn}(eta, c)`
- `derivative`: `dS_{mn}/deta`

### Usage

```julia
using SpheroidalWaves

eta = [-0.5, 0.0, 0.5]
s = smn(0, 2, 20.0, eta; spheroid=:prolate, precision=:double)
```

## Radial Function: `rmn`

### Mathematical Basis

`rmn` computes radial spheroidal functions and first derivatives with respect to `x`.

For real prolate calls at `x=1`, kinds 2–4 return `NaN` for the undefined
value and derivative; `accuracy` reports `-1` there.

Radial `c=0` calls throw `DomainError` in both geometries and precisions, for
all four kinds. The same restriction applies to radial accuracy, Wronskians,
and parameter derivatives. The former real prolate `m=0` Legendre substitution
has been removed because it changed the radial normalization at exactly zero.
Angular functions and eigenvalues continue to support `c=0`.

A prolate radial equation form is:

$$
\frac{d}{dx}\left[(x^2-1)\frac{dR_{mn}}{dx}\right] + \left(c^2x^2 - \lambda_{mn}(c) - \frac{m^2}{x^2-1}\right)R_{mn} = 0.
$$

For oblate functions the radial equation is
`((x^2+1)*R′)′ + (c^2*x^2-lambda+m^2/(x^2+1))*R = 0`.

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

For real `c > 0`, standard radial normalization gives
`W(x) = 1 / (c * (x^2 - 1))` for prolate functions and
`W(x) = 1 / (c * (x^2 + 1))` for oblate functions.
The corresponding scaled quantity `c * (x^2 - 1) * W(x)` or
`c * (x^2 + 1) * W(x)` should be approximately one; the raw Wronskian is
not constant. Exactly zero is outside the supported radial domain.

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

`accuracy` returns integer diagnostics at each evaluation point: `-1` means no
estimate is available, `0` means no reliable decimal digits are reported, and
positive values are backend estimates. These are neither rigorous error bounds
nor statistical confidence intervals, and do not certify coordinate derivatives.

Angular `c=0` calls evaluate the recurrence and return `-1`. Radial `c=0` calls
throw `DomainError`. Input validation is shared with `smn` and `rmn`.
For nonzero radial calls, only `kind=2` reports the native second-kind estimate;
kinds `1`, `3`, and `4` are evaluated and return `-1`. Nonfinite values and invalid
native estimates also return `-1`. Quad diagnostics preserve the same input
precision as function evaluation. Check sensitive results against independent
references, identities, or precision comparisons.

### Usage

```julia
acc_s = accuracy(0, 2, 20.0, [-0.3, 0.0, 0.3]; target=:angular)
acc_r = accuracy(0, 2, 20.0, [1.1, 1.3]; target=:radial, kind=2)
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
using SpheroidalWaves

m, n = 0, 1
c = 60.0
eta = [-0.4, 0.0, 0.4]
x = [1.1, 1.3]

s = smn(m, n, c, eta; spheroid=:prolate, precision=:double)
r = rmn(m, n, c, x; spheroid=:prolate, precision=:double, kind=1)
lambda = eigenvalue(m, n, c; spheroid=:prolate, precision=:double)
acc_s = accuracy(m, n, c, eta; target=:angular)
acc_r = accuracy(m, n, c, x; target=:radial, kind=2)
W = radial_wronskian(m, n, c, x; spheroid=:prolate, precision=:double)

j_lambda = jacobian_eigen(m, n, c; with_metadata=true)
j_s = jacobian_smn(m, n, c, eta; with_metadata=true)
j_r = jacobian_rmn(m, n, c, x; kind=1, with_metadata=true)
```

## Related Documents

- `README.md` for package overview and quick start
- `docs/fortran-api-user-facing.md` for low-level Fortran argument mapping

