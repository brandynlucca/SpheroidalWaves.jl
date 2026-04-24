# Mathematical Basis and Usage

This page summarizes the mathematical model and practical usage for all public operations.

## Shared Separation Structure

Angular and radial spheroidal equations share the separation constant `lambda_mn(c)`.

- `m`: order with `m >= 0`
- `n`: degree with `n >= m`
- `c`: size parameter (real or complex)

## Angular Function `smn`

`smn` computes the angular spheroidal function and its derivative in `eta`.

```math
\frac{d}{d\eta}\left[(1-\eta^2)\frac{dS_{mn}}{d\eta}\right]
+\left(\lambda_{mn}(c)-c^2(1-\eta^2)-\frac{m^2}{1-\eta^2}\right)S_{mn}=0
```

## Radial Function `rmn`

`rmn` computes radial spheroidal functions and derivatives in `x` for selected `kind`.

```math
\frac{d}{dx}\left[(x^2-1)\frac{dR_{mn}}{dx}\right]
+\left(\lambda_{mn}(c)-c^2x^2-\frac{m^2}{x^2-1}\right)R_{mn}=0
```

## Wronskian Diagnostic `radial_wronskian`

The radial diagnostic uses:

```math
\mathscr{W}(x)=R_1(x)R_2'(x)-R_1'(x)R_2(x)
```

In a stable regime, the Wronskian should be approximately constant versus `x`.

## Eigenvalue `eigenvalue`

Returns `lambda_mn(c)` used in both separated equations.

At spherical limit:

```math
\lambda_{mn}(0)=n(n+1)
```

## Accuracy `accuracy`

Returns per-point backend-estimated decimal digits.

- `target=:angular`
- `target=:radial`

Use these values as solver-quality indicators for downstream acceptance criteria.

## Jacobians

The Jacobians are numerical derivatives with respect to `c`.

- `jacobian_eigen`
- `jacobian_smn`
- `jacobian_rmn`

For real `c`, centered difference uses:

```math
f'(c)\approx\frac{f(c+h)-f(c-h)}{2h}
```

For complex `c = a + ib`, partials are computed separately with respect to `a` and `b`.

### Reliability Metadata

With `with_metadata=true`, Jacobian APIs return derivative quality fields:

- `step_used`
- `relative_change_when_halving_step`
- `finite_flag`
- `conditioning_flag`
- `suggested_action`

These fields support robust root-finding and inverse workflows.

## Root Finding

`find_c_for_eigenvalue` solves the scalar inverse problem:

```math
\lambda_{mn}(c) - \lambda_{\text{target}} = 0
```

using a bracketed hybrid method with bisection fallback and optional Jacobian acceleration.

Example:

```julia
lambda_target = eigenvalue(0, 1, 30.0; spheroid=:prolate, precision=:double)
root = find_c_for_eigenvalue(0, 1, lambda_target;
							 bracket=(10.0, 50.0),
							 spheroid=:prolate,
							 precision=:double)
```

## Full Reference

For the complete long-form theory-and-usage narrative, see [docs/mathematical-basis-and-usage.md](https://github.com/brandynlucca/SpheroidalWaveFunctions.jl/blob/main/docs/mathematical-basis-and-usage.md) in the repository.
