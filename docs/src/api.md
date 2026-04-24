# API Reference

This page provides a concise signature-level reference for the public API.

## Core Evaluation

- `smn(m, n, c, eta; spheroid=:prolate, precision=:double, normalize=false)`
- `rmn(m, n, c, x; spheroid=:prolate, precision=:double, kind=1)`

## Analysis and Diagnostics

- `eigenvalue(m, n, c; spheroid=:prolate, precision=:double)`
- `accuracy(m, n, c, arg; spheroid=:prolate, precision=:double, kind=1, target=:radial, normalize=false)`
- `radial_wronskian(m, n, c, x; spheroid=:prolate, precision=:double)`

## Jacobians

- `jacobian_eigen(m, n, c; spheroid=:prolate, precision=:double, h=nothing, with_metadata=false, adaptive=true, rtol=1e-6, atol=1e-10)`
- `jacobian_smn(m, n, c, eta; spheroid=:prolate, precision=:double, normalize=false, h=nothing, with_metadata=false, adaptive=true, rtol=1e-6, atol=1e-10)`
- `jacobian_rmn(m, n, c, x; spheroid=:prolate, precision=:double, kind=1, h=nothing, with_metadata=false, adaptive=true, rtol=1e-6, atol=1e-10)`

## Root Finding

- `find_c_for_eigenvalue(m, n, lambda_target; bracket=(c_lo, c_hi), spheroid=:prolate, precision=:double, atol=1e-10, rtol=1e-8, maxiter=80, use_jacobian=true)`

For mathematical definitions and worked usage examples, see the `Math and Usage` page.
