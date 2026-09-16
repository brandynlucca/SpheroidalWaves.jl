# Mathematical tools

All examples use ordinary numeric arguments. Add `precision=:quad` for quad precision.

## Eigenvalue operators

`eigenvalue(m, n, c; operator=:separation)` selects:

| `operator` | Result | `form` |
|:--|:--|:--|
| `:separation` (default) | Separation constant `λₘₙ(c)` | `:value` |
| `:concentration` | Concentration eigenvalue `Λₙ(c)` | `:value`, `:log`, `:complement` |
| `:fourier` | Finite Fourier eigenvalue `μₙ(c)` | `:value` |

### Separation constant

The angular and radial equations share `λₘₙ(c)`, with
`λₘₙ(0) = n(n+1)`. Both geometries and real or complex parameters are supported.

~~~julia
eigenvalue(1, 2, 1.25; precision=:quad)
eigenvalue(1, 2, 1.25; spheroid=:oblate, precision=:quad)
~~~

### Concentration and finite Fourier eigenvalues

For the first-kind angular mode `Sₙ = S₀ₙ` on `[-1,1]`:

~~~math
\int_{-1}^{1}\frac{\sin(c(x-t))}{\pi(x-t)}S_n(t)\,dt
=\Lambda_n(c)S_n(x),
\qquad
\int_{-1}^{1}e^{-icxt}S_n(t)\,dt=\mu_n(c)S_n(x).
~~~

The sinc kernel is `c/π` at `x=t`. For `c>0`,

~~~math
0<\Lambda_n<1,\qquad
\mu_n=(-i)^n\sqrt{\frac{2\pi\Lambda_n}{c}}.
~~~

These conventions follow [DLMF 30.15](https://dlmf.nist.gov/30.15)
with the interval scaled to `[-1,1]`.

~~~julia
eigenvalue(0, 2, 1.0; operator=:concentration, precision=:quad)
eigenvalue(0, 2, 1.0; operator=:fourier, precision=:quad)
eigenvalue(0, 0, 50.0; operator=:concentration,
    form=:complement, precision=:quad)
~~~

!!! note "Integral-operator domain"
    Requires `m=0`, `n≥0`, `spheroid=:prolate`, and finite real `c≥0`.
    At `c=0`, all `Λₙ=0`; `μ₀=2` and `μₙ=0` for `n>0`.

!!! tip "Small eigenvalues and complements"
    `form=:log` returns `log(Λₙ)`; `form=:complement` returns `1-Λₙ`.
    Use these directly when `Λₙ` rounds to zero or one.
    At `c=0`, they return `-Inf` and `1`, respectively.

## Parameter derivatives

The Jacobians differentiate with respect to `c`; coordinate derivatives are
returned by `smn` and `rmn`.

| Call | Result for real `c` |
|:--|:--|
| `jacobian_eigen(m,n,c)` | `dλ/dc` |
| `jacobian_smn(m,n,c,eta)` | `dvalue_dc = ∂S/∂c`, `dderivative_dc = ∂²S/(∂c∂η)` |
| `jacobian_rmn(m,n,c,x)` | `dvalue_dc = ∂R/∂c`, `dderivative_dc = ∂²R/(∂c∂x)` |

~~~julia
jacobian_smn(1, 2, 1.25, [0.3]; kind=2, precision=:quad)
jacobian_rmn(1, 2, 1.25, [2.0]; kind=2, precision=:quad)
jacobian_eigen(0, 2, 1.0; operator=:concentration, precision=:quad)
jacobian_eigen(0, 2, 1.0; operator=:fourier, precision=:quad)
~~~

`jacobian_eigen` accepts the same operators and forms as `eigenvalue`.
For a real mode `ψₙ` with unit norm on `[-1,1]` and `c>0`:

~~~math
\frac{d\Lambda_n}{dc}=\frac{|\mu_n|^2}{\pi}\psi_n(1)^2,
\qquad
\frac{d\log\Lambda_n}{dc}=\frac{2\psi_n(1)^2}{c},
\qquad
\frac{d(1-\Lambda_n)}{dc}=-\frac{d\Lambda_n}{dc}.
~~~

For complex `c=a+ib`, results contain separate real- and imaginary-direction
partials (`d_dcreal`/`d_dcimag` for eigenvalues;
`dvalue_dcreal`/`dvalue_dcimag` and
`dderivative_dcreal`/`dderivative_dcimag` for functions).
On a local analytic branch, `∂f/∂b = i ∂f/∂a`.

!!! note "Derivative options and limits"
    Defaults use analytic parameter differentiation; supplying `h` selects finite differences.
    `with_metadata=true` adds convergence diagnostics, not error bounds.
    Angular `kind=2` parameter derivatives return `NaN` at `η=±1`.
    Radial calls require `c≠0`.

!!! note "Integral derivatives at zero"
    Right-hand limits are `Λ₀′=2/π`, `Λₙ′=0` for `n>0`,
    and `μ₁′=-2im/3` (all other Fourier derivatives vanish).
    `form=:log` returns `Inf`. An explicit finite-difference step
    requires `h<c` for positive bandwidth; at zero, logarithmic finite differences are unsupported.

## Inverse bandwidth

`find_c_for_eigenvalue` solves for real `c` within a supplied bracket:

~~~math
\lambda_{mn}(c)=\lambda_{\mathrm{target}}
\quad\text{or}\quad
\Lambda_n(c)=\Lambda_{\mathrm{target}}.
~~~

~~~julia
target = eigenvalue(1, 2, 1.25; precision=:quad)
root = find_c_for_eigenvalue(1, 2, target;
    bracket=(0.5, 2.0), precision=:quad)

root = find_c_for_eigenvalue(0, 0, 0.999;
    operator=:concentration, bracket=(0.0, 10.0), precision=:quad)
root.c
root.converged
root.residual
~~~

For concentration, `form` defines the target:

| Form | Target | Allowed range |
|:--|:--|:--|
| `:value` | `Λₙ` | `[0,1)` |
| `:complement` | `1-Λₙ` | `(0,1]` |
| `:log` | `log(Λₙ)` | `[-Inf,0)` |

!!! warning "Check convergence"
    Check `root.converged` before using `root.c`; `root.residual` uses the requested form.
    Unit concentration has no finite bandwidth. Zero concentration requires a bracket
    containing zero. Fourier inversion is unsupported.
    Set `atol` and `rtol` explicitly for tighter separation-constant inversion.

## Eigenvalue sweeps

~~~julia
sweep = eigenvalue_sweep(1, 2, [0.5, 1.0, 1.5]; precision=:quad)
sweep.lambda

path = ComplexF64[2, 2+3im, 3im, 0, 2]
continued = eigenvalue_sweep(0, 0, path; precision=:quad)
continued.selected_n
~~~

Sweeps track separation constants only. Real grids must be monotone; complex
grids specify an ordered piecewise-linear path.
`selected_n` records degree labels and `switched_branch` flags label changes.
`branch_lock=false` evaluates the native label independently at each point.

!!! warning "Complex branches"
    A closed path can exchange eigenmodes. Its final value need not equal its initial value.
    Nearly coalescing eigenvalues can make continuation and parameter derivatives ill-conditioned.

## Multiple degrees

~~~julia
smn(1, 1:24, 1.25, [-0.7, 0.2, 0.8];
    spheroid=:oblate, precision=:quad, normalize=true)
rmn(1, 1:24, 1.25, [0.4, 1.0, 3.0];
    spheroid=:oblate, precision=:quad, kind=3, scaled=true)
~~~

Rows correspond to coordinates; columns correspond to degrees.
A scalar coordinate gives one row. Single-degree calls return vectors,
including length-one vectors for scalar coordinates.

## Scaled values and logarithmic derivatives

~~~math
f=M\,10^E,\qquad \ell=\frac{f'}{f}.
~~~

~~~julia
s = smn(200, 200, 0.0, 0.3; scaled=true, logderivative=true)
s.value.mantissa
s.value.exponent
s.logderivative

r = rmn(1, 2, 1.25, 2.0;
    kind=2, precision=:quad, scaled=true, logderivative=true)
~~~

With `scaled=true`, `value` and `derivative` each contain
`mantissa` and integer `exponent` arrays.
For finite nonzero values, `1≤abs(mantissa)<10`.
Complex mantissas use a common exponent for both components;
zeros and nonfinite values use exponent zero.

!!! warning "Range and zeros"
    Scaling extends the representable magnitude, not the number of accurate digits.
    `logderivative=true` returns `NaN` at a computed zero; nearby ratios are ill-conditioned.
    Unscaled values may overflow even when the logarithmic derivative remains finite.

## Second coordinate derivatives

For `λ=λₘₙ(c)` and `σ=1` (prolate) or `σ=-1` (oblate):

~~~math
S''=\frac{2\eta S'-
\left(\lambda-\sigma c^2\eta^2-\frac{m^2}{1-\eta^2}\right)S}{1-\eta^2},
~~~

~~~math
R''=\frac{-2xR'-
\left(c^2x^2-\lambda-\frac{\sigma m^2}{x^2-\sigma}\right)R}{x^2-\sigma}.
~~~

~~~julia
s = smn(1, 2, 1.25, 0.3; second_derivative=true, precision=:quad)
r = rmn(1, 2, 1.25, 2.0; second_derivative=true, precision=:quad)
s.second_derivative
r.second_derivative
~~~

Both options also work with scaled outputs and degree ranges.

!!! note "Singular endpoints"
    The formulas above apply away from singular endpoints.
    First-kind angular and regular first-kind prolate radial endpoints use one-sided limits:
    finite for `m=0,2,4`, signed infinities for `m=1,3`, and zero for `m>4`.
    Angular second-kind derivatives diverge at `η=±1`.
    Undefined prolate radial kinds 2–4 at `x=1` return `NaN`;
    complex prolate calls require `x>1`.

## Wronskians and accuracy

For real `c>0` and nonsingular radial coordinates:

~~~math
W=R_1R_2'-R_1'R_2,\qquad
\widehat W=c(x^2-\sigma)W=1,\qquad
\varepsilon_W=|\widehat W-1|.
~~~

| `radial_wronskian` form | Output |
|:--|:--|
| `:raw` (default) | `W` |
| `:normalized` | `Ŵ` |
| `:error` | `ε_W` |

~~~julia
radial_wronskian(1, 2, 1.25, [1.5, 2.0, 4.0];
    precision=:quad, form=:error)
accuracy(0, 1, 2.0, [1.5]; target=:radial, kind=2)
accuracy(0, 2, 0.0, [0.3]; target=:angular)
~~~

The prolate identity is [DLMF 30.11.7](https://dlmf.nist.gov/30.11#E7);
the oblate factor follows its radial equation and normalization.

`accuracy` returns estimated decimal digits for function values:
`-1` means unavailable, `0` means no reliable digits reported.
Only radial `kind=2` has a digit estimate where available.
Angular `c=0` and `kind=2` return `-1`.

!!! warning "Diagnostics are not error bounds"
    A small Wronskian error checks consistency of the solution pair.
    Correlated errors can preserve that identity.
    Digit estimates do not certify coordinate derivatives.


