# Mathematical Basis and Usage

This page summarizes the mathematical model and practical usage for all public operations.

## Shared Separation Structure

Angular and radial spheroidal equations share the separation constant `lambda_mn(c)`.

- `m`: order with `m >= 0`
- `n`: degree with `n >= m`
- `c`: size parameter (real or complex)

### Precision

With `precision=:quad`, real angular values and eigenvalues are `BigFloat`;
complex angular/radial values, their coordinate derivatives, and complex
eigenvalues are `Complex{BigFloat}`. Real-parameter radial results also use
`Complex{BigFloat}`. Ordinary numeric arguments and arrays are accepted with
either precision setting; manual conversion is not required.

Quad precision does not recover digits already lost in the supplied inputs,
and it does not guarantee 30 accurate digits in every parameter regime.

### Numerical limits

Large degree and bandwidth can produce values much smaller or larger than
individual expansion terms. Use `scaled=true` to retain decimal exponents and
`logderivative=true` to form derivative/value ratios before output rounding.
Scaling extends the output range; it does not certify accuracy. At a computed
zero, the logarithmic derivative is `NaN`; near a zero it is ill-conditioned.

`accuracy` returns `-1` when no estimate is available. Nearly coalescing
complex eigenvalues and singular endpoints can make evaluation ill-conditioned.

## Angular Function `smn`

`smn` computes the angular spheroidal function and its derivative in `eta`.

```math
\frac{d}{d\eta}\left[(1-\eta^2)\frac{dS_{mn}}{d\eta}\right]
+\left(\lambda_{mn}(c)-\sigma c^2\eta^2-\frac{m^2}{1-\eta^2}\right)S_{mn}=0
```

Here `sigma=1` for prolate and `sigma=-1` for oblate functions, using the
package's separation constant convention.

### Angular phase and spherical limit

This subsection describes the default `kind=1`. For the second kind, see below.

Angular values and their derivatives include the **Condon–Shortley phase**
`(-1)^m`, for both choices of `normalize`. With `normalize=false`, the functions
use Meixner–Schäfke normalization.

At exactly `c=0`, both geometries return associated Legendre (Ferrers) functions
`P_n^m(eta)` and their derivatives for every `0 <= m <= n`.
With `normalize=true`, the squared integral over `[-1,1]` is one; the scale
relative to `P_n^m` is `sqrt((2n+1)(n-m)! / (2(n+m)!))`.

The precise references are [DLMF equation 30.4.2](https://dlmf.nist.gov/30.4#E2)
for the exact-zero identity and the text immediately after
[equation 30.4.1](https://dlmf.nist.gov/30.4#E1) for the phase. DLMF uses
`gamma^2` as its parameter: `gamma^2 = c^2` for prolate functions and
`gamma^2 = -c^2` for oblate functions. When `n-m` is even, the sign at `eta=0`
is `(-1)^((n+m)/2)`; when `n-m` is odd, the derivative there has sign
`(-1)^((n+m-1)/2)` for real `gamma^2`.

For complex `c`, `eigenvalue`, `smn`, and `rmn` follow degree `n` from `real(c)`
along the vertical segment to `c`, preserving the angular normalization and
phase continuously. The branch selection is independent of the requested
coordinates and applies to values, derivatives, and accuracy estimates.

Use `eigenvalue_sweep` with a complex vector to specify an ordered
piecewise-linear path. It initializes its first point using the scalar rule,
then continues from each preceding point. Complex paths need not be monotone;
repeated points and closed loops are allowed. `selected_n` reports native labels,
and `switched_branch` marks changes of those labels. With `branch_lock=false`,
the sweep instead returns independently evaluated native labels.

```julia
path = [2, 2+3im, 3im, 0, 2]
s = eigenvalue_sweep(0, 0, path; precision=:quad)
s.selected_n  # [0, 2, 2, 2, 2]
# The endpoint follows the n=2 branch, despite returning to c=2.
```

Complex eigenvalues have branch points; see
[Skorokhodov and Khristoforov (2006)](https://www.mathnet.ru/eng/zvmmf438).
A closed path can exchange branches, so its endpoint need not equal a fresh
scalar call using the canonical vertical path. This is numerical continuation
along a specified path, not a globally single-valued definition. Paths that
cannot be resolved within the subdivision and candidate limits raise an error.
Passing near a coalescing eigenvalue can be ill-conditioned; successful
continuation does not certify accuracy there or agreement with another system's
branch cuts.

```julia
smn(1, 1, 0.3, 0.0).value  # approximately [-1.001795214382128]
smn(1, 1, 0.0, 0.3).value  # approximately [-0.953939201416946]
smn(0, 2, 0.0, [-1.0, 1.0]).derivative  # [-3.0, 3.0]
```

Endpoint derivatives use their one-sided limits. For `m=1` at `c=0`, these
derivatives are infinite and are returned as signed infinities.

### Angular functions of the second kind

Pass `kind=2` to `smn`:

```julia
q = smn(1, 2, 1.25, [-0.3, 0.3]; kind=2, precision=:quad)
q.value
q.derivative
```

This is the Ferrers-branch function `Qs` on `-1 < eta < 1`, as defined in
[DLMF 30.5](https://dlmf.nist.gov/30.5) and
[30.8(ii)](https://dlmf.nist.gov/30.8#ii). The parameter is `gamma^2 = sigma*c^2`.
At `c=0` it is the associated Ferrers function `Q_n^m`, including the
Condon–Shortley phase. Complex parameters follow the same eigenvalue and
normalization branch as the first kind; complex coordinates are not accepted.

The normalization is fixed by opposite parity and the Wronskian:

```math
Q(-\eta)=(-1)^{n-m+1}Q(\eta),\qquad
(1-\eta^2)(P Q'-P'Q)=\frac{(n+m)!}{(n-m)!}A_n^m A_n^{-m},
```

where the coefficient sums `A` are defined in
[DLMF 30.11.4](https://dlmf.nist.gov/30.11#E4), and `P` is the default,
unnormalized first-kind result. `normalize=true` is rejected for `kind=2`;
the first kind's unit-integral normalization does not define this singular family.
For `m>=1`, the squared magnitude is not integrable at the endpoints. For `m=0`
it is integrable, but its integral differs from the first-kind normalization.
Multiplying `Qs` by a chosen finite scale is valid; it is not the unit-integral
operation represented by this keyword.

Both values and derivatives diverge at `eta=±1`: logarithmically for the value
when `m=0`, and algebraically when `m>0`. Exact endpoints return signed
one-sided infinities (componentwise for complex results). Nearby interior
coordinates remain evaluable in both precisions. `second_derivative=true`,
`scaled=true`, `logderivative=true`, and degree ranges are supported.
The logarithmic derivative returns `NaN` at a computed zero and has limits
`-Inf` and `Inf` at the left and right endpoints, respectively.

Singular or numerically unresolved second-kind normalization raises an error.
`accuracy(...; target=:angular, kind=2)` evaluates the requested function and
returns `-1`, because no calibrated decimal-digit estimate is available.
Parameter Jacobians support both kinds through `jacobian_smn(...; kind=2)`.
Angular zero finding currently concerns the first kind.

## Radial Function `rmn`

`rmn` computes radial spheroidal functions and derivatives in `x` for selected `kind`.

For real prolate calls at `x=1`, kinds 2–4 return `NaN` for the undefined
value and derivative; `accuracy` reports `-1` there.

Radial calls require `c != 0`. Exactly zero throws `DomainError` for all kinds,
both geometries, and both precisions, including complex-valued zero.
`radial_wronskian`, `jacobian_rmn`, and radial `accuracy` enforce the same domain.
Earlier versions substituted Legendre functions at zero for real prolate `m=0`;
that substitution used a different normalization and has been removed.
The angular `smn` and eigenvalue zero-parameter cases remain supported.

For prolate functions:

```math
\frac{d}{dx}\left[(x^2-1)\frac{dR_{mn}}{dx}\right]
+\left(c^2x^2-\lambda_{mn}(c)-\frac{m^2}{x^2-1}\right)R_{mn}=0
```

For oblate functions, with the package's oblate eigenvalue:

```math
\frac{d}{dx}\left[(x^2+1)\frac{dR_{mn}}{dx}\right]
+\left(c^2x^2-\lambda_{mn}(c)+\frac{m^2}{x^2+1}\right)R_{mn}=0.
```

See [Mathematical Tools](mathematical-tools.md) for scaled outputs, logarithmic
derivatives, second derivatives, and coordinate zeros.

## Wronskian Diagnostic `radial_wronskian`

The radial diagnostic uses:

```math
\mathscr{W}(x)=R_1(x)R_2'(x)-R_1'(x)R_2(x)
```

For real `c > 0`, the standard radial normalization gives
`W(x) = 1 / (c * (x^2 - 1))` for prolate functions and
`W(x) = 1 / (c * (x^2 + 1))` for oblate functions.
Multiply by `c * (x^2 - 1)` or `c * (x^2 + 1)`, respectively, to obtain
a diagnostic that should be approximately one. The raw Wronskian varies with `x`.
Exactly zero is outside the supported radial domain.

Use `form=:normalized` for the quantity with target one, or `form=:error` for
its absolute difference from one. The default `form=:raw` preserves the original
return value. These are consistency diagnostics, not accuracy bounds.

## Eigenvalue `eigenvalue`

The default `operator=:separation` returns `lambda_mn(c)` used in both separated equations.

At spherical limit:

```math
\lambda_{mn}(0)=n(n+1)
```

### Concentration and finite Fourier operators

For `m=0`, `spheroid=:prolate`, and real `c>=0`, `eigenvalue` also provides the
integral-operator eigenvalues associated with the same angular modes. With
`S_n(x)=smn(0,n,c,x).value[1]`, their definitions on `[-1,1]` are

```math
\int_{-1}^{1}\frac{\sin(c(x-t))}{\pi(x-t)}S_n(t)\,dt
=\Lambda_n(c)S_n(x),
```

```math
\int_{-1}^{1}e^{-icxt}S_n(t)\,dt=\mu_n(c)S_n(x).
```

The sinc kernel has diagonal value `c/pi`. For `c>0`, the concentration
eigenvalues satisfy `0<Λ_n<1`. The Fourier convention uses the negative
exponential: `μ_n=(-im)^n*sqrt(2pi*Λ_n/c)`. These definitions follow
[DLMF 30.15.3](https://dlmf.nist.gov/30.15#E3) and
[30.15.5](https://dlmf.nist.gov/30.15#E5), with the interval scaled to `[-1,1]`.
Multiplying an angular mode by a nonzero normalization factor does not change
either integral eigenvalue.

```julia
eigenvalue(0, 2, 1.0; operator=:concentration, precision=:quad)
eigenvalue(0, 2, 1.0; operator=:fourier, precision=:quad)

# Compute the small quantity before rounding a near-unit eigenvalue.
eigenvalue(0, 0, 50.0; operator=:concentration,
           form=:complement, precision=:quad) # approximately 1.84851342e-42

# Retain the logarithm when ordinary double output would underflow.
eigenvalue(0, 10, 1e-50; operator=:concentration, form=:log)
```

For concentration, `form=:value` returns `Λ_n`, `form=:log` returns its natural
logarithm, and `form=:complement` returns `1-Λ_n`. The latter two are calculated
at extra working precision before conversion, so they can retain information
lost by calling `log` or subtracting from one on an already rounded result.
Concentration outputs are `Float64` or `BigFloat`; Fourier outputs are
`ComplexF64` or `Complex{BigFloat}`. The Fourier and separation operators accept
only `form=:value`. All outputs remain subject to their type's exponent range;
the logarithmic form is useful for especially small concentration eigenvalues.

At `c=0`, all concentration eigenvalues are exactly zero (`log(Λ_n)=-Inf`,
`1-Λ_n=1`). The Fourier eigenvalue is 2 for `n=0` and zero for higher modes.
These integral limits are supported independently of the radial `c=0` domain.
Oblate geometry, nonzero order, negative bandwidth, and complex bandwidth are
outside the supported integral-operator domain.

The calculation refines the Legendre coefficients and working precision and
checks convergence of both `Λ_n` and `1-Λ_n`. It uses the Fourier identity at
the origin for even modes and its coordinate derivative there for odd modes;
the required integrals reduce to the lowest Legendre coefficient. Failure to
resolve a positive concentration and complement raises an error.

`jacobian_eigen` supports the same operators and forms. `find_c_for_eigenvalue`
supports separation and concentration; `eigenvalue_sweep` uses separation constants.

### Bandwidth derivatives and target concentration

```julia
jacobian_eigen(0, 2, 1.0; operator=:concentration, precision=:quad)
jacobian_eigen(0, 2, 1.0; operator=:fourier, precision=:quad)
jacobian_eigen(0, 0, 50.0; operator=:concentration,
               form=:complement, precision=:quad)

# Bandwidth needed for 99.9% concentration in the lowest mode.
root = find_c_for_eigenvalue(0, 0, 0.999;
    operator=:concentration, bracket=(0.0, 10.0), precision=:quad)
root.c
root.converged

# Specify the tail directly when subtraction from one would lose it.
find_c_for_eigenvalue(0, 0, 1e-40;
    operator=:concentration, form=:complement,
    bracket=(1.0, 60.0), precision=:quad)
```

For a unit-norm real angular mode `ψ_n`, differentiating the sinc kernel and
using the finite Fourier identity gives

```math
\frac{d\Lambda_n}{dc}=\frac{|\mu_n|^2}{\pi}\psi_n(1)^2,
\qquad
\frac{d\log\Lambda_n}{dc}=\frac{2\psi_n(1)^2}{c}.
```

These are consequences of the integral definitions above. The calculation
uses refined coefficients and extra working precision to retain small endpoint
values. Fourier derivatives use the same identity, with coefficient
sensitivities near zero to avoid cancellation. `form=:complement` returns
`-dΛ_n/dc`; `form=:log` returns `d(log Λ_n)/dc`.

At `c=0`, derivatives are right-hand limits: `Λ_0′=2/pi`, `Λ_n′=0` for `n>0`,
`μ_1′=-2im/3`, and all other Fourier derivatives are zero. The logarithmic
derivative tends to positive infinity for every degree and returns `Inf`.
`with_metadata=true` distinguishes analytic identities and endpoint limits;
it flags the infinite limit as singular. Supplying `h` requests finite
differences: a centered stencil with `h<c`, or a forward stencil at zero.
Finite differences of the logarithm at zero are rejected.

Concentration inversion accepts a value in `[0,1)`, a complement in `(0,1]`,
or a logarithm in `[-Inf,0)`. Zero concentration requires a bracket containing
zero; unit concentration has no finite bandwidth. The solver brackets the root
and uses analytic Newton steps with bisection fallback. It checks
`log(Λ/(1-Λ))` rather than an absolute eigenvalue residual alone, so tiny
targets and near-unit targets remain distinguishable. The returned `residual`
uses the requested form. Inspect `converged` before using the returned bandwidth.

Concentration defaults are `atol=1e-12, rtol=1e-10` in double and
`atol=1e-30, rtol=1e-28` in quad. The coordinate tolerance is
`atol+rtol*abs(c)`; the residual in `log(Λ/(1-Λ))` must also be at most `rtol`.
Both inverse operators preserve `BigFloat` coordinates, residuals, and brackets
with `precision=:quad`. Separation keeps its existing, looser default
tolerances; specify tighter tolerances when additional digits are required.

## Accuracy `accuracy`

Returns a vector of integer diagnostics for the requested function values:

- `-1`: no estimate is available.
- `0`: the backend reports no reliable decimal digits.
- Positive values: estimated decimal digits, not guaranteed accuracy.

These diagnostics are neither rigorous error bounds nor statistical confidence
intervals. They do not certify coordinate derivatives. Use independent references,
identities, or precision comparisons when accuracy matters.

Angular `c=0` calls evaluate the associated Legendre recurrence and return `-1`,
because that path has no error estimator. Radial `c=0` calls throw `DomainError`,
just as `rmn` does. Invalid orders, degrees, and coordinates are rejected by both
the evaluation and diagnostic APIs.

For radial functions, only `kind=2` exposes the native second-kind estimate.
Kinds `1`, `3`, and `4` are evaluated but return `-1`: the second-kind diagnostic
does not establish the accuracy of the first kind or of a Hankel combination.
Nonfinite values and unavailable or out-of-range estimates also return `-1`.
Quad diagnostics use the same precision-preserving inputs as function evaluation.

```julia
accuracy(0, 2, 0.0, [0.3]; target=:angular) # [-1]
accuracy(0, 1, 2.0, [1.5]; target=:radial, kind=2)
```

## Jacobians

The Jacobians differentiate with respect to the parameter `c`:

- `jacobian_eigen` differentiates the Legendre coefficient eigenproblem.
- `jacobian_smn` evaluates the differentiated coefficient expansion, preserving
  the requested angular normalization. It reuses the same coefficient calculation.
- With `kind=2`, it differentiates the Qs normalization and propagates the
  differentiated differential equation alongside the function.
- `jacobian_rmn` differentiates the radial spherical-Bessel expansion, including
  its coefficient normalization. Near boundaries it propagates the differentiated
  radial equation from a converged expansion at `x=2`.

For example:

```julia
jacobian_eigen(1, 2, 1.25; precision=:quad)
jacobian_smn(1, 2, 1.25, [0.3]; precision=:quad)
jacobian_smn(1, 2, 1.25, [0.3]; precision=:quad, kind=2)
jacobian_rmn(1, 2, 1.25, [2.0]; precision=:quad, kind=2)
```

Supplying `h` explicitly selects the original centered-difference calculation
for any of these functions. `adaptive=false` then returns the derivative at
that step and reports its consistency with a halved step.
`adaptive`, `rtol`, and `atol` control finite-difference step checks; coefficient
refinement follows the requested `precision`.

Real first-kind derivatives with `abs(c) > n + 1` use the same extra working precision
and tighter coefficient tolerances as `smn`, protecting small interior values
against cancellation. Both angular kinds support real and complex parameters
in double and quad precision.

For `kind=2`, `normalize=true` remains unsupported. Exact singular endpoints
`eta=±1` return `NaN` for both parameter derivatives and report nonfinite metadata;
this does not assert a limit obtained by differentiating only the leading
endpoint singularity. Interior coordinates, including `c=0`, are supported.
At `c=0` the interior parameter derivatives vanish because the functions are
even in `c`. Singular or unresolved Qs normalization raises an error, as for `smn`.

For Qs, metadata uses `method=:differentiated_equation`; its coefficient residuals
describe the initial-data calculation, not an error bound for the propagated
derivative. The differentiated equation is

```math
L_c[\partial_c Q]=-(\lambda_c-2\sigma c\eta^2)Q,
\qquad L_c[y]=(1-\eta^2)y''-2\eta y'
+\left(\lambda-\sigma c^2\eta^2-\frac{m^2}{1-\eta^2}\right)y,
```

with `sigma=1` for prolate and `sigma=-1` for oblate. Initial data also
differentiate the coefficient sums in the [DLMF Wronskian normalization](https://dlmf.nist.gov/30.5#E4).

Radial parameter derivatives support both geometries and all four kinds.
Kinds 3 and 4 differentiate their direct spherical Hankel expansions and
propagate those waves and their sensitivities together. The identities
`R3_c = R1_c + im*R2_c` and `R4_c = R1_c - im*R2_c` still hold; forming them
from separately rounded values can lose the decaying wave. The default method has no step-size
parameter; supplying `h` explicitly retains centered differences for comparison.
The initial data follow the [radial Bessel expansion](https://dlmf.nist.gov/30.11#E3).
For propagation, the differentiated equation is

```math
L_c[\partial_c R]=-(2cx^2-\lambda_c)R,
\qquad L_c[y]=(x^2-\sigma)y''+2xy'
+\left(c^2x^2-\lambda-\frac{\sigma m^2}{x^2-\sigma}\right)y.
```

Here `sigma=1` for prolate and `sigma=-1` for oblate. The oblate origin is
regular and supported. At the prolate boundary `x=1`, real first-kind results
use one-sided limits (the mixed derivative can be infinite); kinds 2–4 return
`NaN`. Complex prolate coordinates still require `x>1`. All radial kinds reject
exactly zero `c`.

Radial metadata uses `method=:differentiated_expansion` or
`:differentiated_equation`, with `step_used=nothing`. `radial_series_tail` checks
the initial Bessel sum; coefficient residuals and this tail are convergence
diagnostics, not error bounds for the complete propagated result.

For complex `c = a + ib`, eigenvalue and angular derivatives use the local
analytic branch selected by the wave-function evaluator: the `b` partial is
`im` times the `a` partial. The coefficient normalization uses a transpose
bilinear product, without complex conjugation. Radial finite-difference stencils
transport the same mode locally. Genuine eigenvalue branch points can make
sensitivities singular or prevent resolution of a mode.

The coefficient recurrence is based on [DLMF 30.8](https://dlmf.nist.gov/30.8),
expressed in an orthonormal Ferrers basis and the package's eigenvalue convention.
For `H(c)v = lambda*v`, differentiating gives

```math
\lambda' = v^T H'v,\qquad
(H-\lambda I)v' = -(H'-\lambda'I)v,\qquad v^Tv=1,\quad v^Tv'=0.
```

At first-kind angular endpoints, mixed coordinate/parameter derivatives use one-sided
limits; these can be infinite even when the parameter derivative of the value
is zero.

### Reliability metadata

With `with_metadata=true`, `method` identifies the calculation:

- `:coefficients`: expansion refinement, coefficient tail, and eigenproblem and
  sensitivity residuals. `step_used` and `relative_change_when_halving_step`
  are `nothing`, because no parameter step was taken.
- `:differentiated_expansion` or `:differentiated_equation`: coefficient
  diagnostics and, for radial calculations, the Bessel-series tail.
  `step_used` and `relative_change_when_halving_step` are `nothing`.
- `:finite_difference`: `step_used` and
  `relative_change_when_halving_step` describe the numerical step comparison.

All methods report `finite_flag`, `conditioning_flag`, and `suggested_action`.
These are numerical consistency diagnostics, not certified error bounds.
Explicit finite differences can lose digits through subtraction even when
the underlying function values are accurate.

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

For the complete long-form theory-and-usage narrative, see [docs/mathematical-basis-and-usage.md](https://github.com/brandynlucca/SpheroidalWaves.jl/blob/main/docs/mathematical-basis-and-usage.md) in the repository.

