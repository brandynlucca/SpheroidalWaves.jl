# Mathematical basis and usage

Use integer order and degree `0≤m≤n`.
The parameter `c` and evaluation coordinates may be real or complex. Write `σ=1` for prolate and `σ=-1` for oblate geometry.

## [Precision and outputs](@id precision-and-outputs)

~~~julia
s = smn(1, 2, 1.25, [-0.3, 0.3]; precision=:quad)
r = rmn(1, 2, 1.25, 2.0; precision=:quad)
s.value
s.derivative
~~~

`value` and `derivative` are vectors; scalar coordinates give length-one vectors.
The derivative is with respect to the coordinate.

| Precision | Real angular values and eigenvalues | Complex values and radial outputs |
|:--|:--|:--|
| `:double` (default) | `Float64` | `ComplexF64` |
| `:quad` | `BigFloat` | `Complex{BigFloat}` |

Input precision is what the arguments retain, calculation precision is the arithmetic used, output precision is what returned numbers store, and accuracy is agreement with the mathematical result.

`precision=:quad` uses fixed native quad arithmetic (113 significand bits, about 34 decimal digits). The caller's `BigFloat` precision and rounding control output components and some Julia calculations. Increasing them does not increase native precision or guarantee more accurate digits. `scaled=true` extends range without adding accuracy.

| Function | Effect of `BigFloat` settings with `precision=:quad` |
|:--|:--|
| `smn`, `rmn` | Value and derivative output rounding, plus working precision in guarded calculations. |
| `spheroidal` | Follows the requested individual functions. |
| `eigenvalue` | Output rounding and guarded integral-eigenvalue refinement. Native separation calculations stay fixed quad. |
| `eigenvalue_sweep` | Parameter conversion, stored eigenvalues and predictor arithmetic. |
| `dmn`, `kmn`, `amn` | Output rounding and guarded working precision. Expansion tolerances remain separate. |
| `jacobian_eigen`, `jacobian_smn`, `jacobian_rmn` | Output rounding and Julia derivative arithmetic. Finite-difference steps and tolerances remain separate. |
| `find_c_for_eigenvalue` | Target and bracket conversion, root arithmetic and outputs. Default stopping tolerances stay unchanged. |
| `radial_wronskian` | Output rounding and arithmetic combining radial results. |
| `accuracy` | Integer estimates remain integers, with `-1` for unavailable estimates. Refinement and diagnostics can depend on the settings. |

Construct decimal inputs from strings within the precision scope. Widening `0.65` after Julia has rounded it to `Float64` cannot recover lost digits:

~~~julia
setprecision(BigFloat, 256) do
    c, eta = parse.(BigFloat, ("0.65", "0.3"))
    s = smn(2, 3, c, eta; precision=:quad)
    # Compare using exactly representable inputs.
    d = only(rmn(0, 1, 1.0, 2.0; precision=:double).value)
    q = only(rmn(0, 1, 1.0, 2.0; precision=:quad).value)
    (output_bits=precision(real(q)), difference=abs(Complex{BigFloat}(d) - q))
end
~~~

The difference measures agreement, not an accuracy guarantee. `accuracy` provides value estimates, not rigorous bounds or derivative certification.

A future arbitrary-precision request could specify a bit count explicitly while preserving `:double` defaults and existing `:quad` behavior. It would require refinement beyond native quad accuracy. This capability is not currently available, and arbitrary precision alone guarantees neither rigorous bounds nor correct rounding.

## Parallel evaluation

Calls are thread-safe. Start Julia with `julia --threads=auto --project`:

~~~julia
parameters = [0.5, 1.0, 1.5, 2.0]
tasks = [Threads.@spawn smn(1, 2, c, [-0.3, 0.3]; precision=:quad)
         for c in parameters]
results = fetch.(tasks)
~~~

Keep shared inputs unchanged and use separate output buffers for each task.

## Angular functions

`smn(m,n,c,eta)` solves

~~~math
\frac{d}{d\eta}\left[(1-\eta^2)S'\right]
+\left(\lambda_{mn}(c)-\sigma c^2\eta^2-\frac{m^2}{1-\eta^2}\right)S=0,
\qquad -1<\eta<1.
~~~

The separation constant satisfies `λₘₙ(0)=n(n+1)`.
The DLMF parameter is `γ²=σc²`;
see [DLMF 30.2](https://dlmf.nist.gov/30.2).

### First kind

`kind=1` is the default. It uses Meixner–Schäfke normalization and includes
the Condon–Shortley phase `(-1)^m`:

~~~math
S_{mn}(0,\eta)=\mathsf P_n^m(\eta),\qquad
S_{mn}(c,-\eta)=(-1)^{n-m}S_{mn}(c,\eta).
~~~

For real `c` and fixed `m`:

~~~math
\int_{-1}^{1}S_{mn}(c,\eta)S_{mk}(c,\eta)\,d\eta
=\frac{2(n+m)!}{(2n+1)(n-m)!}\,\delta_{nk}.
~~~

`normalize=true` multiplies by the reciprocal square root of this norm,
giving unit norm. For complex `c`, normalization is continued analytically;
it is not a unit integral of the squared magnitude.
See [DLMF 30.4](https://dlmf.nist.gov/30.4).

~~~julia
smn(1, 1, 0.0, 0.3).value  # approximately [-0.953939201416946]
smn(1, 2, 1.25, [-0.3, 0.3]; normalize=true, precision=:quad)
~~~

### First-kind endpoints

`smn` accepts $\eta=\pm1$ and returns the one-sided value and coordinate-derivative limits from inside $[-1,1]$. This applies to both geometries and supported real or complex parameters, including $c=0$.

The regular angular solution factors as $S_{mn}(c,\eta)=(1-\eta^2)^{m/2}F(\eta)$, where $F$ is analytic at the endpoints. See [DLMF 30.3(i)](https://dlmf.nist.gov/30.3#i). The first derivative therefore has the following limits:

| Order | Value at $\eta=\pm1$ | First derivative |
|:--|:--|:--|
| $m=0$ | Finite | Finite |
| $m=1$ | Zero | Divergent |
| $m=2$ | Zero | Finite, generally nonzero |
| $m>2$ | Zero | Zero |

For $m=2$, differentiation gives $S'_{2n}(c,1)=-2F(1)$ and $S'_{2n}(c,-1)=2F(-1)$. For example, `smn(2,3,0.65,1)` returns a zero value and a derivative of approximately $-29.53$. At $c=0$, $S_{23}(0,\eta)=15\eta(1-\eta^2)$ has the endpoint derivative $-30$. The singular coefficients of the differential equation do not require the
regular solution's derivative to diverge. Other evaluators may omit these limits or return zero at an endpoint where the slope is nonzero.

Higher coordinate derivatives have their own finite or divergent limits. `accuracy(...; target=:angular)` returns `-1` at endpoints because no calibrated digit estimate is available there.

### Second kind

`kind=2` selects the Ferrers-branch `Qs`:

~~~julia
q = smn(1, 2, 1.25, [-0.3, 0.3]; kind=2, precision=:quad)
~~~

Its parity, zero-parameter limit, and normalization are

~~~math
Q(c,-\eta)=(-1)^{n-m+1}Q(c,\eta),\qquad
Q(0,\eta)=\mathsf Q_n^m(\eta),
~~~

~~~math
(1-\eta^2)(P Q'-P'Q)
=\frac{(n+m)!}{(n-m)!}A_n^m A_n^{-m}.
~~~

Here `P` is the default first-kind function, and `A` denotes the joining sums
in [DLMF 30.11.4](https://dlmf.nist.gov/30.11#E4).
The second-kind convention is specified by
[DLMF 30.5](https://dlmf.nist.gov/30.5) and [30.8(ii)](https://dlmf.nist.gov/30.8#ii).

!!! warning "Second-kind normalization and endpoints"
    `normalize=true` is unsupported for `kind=2`.
    Values and coordinate derivatives diverge at `η=±1` and return signed one-sided
    infinities. Singular or numerically unresolved normalization raises an error.

## Radial functions

`rmn(m,n,c,x)` solves

~~~math
\frac{d}{dx}\left[(x^2-\sigma)R'\right]
+\left(c^2x^2-\lambda_{mn}(c)-\frac{\sigma m^2}{x^2-\sigma}\right)R=0.
~~~

| `kind` | Function |
|:--|:--|
| `1` (default) | First kind `R₁` |
| `2` | Second kind `R₂` |
| `3` | `R₁+iR₂` |
| `4` | `R₁-iR₂` |

~~~julia
rmn(1, 2, 1.25, [1.5, 2.0]; kind=2, precision=:quad)
rmn(1, 2, 1.25, [0.0, 1.0]; spheroid=:oblate, kind=3, precision=:quad)
~~~

!!! warning "Radial domain"
    Real prolate coordinates satisfy `x≥1`, real oblate coordinates `x≥0`.
    At prolate `x=1`, the first kind uses one-sided limits. Kinds 2–4 return `NaN`.
    Standard normalization requires `c≠0`.

`normalization=:static` evaluates $c^{-n}R_1$ or $c^{n+1}R_2$, including their finite limits at `c=0`. It supports kinds 1 and 2 in `rmn`, `jacobian_rmn`, `radial_wronskian`, and radial `accuracy`.

~~~julia
rmn(0, 1, 0, 2; normalization=:static) # value ≈ [2/3]
~~~

## Complex branches

Complex coordinates continue the real-coordinate normalization on these cuts:

| Function | Excluded coordinates |
|:--|:--|
| Angular | Real rays `z < -1` and `z > 1` |
| Prolate radial | Real segment `-1 ≤ z < 1` |
| Oblate radial | Imaginary rays `z=iy`, `abs(y) ≥ 1` |

Angular `z=±1` and prolate `z=1` retain their one-sided endpoint conventions. Oblate continuation passes through `z=0`. Coordinate order does not change the branch. Complex coordinates support the same derivative, scaling, and degree-range options. Use `accuracy(...; diagnostics=true)` for convergence and conditioning estimates.

At `c=0`, exact-reference checks cover angular zeros and stationary points, scaled angular values at `m=n=200`, and both sides of these cuts at distances `2^-12` and `2^-40`. Values and derivatives are checked separately, using absolute errors near zeros. These samples do not establish accuracy throughout the domain or for nonzero `c`.

~~~julia
smn(1, 2, 1.25, 0.3 + 0.2im; derivatives=2)
rmn(1, 2, 1.25, [1.5 + 0.2im, 2.0 - 0.1im]; kind=3)
~~~

For complex `c`, scalar `eigenvalue`, `smn`, and `rmn` follow degree
`n` from `real(c)` along the vertical segment to `c`.
Angular normalization and phase follow that branch.

!!! warning "Branch points"
    Eigenmodes can exchange around a complex branch point.
    Use `eigenvalue_sweep` to specify another continuation path.
    Unresolved paths raise an error.

## Further calculations

[Mathematical tools](mathematical-tools.md) covers eigenvalue operators,
parameter derivatives, inverse bandwidth, sweeps, degree ranges, scaled values,
logarithmic derivatives, higher coordinate derivatives, and accuracy diagnostics.
