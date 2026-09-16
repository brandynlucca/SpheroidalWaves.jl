# Mathematical basis and usage

Use integer order and degree `0≤m≤n`.
The parameter `c` may be real or complex; coordinates must be real.
Write `σ=1` for prolate and `σ=-1` for oblate geometry.

## Precision and outputs

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

!!! note "Quad precision"
    Ordinary numeric arguments work with `precision=:quad`.
    It cannot recover digits already lost in the inputs or guarantee 30 accurate digits
    in every parameter regime.

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
    Requires `c≠0`; exact zero throws `DomainError`.
    Prolate coordinates satisfy `x≥1` for real parameters and `x>1` for complex parameters.
    Oblate coordinates satisfy `x≥0`.
    At `x=1`, real prolate first-kind calls use one-sided limits; kinds 2–4 return `NaN`.
    The same parameter restriction applies to radial Jacobians, Wronskians, and accuracy diagnostics.

## Complex branches

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
logarithmic derivatives, second coordinate derivatives, and accuracy diagnostics.
