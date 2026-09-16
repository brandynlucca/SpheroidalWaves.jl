# Mathematical tools

Use `precision=:quad` to request quad-precision calculations. The same numeric
arguments and arrays work with either precision setting.

## Multiple degrees

```julia
smn(1, 1:24, 1.25, [-0.7, 0.2, 0.8];
    spheroid=:oblate, precision=:quad, normalize=true)
rmn(1, 1:24, 1.25, [0.4, 1.0, 3.0];
    spheroid=:oblate, precision=:quad, kind=3, scaled=true)
```

Rows correspond to coordinates and columns to degrees. For a scalar coordinate,
a degree-range call returns one-row matrices. Single-degree calls return arrays,
including length-one arrays for scalar coordinates.

## Wronskian error

```julia
x = [1.5, 2.0, 4.0]
radial_wronskian(1, 2, 1.25, x; precision=:quad, form=:error)
```

The default `form=:raw` returns `R1*R2′-R1′*R2`. For real positive `c`,
`form=:normalized` multiplies this by `c*(x^2-1)` (prolate) or `c*(x^2+1)`
(oblate), so the target is one. `form=:error` returns the absolute difference
from one. The prolate identity is [DLMF 30.11.7](https://dlmf.nist.gov/30.11#E7).
The oblate factor follows the corresponding radial equation and normalization.

A small Wronskian error checks consistency of the solution pair. Correlated
errors can preserve the Wronskian; this is not an independent accuracy bound.


## Scaled values and logarithmic derivatives

```julia
s = smn(200, 200, 0.0, 0.3; scaled=true, logderivative=true)
s.value.mantissa
s.value.exponent
s.logderivative

r = rmn(1, 2, 1.25, 2.0;
    precision=:quad, kind=2, scaled=true, logderivative=true)
```

With `scaled=true`, `value` and `derivative` each contain `mantissa` and integer
`exponent` arrays, representing `mantissa*10^exponent`. Nonzero finite mantissas
have magnitude from one up to ten; zeros and nonfinite values use exponent
zero. Exponents for values and derivatives are independent. Complex values use
one common exponent for the real and imaginary parts. Mantissas use `Float64`
or `BigFloat` (and their complex equivalents), following `precision`.

Use scaled output when values may overflow or underflow the output type.
Scaling extends the representable magnitude, not the number of accurate digits.

`logderivative=true` appends `S′/S` or `R′/R`. At a computed zero of the
function it returns `NaN`; near a zero the ratio is ill-conditioned. Without `scaled=true`,
ordinary values can still overflow while the appended ratio remains finite.

## Second coordinate derivatives

```julia
s = smn(1, 2, 1.25, 0.3;
    precision=:quad, second_derivative=true)
s.second_derivative
```

Both `smn` and `rmn` accept `second_derivative=true`, including with scaled
outputs and degree ranges. The `second_derivative` field contains the second
derivative with respect to the coordinate.

At first-kind angular endpoints and the regular first-kind prolate radial endpoint,
derivatives are one-sided limits. Orders 0, 2, and 4 have finite
second-derivative limits; orders 1 and 3 diverge and return signed infinities;
orders greater than 4 have zero limits. For complex angular functions, divergent
real and imaginary components retain their respective signs. Undefined prolate
second/Hankel kinds at `x=1` return `NaN`.
For angular `kind=2`, values and both coordinate derivatives diverge at either
endpoint; the returned infinities follow their one-sided signs.


