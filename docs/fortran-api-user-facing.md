# Spheroidal Wave Function Fortran API Reference

## Scope

The four spheroidal Fortran source files (prolate, oblate, complex prolate, complex oblate) each provide a main public subroutine that computes radial and angular spheroidal wave functions and first derivatives. These subroutines include both user-meaningful arguments and internal computational parameters.

This document defines a stable user-facing argument model for the low-level Fortran interface and distinguishes it from backend-internal workspace controls.

## Position in the Documentation Set

This reference is part of the backend-interface documentation.

- Primary high-level documentation:
  - [README.md](../README.md)
  - [docs/mathematical-basis-and-usage.md](mathematical-basis-and-usage.md)
- Primary low-level audience:
  - backend contributors
  - maintainers of custom bindings that call Fortran entry points directly

This file is the canonical mapping between Fortran subroutine signatures and user-meaningful parameter semantics.

---

## Public Entry Subroutine Signatures

| File | Main Subroutine | Geometry | Eigenvalue Type |
|------|-----------------|----------|-----------------|
| prolate_swf.f90 | `profcn()` | Prolate spheroid | Real eigenvalues |
| oblate_swf.f90 | `oblfcn()` | Oblate spheroid | Real eigenvalues |
| complex_prolate_swf.f90 | `cprofcn()` | Prolate spheroid | Complex eigenvalues |
| complex_oblate_swf.f90 | `coblfcn()` | Oblate spheroid | Complex eigenvalues |

---

## User-Facing Input Parameters

All four subroutines share the same **core input parameter set**:

| Argument | Type | Semantics | Units / Range | Example |
|----------|------|-----------|----------------|---------|
| `c` | `real(knd)` | **Size parameter** = $c = ka/2 = \pi \cdot d / \lambda$ where $d$ is interfocal distance. Also written as $kd/2$ in literature. | Positive real | 1.0, 10.0, 100.0 |
| `m` | `integer` | **Order** (azimuthal mode number) | $m \geq 0$ | 0, 1, 2, 5 |
| `lnum` | `integer` | **Degree range**: angular and radial functions are computed for $l = m, m+1, \ldots, m + \text{lnum} - 1$ | $\geq 1$ | 10, 50, 100 |
| `ioprad` | `integer` | **Radial function control flag** | 0 = none; 1 = first kind only; 2 = both kinds | 0, 1, 2 |
| `x1` | `real(knd)` | **Radial coordinate** (prolate: $x_1 = x - 1$ where $x \geq 1$; oblate: $x_1 = x - 1$ where $x \leq -1$) | $x_1 \geq 0$ (prolate); $x_1 \geq 0$ (oblate) | 0.5, 2.0, 10.0 |
| `iopang` | `integer` | **Angular function control flag** | 0 = none; 1 = first kind only; 2 = first kind + derivatives | 0, 1, 2 |
| `iopnorm` | `integer` | **Angular function normalization** | 0 = Meixner-Schäfke norm; 1 = unity norm (scaled) | 0, 1 |
| `narg` | `integer` | **Number of eta values** (angular coordinate samples) | $\geq 1$ | 1, 10, 100 |
| `arg(narg)` | `real(knd)` array | **Angular coordinate values** ($\eta$ in standard notation; physically: $\eta = \cos(\theta)$ where $\theta$ is the polar angle) | $[-1, 1]$ for angular; $[0, 1]$ for prolate radial; $[-1, 0]$ for oblate radial | [0.0], [0.0, 0.5, 1.0] |

### Notes on Input Parameters

- **`c` (size parameter)**: Controls the regime of computation. Small $c$ ($\ll 1$) favors series expansions; large $c$ ($\gg 1$) requires asymptotic methods and higher precision.
- **`m` (order)**: Azimuthal wavenumber in cylindrical coordinates. Spheroidal harmonics are labeled $S_{mn}^{(1/2)}$ (angular) and $R_{mn}^{(1/2)}$ (radial).
- **`lnum` (degree range)**: Each call evaluates `lnum` consecutive degree values. Typical usage: `lnum = 10` to `100`.
- **`ioprad` flag**: Common patterns:
  - `ioprad = 0`: Skip radial functions entirely (compute only angular).
  - `ioprad = 1`: Compute only first-kind radial functions (interior problems).
  - `ioprad = 2`: Compute both kinds (exterior/scattering problems).
- **`x1` (radial coordinate)**: For **prolate**, $x = x_1 + 1$ must satisfy $x \geq 1$ (equivalently $x_1 \geq 0$). For **oblate**, $x = x_1 + 1$ must satisfy $x \leq -1$ (equivalently $x_1 \geq 0$). The subroutine automatically computes $x$ internally.
- **`iopang` flag**: Common patterns:
  - `iopang = 0`: No angular functions (pure radial).
  - `iopang = 1`: Angular functions only (no derivatives).
  - `iopang = 2`: Angular functions and their first derivatives w.r.t. $\eta$.
- **`iopnorm` flag**: Affects normalization of angular functions:
  - `iopnorm = 0`: Uses Meixner-Schäfke normalization (matches associated Legendre functions at boundary).
  - `iopnorm = 1`: Scales angular functions to have unity norm (numerically cleaner, no overflow risk).
- **`arg(narg)` array**: For angular functions, typically `arg = [0, 0.1, 0.5, 1.0]` to sample the polar angle $\theta$. For radial functions in certain methods, can contain auxiliary coordinate values.

---

## User-Facing Output Parameters

### Radial Functions (when `ioprad /= 0`)

| Argument | Type | Semantics | Shape | Notes |
|----------|------|-----------|-------|-------|
| `r1c(lnum)` | `real(knd)` | Characteristics of first-kind radial functions $R_{1}^{(1)}$ | 1D, length `lnum` | Mantissa; multiply by `10^(ir1e)` to get full value. |
| `r1dc(lnum)` | `real(knd)` | Characteristics of derivatives $\frac{dR_{1}^{(1)}}{dx}$ | 1D, length `lnum` | Mantissa form. |
| `ir1e(lnum)` | `integer` | Exponents for `r1c` | 1D, length `lnum` | Power of 10 scaling. |
| `ir1de(lnum)` | `integer` | Exponents for `r1dc` | 1D, length `lnum` | Power of 10 scaling. |
| `r2c(lnum)` | `real(knd)` | Characteristics of second-kind radial functions $R_{2}^{(1)}$ *(if `ioprad == 2`)* | 1D, length `lnum` | Only present when `ioprad == 2`. |
| `r2dc(lnum)` | `real(knd)` | Characteristics of derivatives $\frac{dR_{2}^{(1)}}{dx}$ *(if `ioprad == 2`)* | 1D, length `lnum` | Only present when `ioprad == 2`. |
| `ir2e(lnum)` | `integer` | Exponents for `r2c` | 1D, length `lnum` | Only present when `ioprad == 2`. |
| `ir2de(lnum)` | `integer` | Exponents for `r2dc` | 1D, length `lnum` | Only present when `ioprad == 2`. |
| `naccr(lnum)` | `integer` | Estimated accuracy (decimal digits) of radial functions | 1D, length `lnum` | Typical range: 6–15. Use to assess solution quality. |

### Angular Functions (when `iopang /= 0`)

| Argument | Type | Semantics | Shape | Notes |
|----------|------|-----------|-------|-------|
| `s1c(lnum, narg)` | `real(knd)` | Characteristics of angular functions $S_{1}^{(1)}$ | 2D, shape `(lnum, narg)` | Mantissa; multiply by `10^(is1e)` to get full value. `s1c(i, j)` corresponds to degree $l = m + i - 1$ and $\eta = \text{arg}(j)$. |
| `s1dc(lnum, narg)` | `real(knd)` | Characteristics of derivatives $\frac{dS_{1}^{(1)}}{d\eta}$ | 2D, shape `(lnum, narg)` | Only if `iopang == 2`. |
| `is1e(lnum, narg)` | `integer` | Exponents for `s1c` | 2D, shape `(lnum, narg)` | Power of 10 scaling. |
| `is1de(lnum, narg)` | `integer` | Exponents for `s1dc` | 2D, shape `(lnum, narg)` | Only if `iopang == 2`. |
| `naccs(lnum, narg)` | `integer` | Estimated accuracy (decimal digits) of angular functions | 2D, shape `(lnum, narg)` | Typical range: 6–15. Use to assess solution quality. |

### How to Reconstruct Full Values from Characteristic Form

The output arrays use **characteristic-exponent** (mantissa-exponent) representation to avoid overflow/underflow:

$$\text{Full value} = \text{characteristic} \times 10^{\text{exponent}}$$

Example in Julia:
```julia
# For radial function of first kind at degree l = m + i - 1:
r1_full = r1c[i] * 10.0^ir1e[i]

# For angular function at degree l = m + i - 1, angle arg[j]:
s1_full = s1c[i, j] * 10.0^is1e[i, j]
```

---

## Internal Computational Parameters (NOT for users)

The following arguments are automatically computed and should **never** be exposed to end-users. They are intermediate workspaces and algorithm-control parameters:

### Workspace Dimensions
- `maxd`, `maxdr`, `maxint`, `maxj`, `maxlp`, `maxm`, `maxmp`, `maxn`, `maxp`, `maxpdr`, `maxq`, `maxt`: Array dimension limits computed automatically from `c`, `m`, `lnum`, and required precision.
- `neta`, `ngau`, `jnenmax`: Quadrature and iteration parameters.

### Precision Metadata
- `ndec`: Number of decimal digits available in current working precision.
- `nex`: Maximum exponent in current precision.
- `kindd`, `kindq`: Kind parameters for double/quadruple precision (Fortran internals).

### Workspace Arrays
- `maxp`-sized arrays: Legendre polynomial coefficients (`alpha`, `beta`, `gamma`, `coefa`–`coefe`).
- `maxq`-sized arrays: Legendre function ratios (`qr`, `qdr`).
- Eigenvalue work: `enr`, `bliste`, `gliste` (expansion coefficient storage).
- Quadrature work: `xr`, `wr` (Gauss-Legendre nodes and weights).
- Special function work: `sbesf`, `sbesdf`, `sneun`, `snedf` (Bessel and Neumann ratios).
- Boundary and matching work: `eta`, `wmeta2`, `xbn`, `xln`.

### Algorithm-Control Flags (Internal)
- `mmin`, `minc`, `mnum`: Loop controls for order $m$ iteration.
- `iopleg`, `iopneu`, `iopeta`, `iopint`: Method-selection flags.
- `limcsav`, `limpsv`, `limnsv`: Dimension tracking for series convergence.

---

## Recommended Low-Level Wrapper Interface

The Julia package already exposes high-level entry points (`smn`, `rmn`, `eigenvalue`, `accuracy`, `radial_wronskian`, `jacobian_eigen`, `jacobian_smn`, `jacobian_rmn`).

For non-Julia wrappers or direct backend integrations, the recommended interface is to expose only core scientific inputs and outputs while hiding backend workspace controls.

For clarity and safety, expose **only** the core input/output sets and keep computational workspace parameters internal.

### Minimal User API Call Signature

```python
# Python-style pseudocode
def compute_spheroidal_functions(
    c: float,
    m: int,
    lnum: int,
    ioprad: int = 2,           # 0=skip, 1=first kind, 2=both
    iopang: int = 2,           # 0=skip, 1=first kind, 2=+derivatives
    iopnorm: int = 0,          # 0=Meixner-Schäfke, 1=unity norm
    x1: float = None,          # radial coordinate (x - 1)
    eta_values: array = None   # angular coordinates
):
    """
    Returns named tuples or dicts:
    {
      'radial': {
          'r1': array(lnum),           # first kind
          'r1d': array(lnum),          # first kind derivative
          'r1_exp': array(lnum),       # exponents (power of 10)
          'r1d_exp': array(lnum),
          'r2': array(lnum),           # second kind (if ioprad==2)
          'r2d': array(lnum),
          'r2_exp': array(lnum),
          'r2d_exp': array(lnum),
          'acc': array(lnum)           # accuracy in decimal digits
      },
      'angular': {
          'S1': array(lnum, narg),     # first kind
          'S1d': array(lnum, narg),    # derivatives (if iopang==2)
          'S1_exp': array(lnum, narg),
          'S1d_exp': array(lnum, narg),
          'acc': array(lnum, narg)
      }
    }
    """
```

### Wrapper Design Principles

1. **Do not expose array dimensions** (`max*` parameters). Allocate dimensions internally based on `c`, `m`, `lnum`.
2. **Provide sensible defaults**. Typical production defaults are `ioprad=2` (both radial kinds), `iopang=2` (angular + derivatives), and `iopnorm=0` or `1` by normalization convention.
3. **Provide a utility function for characteristic-exponent reconstruction**:
   ```python
   def reconstruct_from_characteristic(mantissa, exponent):
       return mantissa * 10.0**exponent
   ```
4. **Validate inputs** before calling Fortran. Required checks include:
   - $c > 0$
   - $m \geq 0$
   - $\ell >= 1$
   - $x_1 \geq 0$
   - $-1 \leq \eta \leq 1$
   - Consistency of `ioprad`, `iopang`, `iopnorm` flags.
5. **Return structured output**. Named tuples or dictionary-like structures reduce ambiguity compared with positional return contracts.

---

## Parameter Interaction Summary

| Use Case | Recommended Settings |
|----------|----------------------|
| Interior problem (radial only, $r \to 0$) | `ioprad=1, iopang=0` |
| Interior problem (radial + angular) | `ioprad=1, iopang=2, iopnorm=1` |
| Scattering problem (exterior) | `ioprad=2, iopang=2` |
| Pure eigenvalue / characteristic study | `ioprad=0, iopang=1` |
| Testing with minimal computation | `ioprad=1, iopang=1, lnum=1, narg=1` |
| High-accuracy (quadruple precision) | Recompile with `knd=16` in `param` module |

---

## References

- [Prolate Spheroid Wave Functions Repository](https://github.com/MathieuandSpheroidalWaveFunctions/Prolate_swf)
- [Oblate Spheroid Wave Functions Repository](https://github.com/MathieuandSpheroidalWaveFunctions/Oblate_swf)
- [Complex Prolate Repository](https://github.com/MathieuandSpheroidalWaveFunctions/complex_prolate_swf)
- [Complex Oblate Repository](https://github.com/MathieuandSpheroidalWaveFunctions/complex_oblate_swf)

---

## Implementation Notes

- All arrays are 1-indexed in Fortran. Wrappers that target 0-indexed environments should translate indices at the binding boundary.
- Precision is controlled globally by the `param` module. Double precision (`knd=8`) is standard; quadruple (`knd=16`) requires recompilation.
- The caching infrastructure (pleg_cache, qleg_cache, gauss_cache) is an implementation detail and should never be exposed to users.
