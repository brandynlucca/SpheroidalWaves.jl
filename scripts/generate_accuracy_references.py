"""Regenerate independent fixtures with Python and mpmath 1.3.0.

Without arguments, print the angular zero fixtures. For the wider study:
python scripts/generate_accuracy_references.py --matrix test/accuracy_references.tsv
Reference generation is separate from the Julia test suite.
"""

import argparse
import csv
from decimal import Decimal, localcontext
import mpmath as mp


def legendre_mode(m, n, c, sigma, terms):
    """Return an eigenpair and ordinary Ferrers expansion coefficients.

    The angular equation gives H = diag(l(l+1)) + sigma*c^2*x^2,
    sigma=+1 for prolate and -1 for oblate. Recurrences and normalization:
    https://dlmf.nist.gov/14.10.E3, https://dlmf.nist.gov/14.10.E5,
    https://dlmf.nist.gov/30.4.E1.
    """
    first = m + (n - m) % 2
    degrees = list(range(first, first + 2 * terms, 2))
    matrix = mp.matrix(terms)

    def a(l):
        return mp.sqrt(mp.mpf((l + 1)**2 - m*m) / ((2*l + 1)*(2*l + 3))) if l >= m else mp.mpf(0)

    for j, l in enumerate(degrees):
        matrix[j, j] = l * (l + 1) + sigma*c*c * (a(l) ** 2 + a(l - 1) ** 2)
        if j + 1 < terms:
            matrix[j, j + 1] = matrix[j + 1, j] = sigma*c*c * a(l) * a(l + 1)
    if mp.im(c):
        eigenvalues, eigenvectors = mp.eig(matrix, left=False, right=True)
        # These mildly complex parameters retain the real-axis ordering.
        order = sorted(range(terms), key=lambda j: mp.re(eigenvalues[j]))
    else:
        eigenvalues, eigenvectors = mp.eigsy(matrix)
        order = list(range(terms))
    index = order[(n - m) // 2]
    vector = [eigenvectors[j, index] for j in range(terms)]
    norm2 = lambda l: 2 * mp.factorial(l + m) / ((2*l + 1) * mp.factorial(l - m))
    scale = mp.sqrt(norm2(n) / mp.fsum(v*v for v in vector))
    if mp.re(scale * vector[(n - m) // 2]) < 0:
        scale = -scale
    coefficients = [scale*v/mp.sqrt(norm2(l)) for v, l in zip(vector, degrees)]
    return eigenvalues[index], degrees, coefficients


def angular_mode(n, sigma, terms):
    eigenvalue, degrees, coefficients = legendre_mode(0, n, mp.mpf(1), sigma, terms)

    def evaluate(x):
        values, slopes = [], []
        for l, coefficient in zip(degrees, coefficients):
            p = mp.legendre(l, x)
            dp = l * (mp.legendre(l - 1, x) - x * p) / (1 - x * x) if l else mp.mpf(0)
            values.append(coefficient * p)
            slopes.append(coefficient * dp)
        return mp.fsum(values), mp.fsum(slopes)

    assert mp.sign(evaluate(mp.mpf(0))[n % 2]) == (-1) ** (n // 2)
    return eigenvalue, evaluate


def zero_fixtures():
    print(f"# mpmath {mp.__version__}, m=0, c=1, default angular normalization")
    print("# Refinements: (70 decimal digits, 24 terms), (100, 24), (100, 32)")
    for geometry, sigma in (("prolate", 1), ("oblate", -1)):
        for n in (2, 3):
            with mp.workdps(100):
                _, evaluate = angular_mode(n, sigma, 32)
                root = mp.findroot(
                    lambda x: evaluate(x)[n - 2], (mp.mpf("0.4"), mp.mpf("0.7"))
                )
                numerator = int(mp.floor(root * 2**32))
            snapshots = []
            for digits, terms in ((70, 24), (100, 24), (100, 32)):
                with mp.workdps(digits):
                    eigenvalue, evaluate = angular_mode(n, sigma, terms)
                    points = [mp.mpf(k) / 2**32 for k in (numerator, numerator + 1)]
                    values = [evaluate(x) for x in points]
                    assert values[0][n - 2] * values[1][n - 2] < 0
                    snapshots.append([eigenvalue, *values[0], *values[1]])
            with mp.workdps(100):
                changes = [
                    max(abs(x - y) for x, y in zip(a, b))
                    for a, b in zip(snapshots, snapshots[1:])
                ]
                assert max(changes) < mp.mpf("1e-55"), changes
                print(f"# precision/truncation absolute changes: {', '.join(mp.nstr(x, 5) for x in changes)}")
                print(f"(:{geometry}, {n}, {numerator},")
                for value in snapshots[-1]:
                    print(f'    big"{mp.nstr(value, 50)}",')
                print("),")


def polynomial_mode(m, n, c, sigma, terms):
    """Rodrigues polynomials times the explicit Ferrers square-root factor."""
    eigenvalue, degrees, coefficients = legendre_mode(m, n, c, sigma, terms)
    powers = [mp.mpf(0)] * (degrees[-1] - m + 1)
    for l, d in zip(degrees, coefficients):
        for k in range((l - m) // 2 + 1):
            exponent = l - 2*k - m
            powers[exponent] += d * (-1)**(m + k) * mp.factorial(2*l - 2*k) / (
                2**l * mp.factorial(k) * mp.factorial(l - k) * mp.factorial(exponent)
            )
    def evaluate(z):
        p, dp = mp.polyval(list(reversed(powers)), z, derivative=True)
        if not m:
            return p, dp
        root = mp.sqrt(1 - z) * mp.sqrt(1 + z)
        return root**m*p, root**m*dp - m*z*root**(m - 2)*p
    return eigenvalue, degrees, coefficients, evaluate


def integrate_path(function, points):
    return mp.fsum(mp.quad(lambda t: (b-a)*function(a+(b-a)*t), [0, 1])
                   for a, b in zip(points, points[1:]) if a != b)


def second_angular(coefficients, regular, z):
    # m=n=0: DLMF 30.5.2 and 30.5.4 fix Q(0)=0 and its Wronskian.
    # Reduction of order with contour quadrature is independent of the
    # package's coordinate Taylor propagation. Detours fix the cut side.
    height = mp.mpf("0.5") * (1 if mp.im(z) >= 0 else -1)
    path = [0, 1j*height, mp.re(z)+1j*height, z] if abs(mp.re(z)) >= 1 else [0, z]
    constant = mp.fsum(coefficients)**2
    integral = integrate_path(lambda t: constant/((1-t*t)*regular(t)[0]**2), path)
    p, dp = regular(z)
    return p*integral, dp*integral + constant/((1-z*z)*p)


def radial_mode(c, sigma, terms, n=0):
    _, degrees, coefficients, angular = polynomial_mode(0, n, c, sigma, terms)
    anchor = mp.mpf(4)
    argument = c*anchor
    denominator = mp.fsum(coefficients)
    def bessel(kind):
        f = mp.besselj if kind == 1 else mp.bessely
        return mp.fsum((-1)**((l-n)//2)*d*mp.sqrt(mp.pi/(2*argument))*f(l+mp.mpf("0.5"), argument)
                       for l, d in zip(degrees, coefficients))/denominator
    # DLMF 30.11.3 fixes radial normalization at a convergent exterior anchor.
    factor = 1 if sigma == 1 else -1j
    connection = bessel(1)/angular(factor*anchor)[0]
    def regular(z):
        value, slope = angular(factor*z)
        return connection*value, connection*factor*slope
    initial_ratio = bessel(2)/regular(anchor)[0]
    def second(z):
        height = mp.mpf("0.5") * (1 if mp.im(z) >= 0 else -1)
        if sigma == -1 and abs(mp.im(z)) >= 1:
            side = mp.mpf("0.5") * (1 if mp.re(z) >= 0 else -1)
            path = [anchor, anchor+1j*height, side+1j*height, side+1j*mp.im(z), z]
        else:
            path = [anchor, anchor+1j*height, mp.re(z)+1j*height, z]
        integral = integrate_path(lambda t: 1/(c*(t*t-sigma)*regular(t)[0]**2), path)
        value, slope = regular(z)
        ratio = initial_ratio + integral
        return value*ratio, slope*ratio + 1/(c*(z*z-sigma)*value)
    return regular, second


def matrix_cases():
    """Dyadic inputs cover real/complex parameters and coordinates separately."""
    for geometry in ("prolate", "oblate"):
        for ci in ("0", "0.25"):
            for family in ("angular1", "angular2", "radial1", "radial2"):
                m = n = 1 if family == "angular1" else 0
                points = [("0.25" if family.startswith("angular") else "2", xi, "interior")
                          for xi in ("0", "0.125")]
                for exponent in (12, 40):
                    for side in (-1, 1):
                        with localcontext() as context:
                            context.prec = 80
                            delta = str(Decimal(side) / Decimal(2)**exponent)
                        if family.startswith("angular"):
                            points.append(("1.5", delta, "cut"))
                        elif geometry == "prolate":
                            points.append(("0.5", delta, "cut"))
                        else:
                            points.append((delta, "1.5", "cut"))
                for xr, xi, region in points:
                    yield (family, geometry, m, n, "1", ci, xr, xi, region)
        for m, n, c, region in ((0, 0, "0.000244140625", "small"),
                                (50, 50, "1", "high-order"), (0, 60, "1", "high-degree")):
            yield ("angular1", geometry, m, n, c, "0", "0.25", "0", region)
        for family in ("radial1", "radial2"):
            for xi in ("0", "0.125"):
                yield (family, geometry, 0, 1, "1", "0.25", "2", xi, "normalization")


def matrix_snapshot(cases, digits, terms):
    result, modes = [], {}
    with mp.workdps(digits):
        for family, geometry, m, n, cr, ci, xr, xi, _ in cases:
            c = mp.mpf(cr) if ci == "0" else mp.mpc(cr, ci)
            z = mp.mpc(xr, xi)
            sigma = 1 if geometry == "prolate" else -1
            key = (family.startswith("radial"), geometry, m, n, cr, ci)
            if key not in modes:
                modes[key] = radial_mode(c, sigma, terms+(n-m)//2, n=n) if key[0] else polynomial_mode(m, n, c, sigma, terms+(n-m)//2)
            if key[0]:
                values = modes[key][int(family[-1])-1](z)
            else:
                _, _, coefficients, regular = modes[key]
                values = regular(z) if family == "angular1" else second_angular(coefficients, regular, z)
            result.append(values)
            print(f"{digits} digits/{terms} terms: {family} {geometry} c={c} z={z}", flush=True)
    return result


def write_matrix(path, region=None):
    cases = [case for case in matrix_cases() if region is None or case[-1] == region]
    snapshots = [matrix_snapshot(cases, digits, terms) for digits, terms in ((90, 32), (120, 32), (120, 44))]
    with mp.workdps(120), open(path, "w", newline="", encoding="utf-8") as output:
        output.write(f"# mpmath {mp.__version__}, independent Legendre eigensolves and contour quadrature\n")
        output.write("# Refinements: (90 decimal digits, 32 terms), (120, 32), (120, 44)\n")
        output.write("# Term counts additionally include the mode index (n-m)//2.\n")
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("family", "geometry", "m", "n", "c_re", "c_im", "x_re", "x_im", "region",
                         "value_re", "value_im", "slope_re", "slope_im", "precision_change", "truncation_change"))
        for case, low, precise, refined in zip(cases, *snapshots):
            changes = [max(abs(a-b)/max(abs(b), mp.mpf("1e-50")) for a, b in zip(left, right))
                       for left, right in ((low, precise), (precise, refined))]
            assert max(changes) < mp.mpf("1e-34"), (case, changes)
            components = [part(value) for value in refined for part in (mp.re, mp.im)]
            writer.writerow((*case, *(mp.nstr(x, 45) for x in components), *(mp.nstr(x, 8) for x in changes)))
    print(f"Wrote {len(cases)} converged reference pairs to {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", metavar="PATH", help="write the wider reference matrix as TSV")
    parser.add_argument("--region", choices=("interior", "cut", "small", "high-order", "high-degree", "normalization"),
                        help="generate only one region of the matrix")
    args = parser.parse_args()
    if args.region and not args.matrix:
        parser.error("--region requires --matrix")
    write_matrix(args.matrix, args.region) if args.matrix else zero_fixtures()
