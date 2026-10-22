"""Compare selected SciPy and SpheroidalWaves evaluations with existing fixtures.

Requires NumPy, SciPy, Julia, and the instantiated local Julia project.
Run with --julia /path/to/julia and optionally --output /existing/dir/result.json.
The angular metric S'/S removes amplitude normalization and phase conventions.
These four cancellation stress cases and one radial control are not a domain survey.
"""

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess

import numpy as np
import scipy
from scipy import special


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = next(p for p in Path(__file__).resolve().parents if (p / "Project.toml").is_file())
    angular_path = root / "test" / "oblate_precision.jl"
    radial_path = root / "test" / "quad_precision.jl"
    pattern = r'\((\d+),(\d+),(\d+),\s*big"([^"]+)",big"([^"]+)",\s*big"([^"]+)",\s*big"([^"]+)"\)'
    fixtures = re.findall(pattern, angular_path.read_text(encoding="utf-8"))
    if len(fixtures) != 4:
        raise ValueError("Expected the four documented oblate cancellation fixtures")
    radial_ref = re.search(r'expected = big"([^"]+)"', radial_path.read_text(encoding="utf-8")).group(1)
    cases = [(int(m), int(n), int(c)) for m, n, c, *_ in fixtures]
    cases_julia = ",".join(f"({m},{n},{c})" for m, n, c in cases)
    julia_code = r'''
using SpheroidalWaves
println("META|julia_version|", VERSION)
println("META|package_version|", pkgversion(SpheroidalWaves))
setprecision(BigFloat, 256) do
    for p in (:double, :quad)
        lib = SpheroidalWaves.backend_library(;precision=p)
        lib === nothing && error("Missing $p backend")
        println("META|backend_", p, "|", lib)
        T = p === :quad ? BigFloat : Float64
        for (m,n,c) in CASES
            s = smn(m,n,T(c),T(3)/10;spheroid=:oblate,precision=p)
            println("ANG|", p, "|", m, "|", n, "|", c, "|", s.value[1], "|", s.derivative[1])
        end
        r = rmn(0,1,T(1),T(2);spheroid=:prolate,precision=p,kind=1)
        iszero(imag(r.value[1])) || error("Unexpected imaginary radial component")
        println("RAD|", p, "|", real(r.value[1]))
    end
end
'''.replace("CASES", f"({cases_julia},)")
    completed = subprocess.run(
        [args.julia, f"--project={root}", "--startup-file=no", "-e", julia_code],
        cwd=root, check=True, capture_output=True, text=True, timeout=180,
    )
    metadata = {
        "utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "scipy": scipy.__version__,
        "numpy": np.__version__,
        "julia_bigfloat_bits": 256,
        "fixtures_sha256": {str(p.relative_to(root)): sha256(p) for p in (angular_path, radial_path)},
    }
    angular = {}
    radial = {}
    for line in completed.stdout.splitlines():
        fields = line.split("|")
        if fields[0] == "META":
            metadata[fields[1]] = fields[2]
        elif fields[0] == "ANG":
            convert = (lambda s: Decimal.from_float(float(s))) if fields[1] == "double" else Decimal
            angular[(fields[1], *map(int, fields[2:5]))] = (convert(fields[5]), convert(fields[6]))
        elif fields[0] == "RAD":
            radial[fields[1]] = Decimal.from_float(float(fields[2])) if fields[1] == "double" else Decimal(fields[2])
        else:
            raise ValueError(f"Unexpected Julia output: {line}")
    for p in ("double", "quad"):
        metadata[f"backend_{p}_sha256"] = sha256(Path(metadata[f"backend_{p}"]))
    rows = []
    with localcontext() as context:
        context.prec = 70
        for fixture, case in zip(fixtures, cases):
            reference = Decimal(fixture[-1]) / Decimal(fixture[-2])
            sv, sd = special.obl_ang1(*case, 0.3)
            values = {"scipy": (Decimal.from_float(float(sv)), Decimal.from_float(float(sd)))}
            values.update({p: angular[(p, *case)] for p in ("double", "quad")})
            row = {"m": case[0], "n": case[1], "c": case[2], "eta": "0.3", "reference_logderivative": str(reference)}
            for implementation, (value, derivative) in values.items():
                ratio = derivative / value
                row[implementation] = {
                    "value": str(value), "derivative": str(derivative),
                    "logderivative": str(ratio),
                    "relative_error": str(abs((ratio - reference) / reference)),
                }
            rows.append(row)
        rv, _ = special.pro_rad1(0, 1, 1.0, 2.0)
        radial["scipy"] = Decimal.from_float(float(rv))
        ref = Decimal(radial_ref)
        control = {
            "m": 0, "n": 1, "c": 1, "x": 2, "reference_value": radial_ref,
            "results": {name: {"value": str(value), "relative_error": str(abs((value-ref)/ref))} for name,value in radial.items()},
        }
    result = {"environment": metadata, "oblate_angular_logderivative": rows, "prolate_radial_first_kind_control": control}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
