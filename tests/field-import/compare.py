#!/usr/bin/env python3
"""
Compares the field dumps a run wrote after an import with the dumps it
imported. Exits with 0 if every reimported dump reproduces its original and
with 1 otherwise.
"""
import sys

import h5py
import numpy as np

# The import divides by the factor the dump multiplied by, so the fields come
# back to within rounding. The scalars are copied or recomputed from the same
# inputs.
RTOL = 1e-12

SCALARS = ("wavelength", "slicespacing", "refposition", "gridsize",
           "gridpoints", "slicecount")


def compare(original, copy):
    """Returns the list of differences between two field dumps."""
    problems = []
    with h5py.File(original, "r") as fa, h5py.File(copy, "r") as fb:
        for name in SCALARS:
            a = fa[name][()].item()
            b = fb[name][()].item()
            if not np.isclose(a, b, rtol=RTOL, atol=0.0):
                problems.append(f"{name}: {a!r} in the dump, {b!r} after import")
        if problems:
            return problems

        nslice = int(fa["slicecount"][()].item())
        worst = 0.0
        peak = 0.0
        for i in range(1, nslice + 1):
            for part in ("field-real", "field-imag"):
                a = fa[f"slice{i:06d}/{part}"][()]
                b = fb[f"slice{i:06d}/{part}"][()]
                if not (np.isfinite(a).all() and np.isfinite(b).all()):
                    return [f"{part} of slice {i} is not finite"]
                worst = max(worst, float(np.max(np.abs(a - b))))
                peak = max(peak, float(np.max(np.abs(a))))
        if peak == 0.0:
            problems.append("the dump holds no field, nothing was compared")
        elif worst > RTOL * peak:
            problems.append(f"field differs by {worst:.3e}, peak is {peak:.3e}")
    return problems


def main():
    run = sys.argv[1]
    failed = False
    for suffix in (".fld.h5", ".h3.fld.h5"):
        problems = compare("dump" + suffix, run + suffix)
        if problems:
            failed = True
            print(f"FAIL {run}{suffix}")
            for line in problems:
                print(f"     {line}")
        else:
            print(f"ok   {run}{suffix} reproduces dump{suffix}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
