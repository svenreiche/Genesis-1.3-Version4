#!/usr/bin/env python3
#
# Compares the output files that run.sh produced for the same input deck at
# different MPI rank counts, and reports the result through the exit code:
#
#   0 => OK
#   1 => not OK
#
# The comparison is exact. Every quantity examined here is recorded per slice,
# and a slice is tracked by exactly one rank from the start of the run to the
# end, so distributing the same window over more ranks reorders no arithmetic
# and the results are expected to agree bit for bit. Quantities under /Global
# are deliberately not examined, because they are formed with MPI collectives
# whose summation order does legitimately depend on the rank count.

import sys

import h5py
import numpy as np

# Per-slice quantities to compare, and, for each deck, the one of them which is
# checked to respond to a change of the seed. Without that second check the
# suite would still pass if the randomness were removed altogether.
DATASETS = [
    "Beam/bunching",
    "Beam/energy",
    "Beam/energyspread",
    "Beam/xsize",
    "Field/power",
]

DECKS = {
    "sase": "Beam/bunching",
    "sponrad": "Beam/energyspread",
    "one4one": "Beam/bunching",
}

# The witness quantity has to move by at least this much, relative to its own
# magnitude, when the seed is changed.
MIN_SEED_RESPONSE = 1e-3


def read(rootname):
    path = f"{rootname}.out.h5"
    with h5py.File(path, "r") as f:
        missing = [name for name in DATASETS if name not in f]
        if missing:
            raise KeyError(f"{path} does not contain {', '.join(missing)}")
        return {name: np.array(f[name]) for name in DATASETS}


def relative_difference(a, b):
    scale = np.max(np.abs(a)) if a.size else 0.0
    if scale == 0.0:
        scale = 1.0
    return float(np.max(np.abs(a - b)) / scale)


def main(argv):
    ranks = [int(n) for n in argv[1:]]
    if len(ranks) < 2:
        print("usage: compare.py NRANKS NRANKS [NRANKS ...]")
        return 1

    reference_ranks = ranks[0]
    failures = []

    for deck, witness in DECKS.items():
        reference = read(f"{deck}_n{reference_ranks}")

        for n in ranks[1:]:
            other = read(f"{deck}_n{n}")
            disagreements = []
            for name in DATASETS:
                a, b = reference[name], other[name]
                if a.shape != b.shape:
                    disagreements.append(
                        f"{deck}: {name} has shape {a.shape} on "
                        f"{reference_ranks} rank(s) but {b.shape} on {n}"
                    )
                elif not np.array_equal(a, b):
                    disagreements.append(
                        f"{deck}: {name} differs between {reference_ranks} and "
                        f"{n} rank(s) by {relative_difference(a, b):.3e} "
                        f"relative to its peak value"
                    )
            if disagreements:
                failures.extend(disagreements)
            else:
                print(f"OK: {deck} is identical on {reference_ranks} and {n} rank(s)")

        # Control: the same comparison applied to a run that differs only in its
        # seed has to report a difference, otherwise the agreement above would
        # not mean anything.
        altered = read(f"{deck}_altseed")
        response = relative_difference(reference[witness], altered[witness])
        if response < MIN_SEED_RESPONSE:
            failures.append(
                f"{deck}: changing the seed moved {witness} by only "
                f"{response:.3e}, so this deck is not sensitive to the random "
                f"sequences and the comparison above proves nothing"
            )
        else:
            print(f"OK: {deck} responds to a change of seed ({witness} moves by {response:.3e})")

    if failures:
        print()
        for message in failures:
            print(f"FAILED: {message}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
