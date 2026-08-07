# Test suite

Each subdirectory of this one is a test case. A case consists of the input decks
it needs and a script named `run.sh`, which runs them and decides whether the
result is acceptable. The driver `run_tests.sh` runs every case and reports a
summary.

## Running the suite

From a source tree with a build in `build/`:

```
./tests/run_tests.sh
```

Against some other executable:

```
./tests/run_tests.sh /path/to/genesis4
```

Through CTest, which is the form a packaging recipe is most likely to use:

```
ctest --test-dir build --output-on-failure
```

The generated output files are written to `./genesis4-test-work`, or to
`$WORKDIR` if that is set, and are left in place when a case fails so that they
can be examined.

## Settings

| Variable   | Meaning                             | Default          |
| ---------- | ----------------------------------- | ---------------- |
| `GENESIS4` | executable under test               | `build/genesis4` |
| `MPIEXEC`  | MPI launcher                        | `mpiexec`        |
| `PYTHON`   | Python interpreter, needs `h5py`    | `python3`        |
| `RANKS`    | MPI rank counts to compare          | `1 2 3`          |
| `WORKDIR`  | where the generated files are put   | `./genesis4-test-work` |

## Adding a case

Create a directory containing the input decks and a `run.sh` which exits with

| Code | Meaning                                                        |
| ---- | -------------------------------------------------------------- |
| `0`  | the case passed                                                 |
| `1`  | the case failed                                                 |
| `77` | the case could not be run, for instance because a dependency of it is missing |

The driver picks the new directory up without further registration, as does
`tests/CMakeLists.txt`, which turns each one into a CTest test. A case should
report enough on failure to say how badly the expectation was missed, and should
skip rather than fail when something it depends on is absent, so that a build
without MPI or without `h5py` does not produce a spurious failure.

Keep the runs short. The decks here are chosen to finish in a few seconds and
are sized for the property under test alone; none of them is a physically
meaningful simulation, and none of their numbers should be quoted as a
benchmark.

## Cases

### `rank-independence`

Requires that a run gives the same answer however many MPI ranks it is
distributed over. Three decks cover the three independent sources of randomness:
the shot noise of the beam, the quantum fluctuation of the incoherent radiation,
and the particle loading of a one-to-one simulation. Each is run at every rank
count in `RANKS` and the per-slice output is required to agree exactly.

Each deck is also run once more with a different seed, and the comparison is
required to detect that. A test which cannot distinguish two different random
realizations would pass for the wrong reason.

Every run is additionally required to report the number of MPI ranks it was
asked for. A launcher which does not match the library the executable was linked
against starts several independent single-rank jobs instead, all of which write
the same output, so without that check the case would agree with itself and pass
having tested nothing.

The window holds 24 slices, which divides evenly by every rank count in the
default set. Genesis extends the time window so that every rank holds the same
number of slices, so a window whose length were not a multiple of the rank count
would change between the runs for a second and unrelated reason.

Not yet covered: the beam imported from an SDDS distribution, whose phase
reconstruction and shot noise are seeded separately in `SDDSBeam`. That path
needs a distribution file which the repository does not currently carry.
