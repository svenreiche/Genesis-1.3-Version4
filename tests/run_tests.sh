#!/bin/bash
#
# Runs the Genesis 1.3 test suite. Every immediate subdirectory of this one that
# contains an executable script named run.sh is treated as a test case and is
# run with its own directory as the working directory.
#
# The suite can be run directly:
#
#     ./tests/run_tests.sh [path to genesis4]
#
# or through CTest from a build tree:
#
#     ctest --test-dir build --output-on-failure
#
# Environment variables:
#   GENESIS4  path to the executable under test  (default: build/genesis4)
#   MPIEXEC   MPI launcher                       (default: mpiexec)
#   PYTHON    Python interpreter with h5py       (default: python3)
#   RANKS     rank counts to compare             (default: set by each case)
#   WORKDIR   directory for the generated files  (default: ./genesis4-test-work)
#
# Exit code:
#   0 => every case passed or was skipped
#   1 => at least one case failed

set -u

here="$(cd "$(dirname "$0")" && pwd)"

: "${GENESIS4:=$here/../build/genesis4}"
: "${MPIEXEC:=mpiexec}"
: "${PYTHON:=python3}"
: "${WORKDIR:=$PWD/genesis4-test-work}"

if [ $# -ge 1 ] ; then
	GENESIS4="$1"
fi

# The cases are run from their own directories, so a relative path to the
# executable would no longer resolve.
case "$GENESIS4" in
	/*) ;;
	*) GENESIS4="$PWD/$GENESIS4" ;;
esac
case "$WORKDIR" in
	/*) ;;
	*) WORKDIR="$PWD/$WORKDIR" ;;
esac

export GENESIS4 MPIEXEC PYTHON WORKDIR

if [ ! -x "$GENESIS4" ] ; then
	echo "no genesis4 executable at $GENESIS4"
	echo "pass its path as an argument or set GENESIS4"
	exit 1
fi

echo "testing $GENESIS4"
echo "writing generated files to $WORKDIR"

passed=0
skipped=0
failed=0
failed_names=""

for runner in "$here"/*/run.sh ; do
	[ -f "$runner" ] || continue
	name="$(basename "$(dirname "$runner")")"

	echo
	echo "--- $name ---"
	( cd "$(dirname "$runner")" && bash ./run.sh )
	result=$?

	if [ "$result" -eq 0 ] ; then
		echo "PASS $name"
		passed=$((passed + 1))
	elif [ "$result" -eq 77 ] ; then
		echo "SKIP $name"
		skipped=$((skipped + 1))
	else
		echo "FAIL $name"
		failed=$((failed + 1))
		failed_names="$failed_names $name"
	fi
done

echo
echo "$passed passed, $skipped skipped, $failed failed"
if [ "$failed" -ne 0 ] ; then
	echo "failing cases:$failed_names"
	exit 1
fi
exit 0
