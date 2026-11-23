#!/bin/bash
#
# Writes a fundamental and a third harmonic to field dumps and reads them back
# with &importfield into a run which takes its time window from the harmonic
# file. The reimported fields must be the ones which were written.
#
# Exit code:
#   0  => OK
#   1  => not OK
#   77 => skipped, because a Python interpreter with h5py is unavailable

set -u

srcdir="$(cd "$(dirname "$0")" && pwd)"

: "${GENESIS4:=$srcdir/../../build/genesis4}"
: "${PYTHON:=python3}"
: "${WORKDIR:=$PWD/genesis4-test-work}"

if [ ! -x "$GENESIS4" ] ; then
	echo "no genesis4 executable at $GENESIS4"
	exit 1
fi

work="$WORKDIR/field-import"
rm -rf "$work"
mkdir -p "$work" || exit 1
cp "$srcdir"/*.in "$srcdir"/*.lat "$work"/ || exit 1
cd "$work" || exit 1

# The case runs on a single rank and so needs no launcher.
for deck in generate import_notime ; do
	echo "running $deck"
	if ! "$GENESIS4" "$deck.in" > "log_$deck.txt" 2>&1 ; then
		echo "genesis4 failed on $deck.in:"
		tail -20 "log_$deck.txt"
		echo "output files kept in $work"
		exit 1
	fi
done

# Only the comparison needs h5py, so a failed import is reported without it.
if ! "$PYTHON" -c "import h5py, numpy" > /dev/null 2>&1 ; then
	echo "$PYTHON has no h5py or no numpy, skipping the comparison"
	exit 77
fi

echo
"$PYTHON" "$srcdir/compare.py" import_notime
result=$?
if [ "$result" -ne 0 ] ; then
	echo "output files kept in $work"
fi
exit $result
