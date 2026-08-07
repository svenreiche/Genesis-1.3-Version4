#!/bin/bash
#
# Runs three input decks at several MPI rank counts and requires that each deck
# gives the same answer every time. The three decks cover the three independent
# sources of randomness in a run: the shot noise of the beam, the quantum
# fluctuation of the incoherent radiation, and the particle loading used for a
# one-to-one simulation.
#
# Exit code:
#   0  => OK
#   1  => not OK
#   77 => skipped, because MPI or a Python interpreter with h5py is unavailable

set -u

srcdir="$(cd "$(dirname "$0")" && pwd)"

: "${GENESIS4:=$srcdir/../../build/genesis4}"
: "${MPIEXEC:=mpiexec}"
: "${PYTHON:=python3}"
: "${RANKS:=1 2 3}"
: "${WORKDIR:=$PWD/genesis4-test-work}"

if [ ! -x "$GENESIS4" ] ; then
	echo "no genesis4 executable at $GENESIS4"
	exit 1
fi
if ! command -v "$MPIEXEC" > /dev/null 2>&1 ; then
	echo "no MPI launcher ($MPIEXEC) on the path, skipping"
	exit 77
fi
if ! "$PYTHON" -c "import h5py, numpy" > /dev/null 2>&1 ; then
	echo "$PYTHON has no h5py or no numpy, skipping"
	exit 77
fi

# A continuous integration runner has fewer cores than the suite has ranks.
# Open MPI declines to start in that situation unless it is told that this is
# intended; the first variable is read by version 4 and earlier, the second by
# version 5 and later. MPICH ignores both.
export OMPI_MCA_rmaps_base_oversubscribe=yes
export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe

work="$WORKDIR/rank-independence"
rm -rf "$work"
mkdir -p "$work" || exit 1
cp "$srcdir"/*.in "$srcdir"/*.lat "$work"/ || exit 1
cd "$work" || exit 1

run() {
	local rootname="$1"
	local nranks="$2"
	local deck="$3"

	if ! "$MPIEXEC" -n "$nranks" "$GENESIS4" -o "$rootname" "$deck" \
			> "log_$rootname.txt" 2>&1 ; then
		echo "genesis4 failed on $deck with $nranks rank(s):"
		tail -20 "log_$rootname.txt"
		return 1
	fi

	# A launcher which does not match the MPI library the executable was
	# linked against starts N independent single-rank jobs rather than one
	# job of N ranks, and says nothing about it. Every one of those jobs
	# would write the same output as the single-rank run, so the comparison
	# below would agree and the case would pass while testing nothing.
	# Genesis reports the size of the communicator it was actually given, so
	# require it to be the size that was asked for.
	if ! grep -Eq "^MPI-Comm Size: $nranks nodes?\$" "log_$rootname.txt" ; then
		echo "genesis4 did not run on $nranks rank(s) for $deck; it reported:"
		grep -E "^MPI-Comm Size:" "log_$rootname.txt" | sort | uniq -c
		echo "the MPI launcher $MPIEXEC probably does not match the library $GENESIS4 was linked against"
		return 1
	fi
	return 0
}

for deck in sase sponrad one4one ; do
	for n in $RANKS ; do
		echo "running $deck on $n rank(s)"
		run "${deck}_n${n}" "$n" "${deck}.in" || exit 1
	done

	# The control run for compare.py: same deck, different seed.
	sed 's/^[[:space:]]*seed[[:space:]]*=.*/seed=192837465/' \
		"${deck}.in" > "${deck}_altseed.in" || exit 1
	echo "running $deck on 1 rank with a different seed"
	run "${deck}_altseed" 1 "${deck}_altseed.in" || exit 1
done

echo
"$PYTHON" "$srcdir/compare.py" $RANKS
result=$?
if [ "$result" -ne 0 ] ; then
	echo "output files kept in $work"
fi
exit $result
