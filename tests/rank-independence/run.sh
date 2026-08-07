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
: "${PYTHON:=python3}"
: "${RANKS:=1 2 3}"
: "${WORKDIR:=$PWD/genesis4-test-work}"

# An MPIEXEC given from outside is taken as the one to use. Otherwise a launcher
# is chosen below, because the generic name on the path does not always belong
# to the MPI the executable was linked against.
mpiexec_given="${MPIEXEC:-}"

if [ ! -x "$GENESIS4" ] ; then
	echo "no genesis4 executable at $GENESIS4"
	exit 1
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

# Reports whether the log of a run shows the size of communicator that was
# asked for. Genesis writes "node" for one rank and "nodes" for more.
comm_size_is() {
	grep -Eq "^MPI-Comm Size: $1 nodes?\$" "$2"
}

# The launchers to try, most specific first. The name recorded by CMake belongs
# to the MPI that was linked in, so it is the best guess available. The generic
# names come next, then the implementation-specific ones: a package manager may
# register mpirun for the implementation it installed while leaving mpiexec
# pointing at a different one, which is exactly the case this has to survive.
launcher_candidates() {
	if [ -n "$mpiexec_given" ] ; then
		echo "$mpiexec_given"
		return
	fi

	local cache
	for cache in "$(dirname "$GENESIS4")/CMakeCache.txt" \
	             "$(dirname "$GENESIS4")/../CMakeCache.txt" ; do
		if [ -f "$cache" ] ; then
			sed -n 's/^MPIEXEC_EXECUTABLE:[A-Z]*=//p' "$cache"
		fi
	done

	printf '%s\n' mpiexec mpirun mpiexec.mpich mpirun.mpich mpiexec.hydra \
	               mpiexec.openmpi mpirun.openmpi
}

# The rank count used to test a launcher: the smallest above one that the
# comparison will actually use. With no such count no launcher is needed.
probe_ranks=""
for n in $RANKS ; do
	if [ "$n" -gt 1 ] ; then
		probe_ranks="$n"
		break
	fi
done

launcher=""
reason="it was not found"
while IFS= read -r candidate ; do
	[ -n "$candidate" ] || continue
	command -v "$candidate" > /dev/null 2>&1 || continue

	if [ -z "$probe_ranks" ] ; then
		launcher="$candidate"
		break
	fi

	# The standard input has to be closed off: it carries the list of
	# candidates, and a launcher which drains it leaves nothing to read for
	# the next turn of the loop.
	if "$candidate" -n "$probe_ranks" "$GENESIS4" -o probe sase.in \
			< /dev/null > log_probe.txt 2>&1 ; then
		if comm_size_is "$probe_ranks" log_probe.txt ; then
			launcher="$candidate"
			break
		fi
		reason="it reported $(grep -cE '^MPI-Comm Size:' log_probe.txt) communicator(s) of another size"
	else
		reason="it exited with an error: $(tail -3 log_probe.txt | tr '\n' ' ')"
	fi

	if [ -z "$mpiexec_given" ] ; then
		echo "$candidate cannot start a job of $probe_ranks ranks; $reason"
	fi
done <<EOF
$(launcher_candidates)
EOF

if [ -z "$launcher" ] ; then
	if [ -n "$mpiexec_given" ] ; then
		# The launcher was named explicitly, so a failure to use it is an
		# error rather than something to work around.
		echo "the MPI launcher $mpiexec_given cannot start a job of $probe_ranks ranks; $reason"
		echo "it probably does not match the library $GENESIS4 was linked against"
		exit 1
	fi
	echo "no MPI launcher on this system starts a job of more than one rank, skipping"
	echo "the launcher has to match the MPI library $GENESIS4 was linked against"
	exit 77
fi

MPIEXEC="$launcher"
echo "using MPI launcher $MPIEXEC"

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
	if ! comm_size_is "$nranks" "log_$rootname.txt" ; then
		echo "genesis4 did not run on $nranks rank(s) for $deck; it reported:"
		grep -E "^MPI-Comm Size:" "log_$rootname.txt" | sort | uniq -c
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
