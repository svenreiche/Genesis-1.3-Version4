#!/bin/bash

output=$(./build/genesis4)
if [[ "$output" != *"Usage: genesis4"* ]]; then
  echo "Error: Expected 'Usage: genesis4' in output but got: $output"
  exit 1
else
  echo "Genesis4 seems to have built correctly; 'usage' found in output."
fi

### BASIC TEST ###

# uncomment if you wish to disable the short test simulation run
# exit 0

# get absolute path to our new genesis4 binary
G4=`realpath ./build/genesis4`

###
# Run test case ('Example1-SteadyState': steady-state simulation, completes in approx 10s on CL's desktop PC)
#
# Currently this tests only for successful completion (the output itself
# is analyzed in the next step).
###
SEMAFILE="testrun.sema"
RUNDIR="./examples/Example1-SteadyState/"
INFILE="Example1.in" # in RUNDIR

echo "Entering directory ${RUNDIR}"
cd $RUNDIR
# delete old semaphore file from previous runs (if exists)
rm -f $SEMAFILE

# now run the simulation
$G4 --semaphore-file $SEMAFILE $INFILE
RESULT=$?
echo "exit code of G4 test run: ${RESULT}"
if [[ "$RESULT" -ne "0" ]] ; then
	echo "non-zero exit code indicates issue, stopping"
	exit 1
fi

# test for existence of semaphore file signalling successful completion of G4 run
if [[ ! -f $SEMAFILE ]] ; then
	echo "semaphore file '${SEMAFILE}' not found"
	exit 1
fi

# list directory contents
echo
echo "directory contents after G4 test run ..."
ls -l


# stop here
# exit 0

###
# analysis of simulation results
# DEMO: report and check power at end of undulator
# (reference value for check is coded in TEST_power.py script)
###

echo
echo
./TEST_power.py
RESULT=$?
if [[ "$RESULT" -ne "0" ]] ; then
	echo "test script with non-zero exit code, stopping"
	exit 1
fi
echo
echo

###
