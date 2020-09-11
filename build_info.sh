#!/bin/bash

# C. Lechner, DESY, 2020-Feb-16
#
# Generate C function that returns essential infos
# about source code (git commit ID etc.) and build
# as string.

F=build_info.c

rm -f $F


echo -n > $F
echo "const char *build_info(void) { " >> $F
echo -n "const char *msg = \"" >> $F
echo -n "compiled by " >> $F
echo -n `whoami` >> $F
echo -n " " >> $F
echo -n `date "+%Y%m%dT%H%M"` >> $F
echo -n " from git commit ID " >> $F
git log --format="%H" -n 1 | tr -d '\n' >> $F

# 20190223: "git status -s -u no" does not show the modified files?
# Lines with untracked files begin with '??', ignore those.
if [[ ! -z $(git status -s | grep --invert "^??") ]]; then
	echo -n " (+ uncommitted changes)" >> $F
fi

echo "\";" >> $F
echo "return(msg);" >> $F
echo "}" >> $F
