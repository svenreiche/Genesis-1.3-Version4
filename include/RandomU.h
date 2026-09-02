/*
 *  RandomU.h
 *  Genesis
 *
 *  Created by Sven Reiche on 5/26/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>

#include "Sequence.h"

#ifndef __GENESIS_RANDOMU__
#define __GENESIS_RANDOMU__



class RandomU : public Sequence{
public:
	RandomU(unsigned int = 0 );
	~RandomU();
	void set(unsigned int);
	double getElement();
private:

	int iv[32],iy;
	int iseed,iseed2;
};


// The independent random streams of a run. Every consumer of seedFromIndex
// below names its own, so that two consumers keyed on the same slice do not
// draw the same sequence.
enum class SeedStream : unsigned int {
	shotnoise  = 1,
	incoherent = 2,
	one4one    = 3,
	sddsphase  = 4
};


// Derives a generator seed from the seed the user gave, a global index and the
// stream the caller belongs to.
//
// The point of it is that the index is GLOBAL, typically the slice number in
// the whole time window rather than the position in one core's share of it.
// A stream keyed this way delivers the same numbers to a given slice however
// many cores the run is spread over and whichever core ends up owning that
// slice, which is what makes a seeded run reproducible across core counts.
//
// The mixing is a splitmix64 finaliser, chosen because it is short, has no
// state and is fully specified in integer arithmetic, so it gives the same
// answer on every platform and compiler. The result is reduced to the range
// the ran2 generator in RandomU needs, which is 1 up to but excluding 2147483563.
unsigned int seedFromIndex(unsigned int base, unsigned long index, SeedStream stream);


#endif


