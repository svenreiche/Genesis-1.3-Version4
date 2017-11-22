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


#endif


