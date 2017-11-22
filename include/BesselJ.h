/*
 *  BesselJ.h
 *  Genesis
 *
 *  Created by Sven Reiche on 12/5/11.
 *  Copyright 2011 Paul Scherrer Institut. All rights reserved.
 *
 */



#ifndef __GENESIS_BESSELJ__
#define __GENESIS_BESSELJ__

#include <math.h>

using namespace std;

class BesselJ  {
public:
	BesselJ();
	~BesselJ();
	double value(int, double);
private:
	double BesselJ0(double);
	double BesselJ1(double);
};

#endif
