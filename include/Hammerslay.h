/*
 *  Hammerslay.h
 *  Genesis
 *
 *  Created by Sven Reiche on 5/26/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "Sequence.h"

#ifndef __GENESIS_HAMMERSLAY__
#define __GENESIS_HAMMERSLAY__



class Hammerslay : public Sequence {
 public:
	Hammerslay(unsigned int = 0 );
	~Hammerslay();
	void set(unsigned int);
	double getElement();
private:
	double base;
	unsigned int idx;
};


#endif
