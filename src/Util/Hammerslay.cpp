/*
 *  Hammerslay.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 5/26/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "Hammerslay.h"

// constructor/destructor


Hammerslay::Hammerslay(unsigned int base_in)
{
	double bases[26] = {2,3,5,7,11,13,17,19,23,29,
		31,37,41,43,47,53,59,61,67,
		71,73,79,83,89,97,101};
	
	base=bases[base_in];
	idx=0;
}

Hammerslay::~Hammerslay(){}


void Hammerslay::set(unsigned int i)
{
	idx=i;
   if (idx<0) {
	   idx=0;
   }
}

double Hammerslay::getElement()
{
	double xs=0;
	double xsi=1.0;
	unsigned int    i1;
	
	unsigned int i2=++idx;
	
	do{
		xsi=xsi/base;
		i1=i2/static_cast< unsigned int > (base); 
		xs+=(i2-base*i1)*xsi;
		i2=i1;
	} while(i2>0);
	
	return xs;
}					               

