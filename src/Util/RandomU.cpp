/*
 *  RandomU.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 5/26/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */


// Routine Ran2 from Numerical Receipe

#include "RandomU.h"

// constructor + destructor

RandomU::RandomU(unsigned int istart)
{

const int ia1=40014;
const int im1=2147483563;
const int iq1=53668;
const int ir1=12211;
const int ntab=32;
int k;
iseed = (1 < istart) ? istart : 1 ;
iseed2=iseed;
for (int i=ntab+7;i>=0;i--){
	k=iseed/iq1;
	iseed=ia1*(iseed-k*iq1)-k*ir1;
	if (iseed < 0 ) iseed+=im1;
		if (i < ntab) iv[i]=iseed;
			}
iy=iv[0]; 
}

RandomU::~RandomU(){}



void RandomU::set(unsigned int istart)
{
  return;
}

double RandomU::getElement()
{
	
	const int ia1=40014;
	const int ia2=40692;
	const int im1=2147483563;
	const int im2=2147483399;
	const int imm1=2147483562;
	const int iq1=53668;
	const int iq2=52774;
	const int ir1=12211;
	const int ir2=3791;
	const int ndiv=67108862;
	const double am=1./2147483563;
	const double rnmx=1.-1.2e-40;	
	
	int k=iseed/iq1;
	iseed=ia1*(iseed-k*iq1)-k*ir1;
	if (iseed <0) iseed=iseed+im1;
	k=iseed2/iq2;
	iseed2=ia2*(iseed2-k*iq2)-k*ir2;
	if (iseed2 < 0) iseed2=iseed2+im2;
	int j=iy/ndiv;
	iy=iv[j]-iseed2;
	iv[j]=iseed; 
	if ( iy < 1 ) iy=iy+imm1;
	return (am*iy < rnmx) ? am*iy : rnmx;
}
