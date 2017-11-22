/*
 *  GaussHermite.h
 *  Genesis
 *
 *  Created by Sven Reiche on 5/14/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include <iostream>
#include <complex>
#include <math.h>


#ifndef __GENESIS_GAUSSHERMITE__
#define __GENESIS_GAUSSHERMITE__

extern const double vacimp;
extern const double eev;

typedef struct
{
  double lambda;
  double power;
  double z0;
  double w0;
  double phase;
  double xcen;
  double ycen;
  double xangle;
  double yangle;
  int nx;
  int ny;
  int harm;
} FieldSlice;


using namespace std;

class GaussHermite  {
public:
	GaussHermite();
	~GaussHermite();
	
	void loadGauss(complex<double> *field, FieldSlice *, double, int);
	
private:
	double Hn(double, int);
	int fac(int);
};


#endif
