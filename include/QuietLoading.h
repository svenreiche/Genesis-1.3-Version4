/*
 *  QuietLoading.h
 *  Genesis
 *
 *  Created by Sven Reiche on 5/14/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include <iostream>
#include <complex>
#include <math.h>


#include "Particle.h"
#include "Sequence.h"
#include "RandomU.h"
#include "Hammerslay.h"
#include "Inverfc.h"


#ifndef __GENESIS_QUIETLOADING__
#define __GENESIS_QUIETLOADING__

extern const double vacimp;
extern const double eev;

typedef struct
{
  double gamma;
  double delgam;
  double current;
  double ex;
  double ey;
  double xcen;
  double ycen;
  double pxcen;
  double pycen;
  double betax;
  double betay;
  double alphax;
  double alphay;
  double bunch;
  double bunchphase;
  double emod;
  double emodphase;
} BeamSlice;


using namespace std;

class QuietLoading  {
public:
	QuietLoading();
	~QuietLoading();
        void loadQuiet(Particle *beam, BeamSlice *, int, int, double,int); 
        void init(bool, int *); 
	
	
private:
        Sequence *sx,*sy,*st,*spx,*spy,*sg;
};


#endif
