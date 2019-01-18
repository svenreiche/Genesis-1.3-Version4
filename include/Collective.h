#ifndef __GENESIS_COLLECTIVE__
#define __GENESIS_COLLECTIVE__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <math.h>

#include "mpi.h"
#include "Particle.h"
#include "Undulator.h"

class Beam;

extern bool MPISingle; 
extern const double ce;   

using namespace std;


class Collective{
 public:
   Collective();
   virtual ~Collective();
   void initWake(unsigned int, double, double *, double, bool);
   void apply(Beam *,Undulator *, double );

 private:
   bool transient,hasWake;
   double ztrans;
   double ds;
   unsigned int ns;
   double *wakeres, *current;
};

#endif
