#ifndef __GENESIS_FIELDSOLVERADI__
#define __GENESIS_FIELDSOLVERADI__

#include <vector>
#include <iostream>
#include <string>
#include <complex>


class Field;
class Beam;

#include "Particle.h"
#include "Undulator.h"
#include "FieldSolver.h"

using namespace std;


class FieldSolverADI : public FieldSolver{
 public:
   ~FieldSolverADI();
   void init(double,double,double,int);
   void advance(double, Field *, Beam *, Undulator *);

 private:
   int ngrid {0};
   double delz_save {0};
   complex<double> cstep;
   vector< complex< double > > r,c,cbet,cwet,crsource;

   void ADI(vector<complex< double > > &);
   void tridagx(vector<complex< double > > &);
   void tridagy(vector<complex< double > > &);
   
};

#endif
