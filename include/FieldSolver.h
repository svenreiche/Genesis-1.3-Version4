#ifndef __GENESIS_FIELDSOLVER__
#define __GENESIS_FIELDSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>


class Field;
class Beam;

#include "Particle.h"
#include "Undulator.h"


using namespace std;


class FieldSolver{
 public:
   FieldSolver();
   virtual ~FieldSolver();
   void getDiag(double,double,double,int);
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
