#ifndef __GENESIS_FIELD__
#define __GENESIS_FIELD__

#include <vector>
#include <iostream>
#include <string>
#include <complex>

class Beam;


#include "Undulator.h"
#include "FieldSolver.h"

using namespace std;


extern const double vacimp;
extern const double eev;


class Field{
 public:
   Field();
   virtual ~Field();
   void initDiagnostics(int);
   void diagnostics(bool);
   void init(int, int, double, double, double, double,int);
   bool harmonicConversion(int, bool);
   bool subharmonicConversion(int, bool);
   void track(double, Beam *, Undulator *);
   bool getLLGridpoint(double, double, double *, double *,int *);
   void setStepsize(double);
   void disable(double);
   bool isEnabled();
   int getHarm();
   double getRHarm();

   vector< vector< complex<double> > > field;
   double xlambda,dgrid,xks,gridmax,dz_save;
   int ngrid, first;  // first points to first slice due to slippage
   int harm;
   bool polarization;
   double slicelength,s0;


   vector<double> power,xsig,xavg,ysig,yavg ;  // buffer to accumulate before writing it out
   vector<double> nf_intensity,nf_phi,ff_intensity,ff_phi;


 private:
   int idx;
   bool disabled;
   double rharm;
   double accuslip;
     
   FieldSolver solver;
};


inline void Field::disable(double conv)
{
  if (disabled==false){  // check whether it hasn't been disabled before
    rharm=harm;         // assign current double harmonic with the given harmonic
  }
  rharm=rharm*conv;     // convert to new harmonic, might be even a non-integer.
  disabled=true;         // disbable it.

}

inline bool Field::isEnabled()
{
  return !disabled;
}


inline double Field::getRHarm()
{
  return rharm;
}




inline int Field::getHarm()
{
  return harm;
}

#endif
