#ifndef __GENESIS_FIELD__
#define __GENESIS_FIELD__

#include <vector>
#include <iostream>
#include <string>
#include <complex>


#ifdef FFTW
#include <fftw3.h>
#endif



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
   bool get_global_stat();
   void set_global_stat(bool);
   void setOutput(bool,bool,bool,bool);
   bool outputFFT();
   bool outputSpatial();
   bool outputIntensity();
   bool dumpFieldEnabled();
   int getHarm();
   double getRHarm();


   vector< vector< complex<double> > > field;
   double xlambda,dgrid,xks,gridmax,dz_save;
   int ngrid, first;  // first points to first slice due to slippage
   int harm;
   bool polarization;
   double slicelength,s0;
   double accuslip;
   FieldSolver solver;

   vector<double> power,xsig,xavg,ysig,yavg ;  // buffer to accumulate before writing it out
   vector<double> txsig,txavg,tysig,tyavg ;  // buffer to accumulate before writing it out
   vector<double> nf_intensity,nf_phi,ff_intensity,ff_phi;
   // global variables   - energy is proportional to the mean power
   vector<double> energy,gl_xsig,gl_xavg,gl_ysig,gl_yavg,gl_nf_intensity,gl_ff_intensity;
#ifdef FFTW
   vector<double> gl_txsig, gl_txavg, gl_tysig, gl_tyavg;  
#endif

 private:
   int idx;
   bool disabled;
   double rharm;
   bool out_global, doFFT,doSpatial, doIntensity;
   bool doDumpField; // controls write of field grid to .dfl files (can be OFF, if intensity projects are sufficient)

   complex<double> *in;
   complex<double> *out;
#ifdef FFTW
   fftw_plan p;
#endif
     

};


inline bool Field::outputFFT(){ return doFFT;}
inline bool Field::outputSpatial(){ return doSpatial;}
inline bool Field::outputIntensity(){ return doIntensity;}
inline bool Field::dumpFieldEnabled(){ return doDumpField;}
inline bool Field::get_global_stat(){return out_global;}
inline void Field::set_global_stat(bool in) {out_global=in;}
inline void Field::setOutput(bool nofft_in, bool noSpatial_in, bool noInten_in, bool noDumpField_in) {
  doFFT = !nofft_in;
  doSpatial = !noSpatial_in;
  doIntensity = !noInten_in;
  doDumpField = !noDumpField_in;
}

inline void Field::disable(double conv)
{
  if (disabled==false){  // check whether it hasn't been disabled before
    rharm=harm;         // assign current double harmonic with the given harmonic
  }
  rharm=rharm*conv;     // convert to new harmonic, might be even a non-integer.
  disabled=true;         // disable it.

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
