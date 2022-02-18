#ifndef __GENESIS_FIELDSOLVER__
#define __GENESIS_FIELDSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>

#ifdef FFTW
#include <fftw3.h>
#endif


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
   void init(int);
   void initSourceFilter(bool, double, double, double);

 private:
   int ngrid;
   double delz_save;
   complex<double> cstep;
   vector< complex< double > > r,c,cbet,cwet,crsource;

   bool   difffilter_;
   double filtcutx_, filtcuty_, filtsig_;
   vector<double> sigmoidx_,sigmoidy_;

#ifdef FFTW
    complex<double> *in, *out;
    fftw_plan p,pi;

#endif
   void filterSourceTerm();
   void ADI(vector<complex< double > > &);
   void tridagx(vector<complex< double > > &);
   void tridagy(vector<complex< double > > &);
   
};

#endif
