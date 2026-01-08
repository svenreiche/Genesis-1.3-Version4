#ifndef __GENESIS_FIELDSOLVERFFT__
#define __GENESIS_FIELDSOLVERFFT__

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
#include "FieldSolver.h"

using namespace std;


class FieldSolverFFT : public FieldSolver{
 public:
   ~FieldSolverFFT();
   void init(double,double,double,int);
   void advance(double, Field *, Beam *, Undulator *);
   void initSourceFilter(double,double,double,bool);

 private:
    int ngrid {0} ;
    double delz_save {0};
    double ks {1};
    double dk {1};
    double xc,yc,sig;
    bool hasPlan {false};
    bool doFilter_ {false};
    complex<double> *in, *out;
    vector<complex<double> > uf, sf, k1, k2, k3, k4, K2,sigmoid_;

#ifdef FFTW
    fftw_plan p,ip;
#endif
    vector< complex< double > > crsource;
    void FFT(vector<complex< double > > &);
};

#endif
