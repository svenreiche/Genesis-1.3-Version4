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
    void init(int);
    void initDiffFilter(bool,double,double,double);
    void getDiag(double,double,double,int);
    void advance(double, Field *, Beam *, Undulator *);

 private:
    bool difffilter_;
    double filtcutx_,filtcuty_,filtsig_;
    int ngrid;     // grid size
    double delz_save;   // integration size of previous step in case step size changes

    // FFTW components
#ifdef FFTW
    complex<double> *in, *out;
    fftw_plan p,pi;
#endif
    // vector for the ADI solver
    complex<double> cstep;
    vector< complex< double > > r,c,cbet,cwet,crsource;

    void filterSourceTerm();
    void ADI(vector<complex< double > > &);
    void tridagx(vector<complex< double > > &);
    void tridagy(vector<complex< double > > &);
};

#endif
