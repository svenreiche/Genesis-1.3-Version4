#ifndef __GENESIS_EFIELDSOLVER__
#define __GENESIS_EFIELDSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <mpi.h>

#include "Particle.h"

class Beam;

using namespace std;

extern const double vacimp;
extern const double eev;



class EFieldSolver {
public:
    EFieldSolver();
    virtual ~EFieldSolver();
    void init(double, int, int, int, double, bool);
    void shortRange(vector<Particle> *, double, double, int);
    void longRange(Beam *beam, double gamma, double aw);
    double getEField(unsigned long i);
    bool hasShortRange() const;
    void allocateForOutput(unsigned long nslice);
    double getSCField(int);

private:
    void analyseBeam(vector<Particle> *beam);
    void constructLaplaceOperator();
    void tridiag();

    vector<double> work1, work2, fcurrent, fsize;  // used for long range calculation
    vector<complex<double> > cwork;
    vector<int> idxr;
    vector<double> lmid, rlog, vol, ldig;
    vector<complex<double> > csrc, clow, cmid, cupp, celm, gam; // used for tridiag routine
    vector<double> ez,efield;

    int nz, nphi, ngrid, rank;
    double rmax, ks, xcen, ycen, dr;
    bool longrange;

};

inline double EFieldSolver::getSCField(int islice) {
    return efield[islice]*511000;  // convert from Lorent mass unit to eV /m
}

inline bool EFieldSolver::hasShortRange() const{
    return (nz>0) & (ngrid > 2);
}

#endif
