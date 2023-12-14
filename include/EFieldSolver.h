#ifndef __GENESIS_EFIELDSOLVER__
#define __GENESIS_EFIELDSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <mpi.h>

//#include "readfieldHDF5.h"
//#include "Undulator.h"
//#include "FieldSolver.h"


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
    void shortRange(vector<Particle> *, double, double);
    void longRange(Beam *beam, double gamma, double aw);
    double getEField(double x, double y, double theta);
    [[nodiscard]] bool hasShortRange() const;
    double diagField();
private:
    auto analyseBeam(vector<Particle> *beam);
    void constructLaplaceOperator();

    vector<double> work1, work2, fcurrent, fsize;
    vector<complex<double> > csrc, clow, cmid, cupp, celm, gam;
    vector<complex<double> > csource,cfield;
    vector<double> lupp, lmid, llow, rlog, vol;
    int nz, nphi, ngrid;
    double rmax, ks, xcen, ycen, dr, rbound;
    bool longrange;


};

inline bool EFieldSolver::hasShortRange() const{
    return (nz>0) & (ngrid > 2);
}

inline double EFieldSolver::diagField() {
    return std::abs(cfield[nphi*ngrid]);
}

#endif
