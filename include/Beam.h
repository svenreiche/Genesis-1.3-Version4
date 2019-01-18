#ifndef __GENESIS_BEAM__
#define __GENESIS_BEAM__

#include <vector>
#include <iostream>
#include <string>

#include "Particle.h"
#include "Undulator.h"
#include "BeamSolver.h"
#include "Incoherent.h"
#include "Sorting.h"
#include "Collective.h"

using namespace std;

extern const double ce;

class Beam{
 public:
   Beam();
   virtual ~Beam();
   void initDiagnostics(int);
   void diagnostics(bool,double);
   void diagnosticsStart();
   void init(int, int, double,double, double,bool);
   void initSorting(int,int,bool,bool);
   void initEField(double,int,int,int,double);
   void initIncoherent(int, int, bool,bool);
   void initWake(unsigned int, double, double *, double, bool);
   bool harmonicConversion(int,bool);
   bool subharmonicConversion(int,bool);
   int sort();
   void track(double, vector<Field *> *, Undulator *);


   void setBunchingHarmonicOutput(int harm_in);
   int getBunchingHarmonics();

   vector< vector<Particle> > beam;
   vector<double> current,eloss;
   double reflength,slicelength;   // for conversion of theta in Particle to real position
   double s0;         // starting position of the time-window
   bool one4one;     // flag whether one4one simulation is done
   int nbins;

   // output buffer
   vector<double> zpos,gavg,gsig,xavg,xsig,yavg,ysig,pxavg,pyavg,bunch,bphi,efld;
   vector<double> bx,by,ax,ay,ex,ey,cu;
   vector< vector<double> > bh,ph;  // harmonic bunching and bunching phase
   
 private:
   BeamSolver solver;
   Incoherent incoherent;
   Collective col;
   Sorting sorting;
   int idx;
   int bharm;
};

inline void Beam::initIncoherent(int base, int rank, bool spread, bool loss){
  incoherent.init(base,rank,spread,loss);
  return;
}

inline void Beam::initEField(double rmax, int ngrid, int nz, int nphi, double lambda){
  solver.initEField(rmax,ngrid,nz,nphi,lambda);
  return;
}

inline void Beam::initWake(unsigned int ns, double ds, double *wakeres, double ztrans, bool transient){
  col.initWake(ns, ds, wakeres, ztrans, transient);
}


inline void Beam::setBunchingHarmonicOutput(int harm_in)
{
  bharm=harm_in;
}

inline int Beam::getBunchingHarmonics()
{
  return bharm;
}

#endif
