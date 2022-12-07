#ifndef __GENESIS_BEAMSOLVER__
#define __GENESIS_BEAMSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>



class Field;


#include "Undulator.h"
#include "EFieldSolver.h"
#include "TrackBeam.h"



using namespace std;


class BeamSolver{
 public:
   BeamSolver();
   virtual ~BeamSolver();

   void initEField(double rmax, int ngrid, int nz, int nphi, double lambda, bool longr, bool redLF);

   void advance(double, Beam *, vector< Field *> *, Undulator *);
   void track(double, Beam *, Undulator *,bool);
   void applyR56(Beam *, Undulator *, double);

 private:
 
   complex <double> cpart;

   vector< double > rharm;
   vector< complex <double > > rpart;
   vector<double> esc;
   
   double ez;
   double xks,xku;

   double theta,gamma,btpar;
   double k2gg,k2pp,k3gg,k3pp;

   bool onlyFundamental;
 

   void RungeKutta(double);
   void ODE(double,double);

   EFieldSolver efield;
   TrackBeam tracker;

};


inline void BeamSolver::initEField(double rmax, int ngrid, int nz, int nphi, double lambda, bool longr, bool redLF){
  efield.init(rmax,ngrid,nz,nphi,lambda,longr,redLF);
  return;
}


inline void BeamSolver::track(double dz, Beam *beam, Undulator *und, bool last)
{
  tracker.track(dz,beam,und,last);
  return;
}

inline void BeamSolver::applyR56(Beam *beam, Undulator *und, double reflen){
  tracker.applyR56(beam,und,reflen);
  return;
}

#endif
