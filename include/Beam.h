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
   void initEField(double,int,int,int,double,bool,bool);
   void initIncoherent(int, int, bool,bool);
   void initWake(unsigned int, unsigned int, double, double *, double *, double *,double *, double,double, bool);
   bool harmonicConversion(int,bool);
   bool subharmonicConversion(int,bool);
   int sort();
   double getSize(int);
   void track(double, vector<Field *> *, Undulator *);
   void setOutput(bool,bool,bool,bool);
   void setWriteFilter(bool,int,int,int);
   bool get_WriteFilter_active();
   int  get_WriteFilter_from();
   int  get_WriteFilter_to();
   int  get_WriteFilter_inc();

   void setBunchingHarmonicOutput(int harm_in);
   int getBunchingHarmonics();
   void set_global_stat(bool);
   bool get_global_stat(void);
   bool outputCurrent();
   bool outputAux();
   bool outputEnergy();
   bool outputSpatial();

   vector< vector<Particle> > beam;
   vector<double> current,eloss,longESC;


   double reflength,slicelength;   // for conversion of theta in Particle to real position
   double s0;         // starting position of the time-window
   bool one4one;     // flag whether one4one simulation is done
   int nbins;

   // output buffer
   vector<double> zpos,gavg,gsig,xavg,xsig,yavg,ysig,pxavg,pyavg,bunch,bphi,efld;
   vector<double> bx,by,ax,ay,ex,ey,cu;
   //   vector<unsigned long long> partcount;
   vector< vector<double> > bh,ph;  // harmonic bunching and bunching phase

   //global values
   vector<double> tgavg, tgsig, txavg,txsig,tyavg, tysig,tbun;  // global values, averaging over the entire beam 
   
 private:
   BeamSolver solver;
   Incoherent incoherent;
   Collective col;
   Sorting sorting;
   int idx;
   int bharm;
   bool do_global_stat;
   bool doCurrent, doSpatial, doEnergy, doAux;

   bool beam_write_filter;
   int beam_write_slices_from, beam_write_slices_to, beam_write_slices_inc;
};

inline bool Beam::outputCurrent(){ return doCurrent;}
inline bool Beam::outputSpatial(){ return doSpatial;}
inline bool Beam::outputEnergy(){ return doEnergy;}
inline bool Beam::outputAux(){ return doAux;}

inline void Beam::initIncoherent(int base, int rank, bool spread, bool loss){
  incoherent.init(base,rank,spread,loss);
  return;
}

inline void Beam::initEField(double rmax, int ngrid, int nz, int nphi, double lambda, bool lngr, bool redLR){
  solver.initEField(rmax,ngrid,nz,nphi,lambda,lngr,redLR);
  return;
}

inline void Beam::initWake(unsigned int ns, unsigned int nsNode, double ds, double *wakeext, double *wakeres, double *wakegeo,double *wakerou, double ztrans, double radius, bool transient){
  col.initWake(ns, nsNode, ds, wakeext, wakeres, wakegeo, wakerou, ztrans, radius, transient);
}


inline void Beam::setBunchingHarmonicOutput(int harm_in){bharm=harm_in;}
inline int Beam::getBunchingHarmonics(){return bharm;}
inline void Beam::set_global_stat(bool in){do_global_stat=in;}
inline bool Beam::get_global_stat(void){return(do_global_stat);}
inline void Beam::setOutput(bool noCurrent_in, bool noEnergy_in, bool noSpatial_in, bool noAux_in) {
  doCurrent = !noCurrent_in;
  doSpatial = !noSpatial_in;
  doEnergy = !noEnergy_in;
  doAux = !noAux_in;
}
inline void Beam::setWriteFilter(bool in_active, int in_from, int in_to, int in_inc) {
   beam_write_filter = in_active;
   beam_write_slices_from = in_from;
   beam_write_slices_to = in_to;
   beam_write_slices_inc = in_inc;
}
inline bool Beam::get_WriteFilter_active() { return beam_write_filter; }
inline int  Beam::get_WriteFilter_from()   { return beam_write_slices_from; }
inline int  Beam::get_WriteFilter_to()     { return beam_write_slices_to; }
inline int  Beam::get_WriteFilter_inc()    { return beam_write_slices_inc; }
#endif
