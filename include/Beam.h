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
#include "BeamDiag.h"

using namespace std;

extern const double ce;

class BeamDiag_Std;

class Beam{
 public:
   Beam();
   virtual ~Beam();
   void initDiagnostics(int);
   void diagnostics(bool,double);
   // void diagnosticsStart();
   void init(int, int, double,double, double,bool);
   void initSorting(int,int,bool,bool);
   void initEField(double,int,int,int,double);
   void initIncoherent(int, int, bool,bool);
   void initWake(unsigned int, unsigned int, double, double *, double *, double *,double *, double,double, bool);
   bool harmonicConversion(int,bool);
   bool subharmonicConversion(int,bool);
   int sort();
   void track(double, vector<Field *> *, Undulator *);

   void setWriteFilter(bool,int,int,int);
   bool get_WriteFilter_active();
   int  get_WriteFilter_from();
   int  get_WriteFilter_to();
   int  get_WriteFilter_inc();


   void register_beam_diag(BeamDiag *);
   void clear_beam_diag(void);
   void beam_diag_store_results(hid_t);
   void beam_diag_list_registered(void);
   BeamDiag_Std *bd_std; // this pointer is currently needed for first call to function replacing Beam::diagnosticsStart from Control::init

   vector< vector<Particle> > beam;
   vector<double> current,eloss;
   double reflength,slicelength;   // for conversion of theta in Particle to real position
   double s0;         // starting position of the time-window
   bool one4one;     // flag whether one4one simulation is done
   int nbins;

private:
   BeamSolver solver;
   Incoherent incoherent;
   Collective col;
   Sorting sorting;
   int idx;

   bool beam_write_filter;
   int beam_write_slices_from, beam_write_slices_to, beam_write_slices_inc;

   vector<BeamDiag *> diaghooks;
   bool can_change_diaghooks;
   void beam_diag_do_diag(double);
};

inline void Beam::initIncoherent(int base, int rank, bool spread, bool loss){
  incoherent.init(base,rank,spread,loss);
  return;
}

inline void Beam::initEField(double rmax, int ngrid, int nz, int nphi, double lambda){
  solver.initEField(rmax,ngrid,nz,nphi,lambda);
  return;
}

inline void Beam::initWake(unsigned int ns, unsigned int nsNode, double ds, double *wakeext, double *wakeres, double *wakegeo,double *wakerou, double ztrans, double radius, bool transient){
  col.initWake(ns, nsNode, ds, wakeext, wakeres, wakegeo, wakerou, ztrans, radius, transient);
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
