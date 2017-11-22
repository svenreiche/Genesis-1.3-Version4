#ifndef __GENESIS_UNDULATOR__
#define __GENESIS_UNDULATOR__

#include <iostream>
#include <vector>
#include <math.h>

#include "hdf5.h"

#include "HDF5base.h"
#include "BesselJ.h"

using namespace std;

class Undulator: public HDF5Base{
 public:
   Undulator();
   virtual ~Undulator();

   bool init(hid_t);

   bool advance(int);
   bool inUndulator();
   bool isHelical();
   double steplength();
   double slippage();
   double autophase();
   int outlength();
   bool outstep();
   int getMarker();

   void updateMarker(int, int, int, double);
   void getUndulatorParameters(double *, double *, double *, double *, double *, double *);
   void getQuadrupoleParameters(double *, double *, double *);
   void getCorrectorParameters(double *, double *);
   void getChicaneParameters(double *, double *, double *, double *);
   double getGammaRef();

   double getaw();
   double getku();
   double getz();
   double faw2(double, double); // should be replace in future versions
   double faw(double, double);  
   double fc(int);

   int getStep();

 private: 
   vector<double> aw,ax,ay,ku,kx,ky,cx,cy,gradx,grady;
   vector<double> qf,qx,qy,z,dz,slip,phaseshift; 
   vector<double> chic_angle,chic_lb,chic_ld,chic_lt; 
   vector<double> paw,pkx,pky,pgradx,pgrady,pphase; // perpendicular undulator parameters
   vector<int> helical,marker;
   vector<bool> out;

   double gammaref,zstop;  
   int istepz,nstepz,nout;
};

inline int Undulator::getStep(){
  return istepz;
}

inline int Undulator::getMarker(){
  return marker[istepz]; 
}

inline bool Undulator::isHelical(){
  if (helical[istepz]==0){
    return false;
  }else {
    return true;
  }
}


inline int Undulator::outlength(){
     return nout;
}
inline bool Undulator::outstep(){
    return out[istepz+1];
}

inline double Undulator::getz(){
  return z[istepz]+dz[istepz];
}

inline double Undulator::getaw(){
  return aw[istepz];
}

inline double Undulator::getku(){
  return ku[istepz];
}

inline double Undulator::steplength(){
  return dz[istepz];
}


inline double Undulator::slippage(){
  return slip[istepz];
}

inline double Undulator::autophase(){
  return phaseshift[istepz];
}

inline double Undulator::getGammaRef(){
  return gammaref;
}


inline bool Undulator::inUndulator(){
  if (aw[istepz]>0) {
    return true;
  }
  return false;
}

inline void Undulator::getUndulatorParameters(double *iaw, double *iax, double *iay, double *iku, double *ikx, double *iky)
{
  if (this->inUndulator()){
    *iaw=aw[istepz];
    *iax=ax[istepz];
    *iay=ay[istepz];
    *iku=ku[istepz];
    *ikx=kx[istepz];
    *iky=ky[istepz];
  } else {
    *iaw=0;
    *iax=0;
    *iay=0;
    *iku=0;
    *ikx=0;
    *iky=0;
  }
  return;
}

inline void Undulator::getQuadrupoleParameters(double *iqf, double *iqx, double *iqy)
{
  *iqf=qf[istepz];
  if (*iqf==0) {
    *iqx=0;
    *iqy=0;
  } else {
    *iqx=qx[istepz];
    *iqy=qy[istepz];
  }
  return;
}


inline void Undulator::getCorrectorParameters(double *icx, double *icy)
{
  *icx=cx[istepz];
  *icy=cy[istepz];
  return;
}


inline void Undulator::getChicaneParameters(double *iangle, double *ilb, double *ild, double *ilt)
{
  *iangle=chic_angle[istepz];
  if (*iangle==0){
    *ilb=0;
    *ild=0;
    *ilt=0;
  } else{
    *ild=chic_ld[istepz];
    *ilb=chic_lb[istepz];
    *ilt=chic_lt[istepz];
  }

  return;
}
#endif
