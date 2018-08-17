#ifndef __GENESIS_SETUP__
#define __GENESIS_SETUP__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <cctype>

#include "StringProcessing.h"
#include "Lattice.h"

using namespace std;

class Lattice;

class Setup: public StringProcessing{
 public:
   Setup();
   virtual ~Setup();
   bool init(int, map<string,string> *, Lattice *,string,bool);
   //   bool writeGlobal(hid_t, double, int, int,int,int,double, double, double, bool, bool);
   double getReferenceLength();
   void   setReferenceLength(double);
   double getReferenceEnergy();
   double getStepLength();
   void   setStepLength(double);
   bool   getOne4One();
   bool   getShotNoise();
   int    getNpart();
   int    getNbins();
   int    getSeed();
   //   bool   getInputFileNameBeam(string *);
   //   bool   getInputFileNameField(string *);
   bool   getRootName(string *);
   void   setRootName(string *);
   void incrementCount();
   string getLattice();

 private:
   void usage();
   string rootname,lattice,beamline,partfile,fieldfile;
   double gamma0,lambda0,delz;
   bool one4one,shotnoise;
   int seed, rank,npart,nbins,runcount;
};

inline string Setup::getLattice(){return lattice;}
inline double Setup::getReferenceLength(){ return lambda0; }
inline void   Setup::setReferenceLength(double lam){ lambda0=lam; return; }
inline double Setup::getReferenceEnergy(){ return gamma0; }
inline bool   Setup::getOne4One(){return one4one;}
inline bool   Setup::getShotNoise(){return shotnoise;}
inline int    Setup::getNpart(){return npart;}
inline int    Setup::getNbins(){return nbins;}
inline int    Setup::getSeed(){return seed;}
inline double Setup::getStepLength(){return delz;}
inline void   Setup::setStepLength(double din){delz=din;return;}
inline void   Setup::incrementCount(){runcount++; return;}
inline void   Setup::setRootName(string *newname){rootname=*newname; runcount=0; return;}

#endif
