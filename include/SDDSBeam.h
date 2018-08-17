#ifndef __GENESIS_SDDSBEAM__
#define __GENESIS_SDDSBEAM__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <cctype>

#include "StringProcessing.h"
#include "Beam.h"
#include "Setup.h"
#include "Time.h"
#include "Lattice.h"
#include "Sorting.h"
#include "hdf5.h"
#include "Particle.h"
#include "RandomU.h"
#include "ShotNoise.h"
#include "Output.h"

using namespace std;

extern const double eev;
extern const double ce;

class SDDSBeam: public StringProcessing, HDF5Base{
 public:
   SDDSBeam();
   virtual ~SDDSBeam();
   bool init(int,int,map<string,string> *, Beam *, Setup *, Time *, Lattice *);
 private:
   void usage();
   void removeParticles(vector< Particle > *, int);
   void addParticles(vector< Particle > *, int);
   void analyse(double,int);
   void initRandomSeq(int);
   double distance(Particle, Particle);
   RandomU *ran;
   int rank,size;
   double betax,betay,alphax,alphay,charge;
   double xcen,ycen,pxcen,pycen,gamma;
   double ds,matchs0,matchs1,aligns0,aligns1;
   bool center,match,settime;
   int align;
   string file;

   double gavg,xavg,yavg,pxavg,pyavg,xvar,yvar,pxvar,pyvar,xpx,ypy,ex,ey,ax,ay,bx,by;
   vector<double> t,g,x,y,px,py;

};



#endif
