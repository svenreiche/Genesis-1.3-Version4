#ifndef __GENESIS_LOADBEAM__
#define __GENESIS_LOADBEAM__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <complex>

#include "StringProcessing.h"
#include "Setup.h"
#include "Time.h"
#include "Profile.h"
#include "Lattice.h"
#include "QuietLoading.h"
#include "ShotNoise.h"
#include "Beam.h"



using namespace std;

extern const double ce;

class LoadBeam : public StringProcessing{
 public:
   LoadBeam();
   virtual ~LoadBeam();
   bool init(int, int, map<string,string> *,Beam *,Setup *, Time *, Profile *, Lattice *);

 private:
   void usage();
   double gamma,delgam,ex,ey,betax,betay,alphax,alphay;
   double xcen,ycen,pxcen,pycen,current,bunch,bunchphase;
   double emod,emodphase;
   string gammaref,delgamref,exref,eyref,betaxref,betayref;
   string alphaxref,alphayref,xcenref,ycenref,pxcenref,pycenref;
   string currentref,bunchref,bunchphaseref,emodref,emodphaseref; 

};


#endif
