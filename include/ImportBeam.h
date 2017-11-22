#ifndef __GENESIS_IMPORTBEAM__
#define __GENESIS_IMPORTBEAM__

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

class ImportBeam : public StringProcessing{
 public:
   ImportBeam();
   virtual ~ImportBeam();
   bool init(int, int, map<string,string> *,Beam *,Setup *, Time *);

 private:
   void usage();
   string file;
   double offset;
   bool dotime;
};


#endif
