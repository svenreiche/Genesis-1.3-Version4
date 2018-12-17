#ifndef __GENESIS_ALTERSETUP__
#define __GENESIS_ALTERSETUP__

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
#include "Time.h"
#include "Setup.h"
#include "Beam.h"
#include "Field.h"


using namespace std;

class Lattice;

class AlterSetup: public StringProcessing, public HDF5Base {
 public:
   AlterSetup();
   virtual ~AlterSetup();
   bool init(int, map<string,string> *, Setup *, Lattice *, Time *, Beam *,vector<Field *> *);


 private:
   void usage();
   string rootname,lattice,beamline;
   double delz;
   bool resample,disable;
   int harmonic,subharmonic,rank;
};


#endif
