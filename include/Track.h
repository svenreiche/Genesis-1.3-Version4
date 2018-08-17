#ifndef __GENESIS_TRACK__
#define __GENESIS_TRACK__

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <map>
#include <stdlib.h>

#include "hdf5.h"
#include "mpi.h"

#include "Setup.h"
#include "Lattice.h"
#include "AlterLattice.h"
#include "Time.h"
#include "Beam.h"
#include "Field.h"


using namespace std;

class Track{
 public:
   Track();
   virtual ~Track();
   bool init(int, int, map<string,string> *,Beam *, vector<Field *> *,Setup *, Lattice *, AlterLattice *, Time *,bool);
 private:
   void usage();
   double zstop,slen,s0;
   int output_step,dumpFieldStep,dumpBeamStep,sort_step,bunchharm;
   int rank, size;
};


#endif
