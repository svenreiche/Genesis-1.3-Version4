#ifndef __GENESIS_ALTERLATTICE__
#define __GENESIS_ALTERLATTICE__

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <map>
#include <stdlib.h>

#include "Lattice.h"
#include "Setup.h"

class Setup;
class Lattice;

using namespace std;

class AlterLattice{
 public:
   AlterLattice();
   virtual ~AlterLattice();
   bool init(int, int, map<string,string> *, Lattice *, Setup *);
 private:
   void usage();
   double err_aw,err_ax,err_ay,err_qx,err_qy;
   double zmatch;
   int rank, size,nlat;
};


#endif
