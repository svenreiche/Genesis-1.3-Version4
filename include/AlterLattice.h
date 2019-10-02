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
#include "Series.h"

class Setup;
class Lattice;

using namespace std;

class AlterLattice:  public StringProcessing{
 public:
   AlterLattice();
   virtual ~AlterLattice();
   bool init(int, int, map<string,string> *, Lattice *, Setup *, Series *);
 private:
   void usage();
   string element,field;
   double value;
   int instance;
   bool resolve,add;
   double zmatch;
   int rank, size;
   string valueref;
};

#endif
