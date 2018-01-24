#ifndef __GENESIS_CONTROL__
#define __GENESIS_CONTROL__

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string>

#include "mpi.h"
#include "Field.h"
#include "Beam.h"
#include "Undulator.h"
#include "HDF5base.h"
#include "Output.h"


using namespace std;

class Control : public HDF5Base{
 public:
   Control();
   virtual ~Control();
   void applySlippage(double, Field *);
   bool applyMarker(Beam *, vector<Field *> *, Undulator *);
   bool init(int, int, const char *, Beam *, vector<Field *> *, Undulator *,bool,bool);
   void output(Beam *, vector<Field*> *,Undulator *);

 private:
   bool timerun,scanrun,one4one;
   int nslice,ntotal,noffset;
   int rank, size;
   double accushift;
   double sample,reflen,slen;
   int nzout;
   int nwork;
   double *work;
   string root;
};


#endif
