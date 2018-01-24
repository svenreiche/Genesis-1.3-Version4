#ifndef __GENESIS_OUTPUT__
#define __GENESIS_OUTPUT__

#include <vector>
#include <iostream>
#include <string>
#include <math.h>

#include "hdf5.h"
#include "HDF5base.h"
#include "Particle.h"
#include "Beam.h"
#include "Field.h"
#include "Undulator.h"

using namespace std;

extern const double vacimp;
extern const double eev;

class Output : public HDF5Base {
 public:
   Output();
   virtual ~Output();
   void open(string,int,int);
   void close();
   void writeFieldBuffer(Field *);
   void writeBeamBuffer(Beam *);
   void writeLattice(Beam *, Undulator *);
   void writeGlobal(double,double,double,double,bool,bool,bool);

 private:
   void write(hsize_t,string,string,double *);
   void writeBuffer(hid_t,string,vector<double> *);
   void writeSingleNode(hid_t,string,vector<double> *);

   hid_t fid;   
   hsize_t s0, ds;
};


#endif
