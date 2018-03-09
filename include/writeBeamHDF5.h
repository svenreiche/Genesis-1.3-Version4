#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <complex>


#ifndef __GEN_WRITEBEAMHDF5__
#define __GEN_WRITEBEAMHDF5__

#include "mpi.h"
#include "hdf5.h"
#include "HDF5base.h"
#include "Beam.h"

using namespace std;

class WriteBeamHDF5 : public HDF5Base {
 public:
  WriteBeamHDF5();
  virtual ~WriteBeamHDF5();
  void write(string, Beam *);
 

 private:
  int rank,size;
  hid_t fid;
  void writeGlobal(int,bool,double,double,double,int);

};



#endif

