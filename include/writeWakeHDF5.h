#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <complex>


#ifndef __GEN_WRITEWAKEHDF5__
#define __GEN_WRITEWAKEHDF5__

#include <mpi.h>
#include <hdf5.h>
#include "HDF5base.h"

using namespace std;

class WriteWakeHDF5 : public HDF5Base {
 public:
  WriteWakeHDF5();
  virtual ~WriteWakeHDF5();
  bool write(string, double *, double *, double *, double *, unsigned long);

 private:
  hid_t fid;
};



#endif

