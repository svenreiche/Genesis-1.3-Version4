#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <complex>


#ifndef __GEN_WRITEBEAMHDF5__
#define __GEN_WRITEBEAMHDF5__

#include "hdf5.h"
#include "HDF5base.h"


using namespace std;

class WriteBeamHDF5 : public HDF5Base {
 public:
  WriteBeamHDF5();
  virtual ~WriteBeamHDF5();
  void open(string);
  void close();
  void writeSlice(int, int, complex< double > *);

 private:
  hid_t fid;
  int  nwork;
  double *work;
};



#endif

