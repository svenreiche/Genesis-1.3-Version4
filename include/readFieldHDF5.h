#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <complex>


#ifndef __GENESIS_READFIELDHDF5__
#define __GENENIS_READFIELDHDF5__

#include "hdf5.h"
#include "HDF5base.h"
#include "Setup.h"
#include "Time.h"

using namespace std;

class ReadFieldHDF5 : public HDF5Base {
 public:
  ReadFieldHDF5();
  virtual ~ReadFieldHDF5();
  //  bool open(char *, int *, double *, double *);
  bool readGlobal(int, int, string, Setup *,Time *, double, bool);
  void close();



 private:
  hid_t fid;
  double s0,slicelen;
  int  nwork;
  double *work;
};



#endif

