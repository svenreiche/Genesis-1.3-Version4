#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <complex>


#ifndef __GEN_WRITEFIELDHDF5__
#define __GEN_WRITEFIELDHDF5__

#include "mpi.h"
#include "hdf5.h"
#include "HDF5base.h"
#include "Field.h"

using namespace std;


extern const double vacimp;
extern const double eev;

class WriteFieldHDF5 : public HDF5Base {
 public:
  WriteFieldHDF5();
  virtual ~WriteFieldHDF5();
  void write(string fileroot, vector<Field *> *field);


 private:
  void writeMain(string fileroot, Field *field);
  void writeGlobal(double, double, double, double, int, int);
  hid_t fid;
  int rank, size;
};



#endif

