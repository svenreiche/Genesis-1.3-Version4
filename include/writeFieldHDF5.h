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

class WriteFieldHDF5 : public HDF5Base {
 public:
  WriteFieldHDF5();
  virtual ~WriteFieldHDF5();
  void writeMain(string fileroot, Field *field,int);
  void write(string fileroot, vector<Field *> *field);
  void writeGlobal(hid_t,double,double,double,double,int);
  int writeSlice(hid_t,Field *field,int);

 private:
  int rank, size;
};



#endif

