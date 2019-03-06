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

extern const double vacimp;
extern const double eev;

class ReadFieldHDF5 : public HDF5Base {
 public:
  ReadFieldHDF5();
  virtual ~ReadFieldHDF5();
  bool readGlobal(int, int, string, Setup *,Time *, int, bool);
  bool readSlice(double, vector<complex<double> >*);
  void close();
  int getNGrid();
  double getDGrid();


 private:
  hid_t fid;
  double s0,slicelen,slen,dgrid,scl;
  int  nwork,count,ngrid;
  double *work;
};

inline int ReadFieldHDF5::getNGrid(){
  return ngrid;
}

inline double ReadFieldHDF5::getDGrid(){
  return dgrid;
}

#endif

