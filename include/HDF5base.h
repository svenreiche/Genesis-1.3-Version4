
#ifndef __GEN_HDF5BASE__
#define __GEN_HDF5BASE__

#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "hdf5.h"


using namespace std;

class HDF5Base{
 public:
  HDF5Base();
  virtual ~HDF5Base();
 protected:
  bool isOpen;
  int nwork;
  double *work;
  int s0;
  hsize_t ds; // size in s for 2D array inwrite buffer

  
  // function from output and wrideHDF5beam/field
  void writeBuffer(hid_t,string,string,vector<double> *);
  void writeSingleNode(hid_t,string,string,vector<double> *);
  void writeSingleNodeString(hid_t,string, string *);
  void writeSingleNodeInt(hid_t, string,vector<int> *);



  void createExpDataset(hid_t fid, char *name, hsize_t nz, hsize_t ns);
  void expandDataset(hid_t fid, vector<double> *rec, int pos, hsize_t recsize, hsize_t slice, char *name);
  void writeDataDouble(hid_t fid, const char *name, const double *data, int size);
  void writeDataInt(hid_t fid, const char *name, const int *data, int size);
  void writeDataChar(hid_t fid, const char *name, const char *data, int size);
  void readDataDouble(hid_t fid, char *name, double *data, int size);
  void readDataChar(hid_t fid, char *name, char *data, int size);
  void readDataInt(hid_t fid, char *name, int *data, int size);
  int getDatasetSize(hid_t fid, char *name);
  bool checkForLink(hid_t fid, string name);



  void writeDouble1DExist(hsize_t datsize, double *data, hid_t,  string);
  void writeDouble1D(hsize_t reclen, hsize_t datsize, double *data, hid_t,  string);
  void writeInt1D(hsize_t reclen, hsize_t datsize, int *data, hid_t,  string);
  void writeChar1D(hsize_t reclen, hsize_t datsize, const char *data, hid_t,  string);

  void readDouble1D(hid_t fid, const char *name, double *data, hsize_t dat, hsize_t chunk);

  bool simpleReadDouble1D(const string &path, vector<double> *data);
  bool browseFile(const string &path, vector<string> *);

};

extern "C"  herr_t file_info(hid_t loc_id,const char *name, const H5L_info_t *linfo, void *opdata); 


#endif

