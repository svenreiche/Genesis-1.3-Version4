#include <iostream>
//#include <iomanip>
//#include <fstream>
#include <string>
//#include <math.h>
//#include <complex>


#ifndef __GEN_HDF5_CWC__
#define __GEN_HDF5_CWC__

//#include <mpi.h>
#include <hdf5.h>
//#include "HDF5base.h"
//#include "Field.h"

using namespace std;

class HDF5_CollWriteCore {
public:
	HDF5_CollWriteCore();
	~HDF5_CollWriteCore();
	void create_and_prepare(hid_t gid, string dataset, string unit, vector<hsize_t> *totalsize, unsigned long ndim);
	void write(vector<double> *data, vector<hsize_t> *count, vector<hsize_t> *offset);
	void close(void);

private:
	unsigned long ndim_;
	hid_t did_;
	bool did_valid_;
};

#endif

