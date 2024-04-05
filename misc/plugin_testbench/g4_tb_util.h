#ifndef __TB_UTIL_H
#define __TB_UTIL_H

// class holding the configuration data
// FIXME: currently everything is public, getter functions need to be implemented
class TB_Cfg {
public:
	TB_Cfg();

	int nslice    {4}; // !number of slices *per* process on the MPI communicator!
	int nz        {3};
	double dgrid  {100e-6};
	int ngrid     {151};
	double lambda {1e-10};
	int sample    {1};
	double s0     {0};
	int harm      {1};

	double power  {1e6};
	double w0     {20e-6};
};

#endif // __TB_UTIL_H
