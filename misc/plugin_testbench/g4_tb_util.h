#ifndef __TB_UTIL_H
#define __TB_UTIL_H

#include <fstream>
#include <string>

// class holding the configuration data
// FIXME: currently everything is public, getter functions need to be implemented
class TB_Cfg {
public:
	TB_Cfg();
	bool update_from_stream(std::ifstream&);


	std::string libfile   {"./libdemo.so"};
	std::string parameter {""};

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

private:
	bool update_param(const std::string, const std::string);
	void eat_whitespaces(std::string &);
};

#endif // __TB_UTIL_H
