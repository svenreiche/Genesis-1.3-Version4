#ifndef __GENESIS_FIELDMANIPULATOR_H
#define __GENESIS_FIELDMANIPULATOR_H

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <complex>

#include "StringProcessing.h"
#include "Setup.h"
#include "GenTime.h"
#include "GenProfile.h"
#include "GaussHermite.h"
#include "Field.h"

using namespace std;

class FieldManipulator: public StringProcessing {
public:
	FieldManipulator();
	~FieldManipulator();

	bool init(int, int, map<string,string> *,vector<Field *> *, Setup *, Time *, Profile *);

private:
	void usage(void);

	// functions performing field manipulations
	bool scale(Field *, Time *, int, double);
	bool apply_SPP(Field *, Time *, int, double, int, double);

	// MPI related members
	int rank_ {0};
	int size_ {1};
};

#endif // __GENESIS_FIELDMANIPULATOR_H
