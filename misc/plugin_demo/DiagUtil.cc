#include <cmath>
#include <cstdlib>
#include <iostream>

#include "DiagUtil.h"
#include "DiagnosticHookS.h"

using namespace std;

bool DiagUtil::verify_datastructure(DiagBeamHookData *pd)
{
	/* verify if version number provided by interface compiled into "GENESIS 1.3" v4 matches */
	if (pd->version!=DIAGFIELD_DATA_STRUCTVERSION) {
		cout << "Plugin: Mismatch of version numbers in data structure: GENESIS provided="
		     << pd->version << " , plugin expects=" << DIAGFIELD_DATA_STRUCTVERSION
		     << ". This is fatal." << endl;

		// crash the program (and produce core dump if user has it enabled)
		abort();
	}
	return(true);
}
bool DiagUtil::verify_datastructure(DiagFieldHookData *pd)
{
	/* verify if version number provided by interface compiled into "GENESIS 1.3" v4 matches */
	if (pd->version!=DIAGFIELD_DATA_STRUCTVERSION) {
		cout << "Plugin: Mismatch of version numbers in data structure: GENESIS provided="
		     << pd->version << " , plugin expects=" << DIAGFIELD_DATA_STRUCTVERSION
		     << ". This is fatal." << endl;

		// crash the program (and produce core dump if user has it enabled)
		abort();
	}
	return(true);
}

int DiagUtil::xlat(int ngrid, int ix, int iy)
{
	int idx = iy*ngrid+ix; // see for instance src/Core/Diagnostic.cpp, function 'DiagField::getTags' (commit id efcc090)

	return(idx);
}

// scale the (sum of) abs^2 of field values to power
double DiagUtil::scale_to_power(DiagFieldHookData *pd, double sum)
{
	// constants from GenMain.cpp
	const double eev = 510999.06;
	const double vacimp = 376.73;

	double ks=4.*asin(1)/pd->xlambda;
	double scl=pd->dgrid*eev/ks;

	// Note 2023-01-07: compiling this and G4 with debug settings (and gcc-9.3), the result of the following computation exhibits sometimes small relative deviations (on the order of 1E-15 relative) from the result in /Field/power:
	// double power = sum * scl*scl/vacimp; // scale to Watt (Diagnostic.cpp, line 627)

	// This gives binary-identical results
	double power = sum;
	power *= scl*scl/vacimp; // scale to Watt (Diagnostic.cpp, line 627)

	return(power);
}

// scale the (sum of) abs^2 of field values to power
double DiagUtil::scale_to_intensity(DiagFieldHookData *pd, double sum)
{
	// constants from GenMain.cpp
	const double eev = 510999.06;
	const double vacimp = 376.73;

	double ks=4.*asin(1)/pd->xlambda;
	double inten = sum * eev*eev/ks/ks/vacimp; // scale to W/m^2 (Diagnostic.cpp, line 632)

	return(inten);
}
