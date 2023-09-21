#include "DiagUtil.h"
#include "DiagnosticHookS.h"

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
