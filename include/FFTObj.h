#ifndef GENESIS_1_3_VERSION4_FFTOBJ_H
#define GENESIS_1_3_VERSION4_FFTOBJ_H

#include <complex>
#include <map>

#ifdef FFTW
#include <fftw3.h>
#endif

// For field diagnostics: helper class for management of FFT-related resources
#ifdef FFTW
class FFTObj {
public:
	FFTObj(int);
	~FFTObj();

	// The fftw_plan is a pointer, delete copy constructor and copy
	// assignment operator.
	// This avoids dangling pointers in object copies (fftw_plan contains
	// the data locations that were used to prepare the plan).
	FFTObj(const FFTObj&) = delete;
	FFTObj& operator= (const FFTObj&) = delete;

	std::complex<double> *in_ {nullptr};
	std::complex<double> *out_ {nullptr};

	// according to FFTW documentation, fftw_plan is an "opaque pointer type" (https://www.fftw.org/fftw3_doc/Using-Plans.html , 19.01.2023)
	fftw_plan p_;

private:
	int rank_;
	int ngrid_;
};
// map holding FFT-related resources for the values of ngrid
typedef std::map<int,FFTObj *> map_fftobj;
#endif

#endif // GENESIS_1_3_VERSION4_FFTOBJ_H
