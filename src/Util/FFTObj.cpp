#include <iostream>
#include <complex>
#include <map>
#include <mpi.h>

#ifdef FFTW
#include <fftw3.h>
#endif

#include "FFTObj.h"

using namespace std;

FFTObj::FFTObj(int ngrid)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

	// this is severe error: print on *all* ranks
	if(ngrid<0) {
		cerr << "Error: enforcing ngrid>0" << endl;
		ngrid=1;
	}
	in_ = new complex<double>[ngrid * ngrid];
	out_ = new complex<double>[ngrid * ngrid];
	p_ = fftw_plan_dft_2d(ngrid, ngrid,
	    reinterpret_cast<fftw_complex *>(in_),
	    reinterpret_cast<fftw_complex *>(out_),
	    FFTW_FORWARD, FFTW_ESTIMATE /* FFTW_MEASURE */);
	ngrid_ = ngrid;
}

FFTObj::~FFTObj() {
	if(rank_==0) {
		cout << "~FFTObj" << endl;
	}

	fftw_destroy_plan(p_);
	delete [] in_;
	delete [] out_;
}

