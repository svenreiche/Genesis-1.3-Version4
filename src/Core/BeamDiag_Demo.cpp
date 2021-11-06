#include <mpi.h>
#include <hdf5.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include "BeamDiag_Demo.h"

BeamDiag_Demo::BeamDiag_Demo() {
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);

	// force call to class-specific function 'init' before 'do_diag'
	is_initialized_=false;
	// force call to class-specific function 'configure' before 'do_diag'
	is_configured_=false;

	// these are overwritten before the actual tracking begins
	param_x_ = 0;

	verbose_ = false;

	ns_=0;
	idx_=0;
}

BeamDiag_Demo::~BeamDiag_Demo() {

}

std::string BeamDiag_Demo::to_str(void) const {
	stringstream ss;
	ss << "BeamDiag_Demo(verbose=" << verbose_ << ")";
	return(ss.str());
}

void BeamDiag_Demo::init(int nz, int ns) {
	demodiag_data_.resize(nz*ns);
	ns_ = ns;
	idx_ = 0;
	is_initialized_=true;
}

/* configuration functions */
void BeamDiag_Demo::config(unsigned long long param_x_in) {
	param_x_ = param_x_in;

	is_configured_=true;
}
void BeamDiag_Demo::set_verbose(bool v_in) {
	verbose_ = v_in;
}

/* Called when the electron beam is diagnosed, typically after every integration step */
void BeamDiag_Demo::do_diag(Beam *beam) {
	int ioff=idx_*ns_;

	if((is_configured_==false) || (is_initialized_==false)){
		cout << "Error (BeamDiag_Demo): do_diag called, but class instance is not configured" << endl;
	}

	if(verbose_ && (0==my_rank_)) {
		cout << "--> BeamDiag_Demo::do_diag called" << endl;
	}

	/* loop over local slices */	
	for (int is=0; is<ns_; is++) {
		// obtain number of particles in this slice
		unsigned int npart = beam->beam.at(is).size();

		/* here you could process the particles */

		demodiag_data_.at(ioff+is) = npart+param_x_;
	}

	idx_++;
}

/* This function is called when the .out.h5 file is written */
void BeamDiag_Demo::output(hid_t parentobj) {
	if(verbose_ && (0==my_rank_)) {
		cout << "--> BeamDiag_Demo::output called" << endl;
	}

	// Configure HDF5 write using members inherited from HDF5Base class (FIXME)
	ds = ns_;
	s0 = my_rank_*ns_;
	writeBufferULL(parentobj, "diag_demo", " ", &demodiag_data_);
}
