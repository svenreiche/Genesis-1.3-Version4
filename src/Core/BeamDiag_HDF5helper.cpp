#include "BeamDiag_HDF5helper.h"
#include "HDF5base.h"

BeamDiag_HDF5Helper::BeamDiag_HDF5Helper() {
	/* avoid random garbage values in these most important variables controlling HDF5 I/O */
	s0 = 1;
	ds = 1;
}
BeamDiag_HDF5Helper::~BeamDiag_HDF5Helper() {}

void BeamDiag_HDF5Helper::set_s0(int s0_in) {
	s0 = s0_in;
}
void BeamDiag_HDF5Helper::set_ds(int ds_in) {
	ds = ds_in;
}
