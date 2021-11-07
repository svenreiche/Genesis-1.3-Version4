#ifndef __GEN4_BEAMDIAG_HDF5HELPER
#define __GEN4_BEAMDIAG_HDF5HELPER

#include "HDF5base.h"

using namespace std;

class BeamDiag_HDF5Helper: public HDF5Base {
public:
	BeamDiag_HDF5Helper();
	~BeamDiag_HDF5Helper();

	void set_s0(int);
	void set_ds(int);

	using HDF5Base::writeSingleNode;
	using HDF5Base::writeBuffer;
};

#endif //  __GEN4_BEAMDIAG_HDF5HELPER
