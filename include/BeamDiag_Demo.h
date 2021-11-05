#ifndef __GEN_BEAMDIAG_DEMO

#include <vector>

#include "BeamDiag.h"
#include "Beam.h"

using namespace std;

// Note that class HDF5Base adds lots of objects
//                                    vvvvvvvvvvvvvvv
class BeamDiag_Demo: public BeamDiag, public HDF5Base {
public:
	BeamDiag_Demo();
	~BeamDiag_Demo();

	void init(int nz, int ns);
	void do_diag(Beam *);
	void output(hid_t parentobj);

	// functions specific to this BeamDiag implementation
	void config(unsigned long long);
	void set_verbose(bool);

private:
	vector<unsigned long long> demodiag_data_;

	int ns_;
	int idx_;

	bool is_initialized_;
	bool is_configured_;
	unsigned long long param_x_;

	bool verbose_;

	int my_rank_;
};

#endif // __GEN_BEAMDIAG_DEMO
