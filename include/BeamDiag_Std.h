#ifndef __GEN_BEAMDIAG_STD
#define __GEN_BEAMDIAG_STD

#include <vector>

#include "BeamDiag.h"
#include "Beam.h"

using namespace std;

// Note that class HDF5Base adds lots of objects
//                                    vvvvvvvvvvvvvvv
class BeamDiag_Std: public BeamDiag {
public:
	BeamDiag_Std();
	~BeamDiag_Std();

	void init(int nz, int ns);
	void do_diag(const Beam *, double);
	void output(hid_t parentobj);
	std::string to_str(void) const;

	// functions specific to this BeamDiag implementation
	void do_initial_diag(const Beam *);

	void setBunchingHarmonicOutput(int);
	int getBunchingHarmonics(void);
	void set_global_stat(bool);
	bool get_global_stat(void);
	void setOutput(bool, bool, bool, bool);


	/* former member of class Beam (only member that is 'public' because Output::writeLattice needs access) */
	vector<double> zpos;

private:
	/* former members of class Beam */
	// output buffer
	vector<double> gavg,gsig,xavg,xsig,yavg,ysig,pxavg,pyavg,bunch,bphi,efld;
	vector<double> bx,by,ax,ay,ex,ey,cu;
	//   vector<unsigned long long> partcount;
	vector< vector<double> > bh,ph;  // harmonic bunching and bunching phase
	//global values
	vector<double> tgavg, tgsig, txavg,txsig,tyavg, tysig,tbun;  // global values, averaging over the entire beam 


   int bharm;
   bool do_global_stat;
   bool doCurrent, doSpatial, doEnergy, doAux;



	int ns_;
	int idx_;

	bool is_initialized_;
	bool is_configured_;
	bool verbose_;

	int my_rank_;
};

#endif // __GEN_BEAMDIAG_STD
