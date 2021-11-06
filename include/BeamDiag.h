#ifndef __GEN_BEAMDIAG_H
#define __GEN_BEAMDIAG_H

#include <iostream>
#include <hdf5.h>

// including "Beam.h" does not give working definition of class Beam,
// as this include file includes this file in turn
class Beam;

class BeamDiag {
public:
	BeamDiag();
	virtual ~BeamDiag();

	virtual void init(int nz, int ns)=0;
	virtual void do_diag(Beam *)=0;
	virtual void output(hid_t parentobj)=0;

	virtual std::string to_str(void) const;
};
std::ostream& operator<< (std::ostream&, const BeamDiag&);

#endif // __GEN_BEAMDIAG_H
