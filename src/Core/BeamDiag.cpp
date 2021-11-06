#include "BeamDiag.h"

BeamDiag::BeamDiag() {
}

BeamDiag::~BeamDiag() {
}

std::string BeamDiag::to_str(void) const {
	return "BeamDiag object (consider implementing a 'to_str' member function)";
}

// This operator is not member of the class, see:
// Item 24 in S. Meyers: Effective C++, 3rd Ed., Pearson Education, 2005
std::ostream& operator<< (std::ostream& os, const BeamDiag& bd) {
	return (os << bd.to_str());
}
