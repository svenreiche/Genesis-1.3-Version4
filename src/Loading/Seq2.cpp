#include <iostream>
#include <vector>
#include <Sequence.h>

using namespace std;

SequenceList::SequenceList(std::vector<double>& vinit) {
	v_ = vinit;
	pos_ = 0;
}

double SequenceList::getElement()
{
	if(pos_<v_.size())
	{
		double t=v_.at(pos_);
		pos_++;
		return(t);
	}

	// data set is exhausted
	cout << "Warning: SequenceList: loaded data is exhausted, using default value" << endl;
	return(default_value_);
}
