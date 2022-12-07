#include <cmath>
#include "FieldManipulator.h"

FieldManipulator::FieldManipulator() { }
FieldManipulator::~FieldManipulator() { }


void FieldManipulator::usage(){
  cout << "List of keywords for FIELD_MANIPULATOR" << endl;
  cout << "&field_manipulator" << endl;
  cout << " int harm = 1" << endl;
  cout << " double scale_power = 1.0" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool FieldManipulator::init(int rank, int size, map<string,string> *arg, vector<Field *> *fieldin,  Setup *setup, Time *time, Profile *prof)
{
	map<string,string>::iterator end=arg->end();
	int harm = 1;
	double scale_P = 1.0;

	// extract parameters
	if (arg->find("harm")!=end)         {harm    = atoi(arg->at("harm").c_str());        arg->erase(arg->find("harm"));}
	if (arg->find("scale_power")!=end)  {scale_P = atof(arg->at("scale_power").c_str()); arg->erase(arg->find("scale_power"));}

	if (arg->size()!=0){
		if (rank==0) {
			cout << "*** Error: Unknown elements in &field_manipulator" << endl;
			this->usage();
		}
		return(false);
	}


	// check parameters
	if(!(harm>0)) {
		if (rank==0) {
			cout << "*** Error in field_manipulator: harmonic must be positive" << endl;
		}
		return(false);
	}
	if(scale_P<0) {
		if (rank==0) {
			cout << "*** Error in field_manipulator: scale_power must not be negative" << endl;
		}
		return(false);
	}


	// identify the harmonic field to be manipulated
	// (and stop if there is no field for requested harmonic)
	Field *p_fld=NULL;
	for (int k=0; k<fieldin->size(); k++) {
		if(fieldin->at(k)->harm == harm) {
			p_fld = fieldin->at(k);
			break;
		}
	}
	if(NULL==p_fld) {
		if (rank==0) {
			cout << "*** Error in field_manipulator: there is no field for requested harmonic " << harm << endl;
		}
		return(false);
	}


	// *** process field ***
	if(rank==0) {
		cout << "Scaling power in field harmonic=" << harm << " by factor " << scale_P << " ..." << endl;
	}
	int ngrid = p_fld->ngrid;
	double scale_field = sqrt(scale_P);
	for(int idslice=0; idslice < time->getNodeNSlice(); idslice++) {
		for(int k=0; k<ngrid*ngrid; k++) {
			p_fld->field[idslice].at(k) *= scale_field;
		}
	}

	return(true);
}
