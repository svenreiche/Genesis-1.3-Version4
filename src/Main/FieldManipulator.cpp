#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "FieldManipulator.h"

FieldManipulator::FieldManipulator() {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
	MPI_Comm_size(MPI_COMM_WORLD, &size_);
}
FieldManipulator::~FieldManipulator() { }


void FieldManipulator::usage(){
  cout << "List of keywords for FIELD_MANIPULATOR" << endl;
  cout << "&field_manipulator" << endl;
  cout << " int harm = 1" << endl;
  cout << " double scale_power = 1.0" << endl;
  cout << " double spp_l = 0.0" << endl;
  cout << " int spp_nsect = 0" << endl;
  cout << " double spp_phi0 = 0.0" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool FieldManipulator::init(int rank, int size, map<string,string> *arg, vector<Field *> *fieldin,  Setup *setup, Time *time, Profile *prof)
{
	map<string,string>::iterator end=arg->end();
	int harm = 1;

	double scale_P = 1.0;
	bool do_scale = false;

	double spp_l    = 0.0;
	int spp_nsect   = 0;
	double spp_phi0 = 0.0;
	bool do_spp = false;

	// extract parameters
	if (arg->find("harm")!=end)         {harm    = atoi(arg->at("harm").c_str());        arg->erase(arg->find("harm"));}
	if (arg->find("scale_power")!=end) {
		do_scale = true; // specifying the parameter switches on code for scaling of the field
		scale_P = atof(arg->at("scale_power").c_str());
		arg->erase(arg->find("scale_power"));
	}

	if (arg->find("spp_nsect")!=end) {
		do_spp = true; // specifying an SPP parameter switches the corresponding code block
		spp_nsect = atoi(arg->at("spp_nsect").c_str());
		arg->erase(arg->find("spp_nsect"));
	}
	if (arg->find("spp_l")!=end) {
		do_spp = true; // specifying an SPP parameter switches the corresponding code block
		spp_l = atof(arg->at("spp_l").c_str());
		arg->erase(arg->find("spp_l"));
	}
	if (arg->find("spp_phi0")!=end) {
		do_spp = true; // specifying an SPP parameter switches the corresponding code block
		spp_phi0 = atof(arg->at("spp_phi0").c_str());
		arg->erase(arg->find("spp_phi0"));
	}



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


	/*** process field ***/
	if(do_scale)
		scale(p_fld, time, harm, scale_P);
	if(do_spp)
		apply_SPP(p_fld, time, harm, spp_l, spp_nsect, spp_phi0);

	return(true);
}

/*****************************************************/
/*** FUNCTIONS PERFORMING THE ACTUAL MANIPULATIONS ***/
/*****************************************************/

bool FieldManipulator::scale(Field *p_fld, Time *time, int harm, double scale_P)
{
	if(rank_==0) {
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


bool FieldManipulator::apply_SPP(Field *p_fld, Time *time, int harm,
	double spp_l, int spp_nsect, double spp_phi0)
{
	bool dbg_dump_phi=true;

	if(rank_==0) {
		cout << "Applying SPP to field harmonic=" << harm
		     << ", parameters: spp_l=" << spp_l << ", spp_nsect=" << spp_nsect << ", spp_phi0=" << spp_phi0 << endl;
	}

	const int ngrid = p_fld->ngrid;

	const int icenter = (ngrid-1)/2;
	assert((ngrid%2)==1); // In GENESIS, ngrid is always odd. Just to be sure.

	int iy, ix;
	for(int idslice=0; idslice < time->getNodeNSlice(); idslice++) {
		bool do_dump = dbg_dump_phi && (rank_==0) && (idslice==0);

		ofstream ofs;

		if(do_dump) {
			ofs.open("dump_phi.txt", ofstream::out);
		}

		for(iy=0; iy<ngrid; iy++) {
			for(ix=0; ix<ngrid; ix++) {
				int idx = iy*ngrid+ix; // see for instance src/Core/Diagnostic.cpp, function 'DiagField::getTags' (commit id efcc090)
				
				int dix = ix-icenter;
				int diy = iy-icenter;
				
				// determine SPP phase for current grid cell
				double phi=0.;
				if ((dix!=0) || (diy!=0))
					phi += atan2(diy,dix); // if both arguments to atan2 are zero this would give a 'domain error'
				phi = spp_l*phi + spp_phi0;
				
				// apply phase factor to complex field in current grid cell
				complex<double> phase_factor = polar(1., phi);
				p_fld->field[idslice].at(idx) *= phase_factor;

				if(do_dump)
					ofs << phi << " ";
			}
		}

		if(do_dump) {
			ofs.close();
		}
	}

	return(true);
}
