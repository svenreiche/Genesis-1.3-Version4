#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "FieldManipulator.h"

// #define DBG_SPP

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
	if(do_spp) {
		FieldManipulator_SPP_Params p;

		p.harm = harm;
		p.spp_l = spp_l;
		p.spp_nsect = spp_nsect;
		p.spp_phi0 = spp_phi0;
		apply_SPP(p_fld, time, p); // harm, spp_l, spp_nsect, spp_phi0);
	}

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


bool FieldManipulator::apply_SPP(Field *p_fld, Time *time, FieldManipulator_SPP_Params &p)
{
	bool dbg_dump_phi=true;

	if(rank_==0) {
		cout << "Applying SPP to field harmonic=" << p.harm
		     << ", parameters: spp_l=" << p.spp_l << ", spp_nsect=" << p.spp_nsect << ", spp_phi0=" << p.spp_phi0 << endl;
	}

	const int ngrid = p_fld->ngrid;

	const int icenter = (ngrid-1)/2;
	assert((ngrid%2)==1); // In GENESIS, ngrid is always odd. Just to be sure.

	int iy, ix;
	for(int idslice=0; idslice < time->getNodeNSlice(); idslice++)
	{
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
				double phi = apply_SPP_getphase(dix, diy, p);
				
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

double FieldManipulator::apply_SPP_getphase(int dix, int diy, FieldManipulator_SPP_Params &p)
{
	const double twopi = 8.*atan(1.0);
	double phi=0.;

#ifdef DBG_SPP
	// test case to check indexing and phase conventions
	phi = 0;
	if((dix<=0) && (diy<=0))
		phi = 2.*atan(1.);
	return(phi);
#endif

	if (0==p.spp_l)
		return(p.spp_phi0);



	if ((dix!=0) || (diy!=0))
		phi = atan2(diy,dix); // if both arguments to atan2 are zero this would give a 'domain error'

	// sectorize SPP, if requested
	if(p.spp_nsect>0) {
		/*
		 * Range of 'atan2' return values to expect:
		 * For diy=0 and dix<0, we expect atan2(diy,dix)=+pi.
		 * For diy=-1 and dix<0, we expect atan2(diy,dix)>-pi.
		 */
		const double phi_in = phi;
		double tmp_phi = (0.5+phi/twopi);
		tmp_phi *= p.spp_nsect;
		// tmp_phi -= 1e-6; // subtract small number, shifting the boundary between neighboring sectors. This makes symmetric cases (such as nsect=4) look more "natural"

		tmp_phi += p.spp_phi0*p.spp_nsect/twopi; // FIXME: divide phase change by spp_l here?

		/* (i) round */
		tmp_phi = floor(tmp_phi);

		/* (ii) enforce value range (phase offsets; in addition, there could be rounding issues with 'floor' at the edges of the range) */
		while(tmp_phi<0) {
			tmp_phi += p.spp_nsect; // periodic
		}
		while(tmp_phi>p.spp_nsect) {
			tmp_phi -= p.spp_nsect; // periodic
		}

		// transform back to original range -pi..+pi
		tmp_phi /= p.spp_nsect;
		phi = twopi * (tmp_phi-0.5);

		phi *= p.spp_l;
		return(phi);
	}


	// continuous SPP
	phi = p.spp_l*phi + p.spp_phi0;
	return(phi);
}
