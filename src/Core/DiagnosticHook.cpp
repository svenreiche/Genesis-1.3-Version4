/*
 * Interface for plugins (field diagnostics)
 *
 * Initial code by C. Lechner, EuXFEL
 */

#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <mpi.h>

#include "DiagnosticHook.h"
#include "DiagnosticPlugin.h"
#include "DiagnosticHookS.h"

DiagFieldHook::DiagFieldHook()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size_);
}

DiagFieldHook::~DiagFieldHook()
{
	if(my_rank_==0) {
		cout << "DiagFieldHook::~DiagFieldHook()" << endl;
	}

	li_.close_lib();
}

bool DiagFieldHook::init(DiagFieldPluginCfg *pin)
{
	if(my_rank_==0) {
		cout << "DiagFieldHook::init" << endl;
	}

	/* Copy all needed infos from the configuration data */
	// li_.libfile_           = pin->libfile;
	//li_.obj_prefix_        = pin->obj_prefix;
	li_.parameter_         = pin->parameter;
	// li_.lib_verbose_       = pin->lib_verbose;
	// li_.interface_verbose_ = pin->interface_verbose;

	obj_prefix_        = pin->obj_prefix;
	// parameter_         = pin->parameter;
	lib_verbose_       = pin->lib_verbose;
	interface_verbose_ = pin->interface_verbose;
	
	return(li_.init_lib(pin->libfile));
}

void DiagFieldHook::set_runid(int runid_in)
{
	runid_ = runid_in;
}
const std::string& DiagFieldHook::get_info_txt() const
{
	return(li_.get_info_txt());
}

bool DiagFieldHook::update_data(std::map<std::string,std::vector<double> > &val, string key, size_t idx, double v)
{
	if (val.find(key)==val.end())
		return(false);

	val[key].at(idx) = v;
	return(true);
}





/***************************************/
/***** INTERFACE TO GENESIS 1.3 v4 *****/
/***************************************/

std::map<std::string,OutputInfo> DiagFieldHook::getTags(FilterDiagnostics &filter_in)
{
	// filter_in configures what is active and is provided by the caller

	tags.clear();
	filter.clear();

#if 0
	/*
	*  register resource that will appear in the .out.h5 file under
	*  /Field/hook, /Field3/hook, etc.
	*/
	tags["hook"] = {false,false," " /* <- if no unit, you must provide a SPACE */};
#endif

	for(auto const& the_obj: li_.obj_names_) {
		stringstream ss;
		ss << obj_prefix_;
		ss << "/";
		ss << the_obj;
		tags[ss.str()] = {false,false," " /* <- if no unit, you must provide a SPACE */};
	}

	return tags; // must be returned
}

void DiagFieldHook::getValues_worker(Field *field, std::map<std::string,std::vector<double> >&val, int iz)
{
	// remark: data management from 'DiagFieldUser::getValues'

	const int ns = field->field.size();
	int is0 = 0;

	// loop over field
	for (auto const &slice: field->field) {
		vector<double> dataout(li_.obj_names_.size(),-1.);

		int is = (ns + is0 - field->first) % ns;

		DiagFieldHookData hd;
		hd.version = DIAGFIELD_DATA_STRUCTVERSION; // version number that is incremented when structure layout changes
		hd.verbose = lib_verbose_;
		hd.mpi_rank = my_rank_;
		hd.mpi_size = comm_size_;
		hd.runid = runid_;
		hd.iz = iz;
		hd.is = is;
		hd.ns = ns;
		hd.harm = field->harm;
		hd.ngrid = field->ngrid;
		hd.dgrid = field->dgrid;
		hd.xlambda = field->xlambda;
		// C++11 array-oriented access, see for instance https://en.cppreference.com/w/cpp/numeric/complex
		// hd.datain = reinterpret_cast<const double *>(&slice);
		hd.datain = &slice;
		hd.dataout = &dataout;
		hd.do_multi = false;

		/* execute the code in the library */
		li_.pdiagfield_->doit(&hd);
		if(interface_verbose_ && (my_rank_==0)) {
			cout << "Rank 0: is=" << is << ", data from plugin=(";
			for(int kk=0; kk<dataout.size(); kk++) {
				cout << dataout[kk];
				if((1+kk)<dataout.size())
					cout << ", ";
			}
			cout << ")" << endl;
		}


		/* save the result into the provided array */
		int idx = iz*ns+is;         // compute index for saving the data
		for(int kk=0; kk<li_.obj_names_.size(); kk++) {
			stringstream ss;
			ss << obj_prefix_;
			ss << "/";
			ss << li_.obj_names_[kk];
			update_data(val, ss.str(), idx, dataout.at(kk));
		}

		// increase the slice counter
		is0++;
	}
}

void DiagFieldHook::getValues_multiworker(Field *field, std::map<std::string,std::vector<double> >&val, int iz)
{
	const int ns = field->field.size();
	vector<const vector<complex <double> > *> multi_datain(ns, nullptr);
	vector<vector<double> > multi_dataout(ns);
	for(int k=0; k<ns; k++)
		multi_dataout.at(k).resize(li_.obj_names_.size(),-1.);


	/* reorganize the slices */
	int is0 = 0;
	for (auto const &slice: field->field) {
		int is = (ns + is0 - field->first) % ns;
		
		multi_datain.at(is) = &slice;

		// increase the slice counter
		is0++;
	}
	
	/* sanity check that would only fail if two slices have the same local slice number (must never happen) */
	bool ok = all_of(multi_datain.begin(), multi_datain.end(), [](auto *p){return(p!=nullptr);});
	if(!ok) {
		abort();
	}
	
	DiagFieldHookData hd;
	hd.version = DIAGFIELD_DATA_STRUCTVERSION; // version number that is incremented when structure layout changes
	hd.verbose = lib_verbose_;
	hd.mpi_rank = my_rank_;
	hd.mpi_size = comm_size_;
	hd.runid = runid_;
	hd.iz = iz;
	hd.is = -1; // has no meaning when all slices are sent at once
	hd.ns = ns;
	hd.harm = field->harm;
	hd.ngrid = field->ngrid;
	hd.dgrid = field->dgrid;
	hd.xlambda = field->xlambda;
	hd.do_multi = true;
	hd.datain = nullptr;
	hd.dataout = nullptr;
	hd.multi_datain = &multi_datain;
	hd.multi_dataout = &multi_dataout;

	/* execute the code in the library */
	li_.pdiagfield_->doit(&hd);
	if(interface_verbose_ && (my_rank_==0)) {
		if(multi_dataout.at(0).size()>0) {
			cout << "Rank 0: data from plugin=(";
			for(int is=0; is<ns; is++) {
				cout << "(";
				for(int kk=0; kk<multi_dataout.at(is).size(); kk++) {
					cout << multi_dataout.at(is).at(kk);
					if((1+kk)<multi_dataout.at(is).size())
						cout << ",";
				}
				cout << ")";
				if(is<(ns-1))
					cout << "; ";
			}
			cout << ")" << endl;
		} else {
			cout << "Rank 0: data from plugin=(plugin does not provide output data)" << endl;
		}
	}
	
	/* save the result into the provided array */
	for(int is=0; is<ns; is++) {
		int idx = iz*ns+is;         // compute index for saving the data
		for(int kk=0; kk<li_.obj_names_.size(); kk++) {
			stringstream ss;
			ss << obj_prefix_;
			ss << "/";
			ss << li_.obj_names_[kk];
			update_data(val, ss.str(), idx, multi_dataout.at(is).at(kk));
		}
	}
}

void DiagFieldHook::getValues(Field *field, std::map<std::string,std::vector<double> >&val, int iz)
{
	if(interface_verbose_ && (my_rank_==0)) {
		cout << "DiagFieldHook::getValues for iz=" << iz << ", h=" << field->harm << " (obj_prefix=" << obj_prefix_ << ")" << endl;
	}
	
	if(li_.supports_multimode()) {
		if(interface_verbose_ && (my_rank_==0)) {
			cout << "multimode enabled" << endl;
		}
		getValues_multiworker(field, val, iz);
	} else {
		getValues_worker(field, val, iz);
	}
}

