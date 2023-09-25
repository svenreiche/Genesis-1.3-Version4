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
#include <dlfcn.h>

#include "DiagnosticHook.h"
#include "DiagnosticPlugin.h"
#include "DiagnosticHookS.h"

// some info sources:
// https://stackoverflow.com/questions/496664/c-dynamic-shared-library-on-linux/497158#497158
// https://tldp.org/HOWTO/C++-dlopen/thesolution.html
// https://0x00sec.org/t/c-dynamic-loading-of-shared-objects-at-runtime/1498


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
	if(libok && (destroyer_!=nullptr)) {
		if(my_rank_==0) {
			cout << "calling destroyer function" << endl;
		}
		destroyer_(pdiag_);
		pdiag_ = nullptr;
		libok=false;
	}
}

bool DiagFieldHook::init(DiagFieldPluginCfg *pin)
{
	/* Copy all needed infos from the configuration data */
	libfile_           = pin->libfile;
	obj_prefix_        = pin->obj_prefix;
	parameter_         = pin->parameter;
	lib_verbose_       = pin->lib_verbose;
	interface_verbose_ = pin->interface_verbose;

	if(my_rank_==0) {
		cout << "DiagFieldHook::init" << endl;
	}

	bool my_status = get_shared_lib(&h_lib_);
	bool ok = get_shared_lib_diag(my_status, "Unable to load shared library");
	if(!ok) {
		exit(1);
		return(false);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank_==0) {
		cout << "Got the library" << endl;
	}

	my_status = get_shared_lib_objs(&h_lib_);
	ok = get_shared_lib_diag(my_status, "Error when obtaining the required resources from the shared library");
	if(!ok) {
		exit(1);
		return(false);
	}

	libok = true;
	return(true);
}

void DiagFieldHook::set_runid(int runid_in)
{
	runid_ = runid_in;
}


bool DiagFieldHook::get_shared_lib_diag(bool my_status, const char *msg)
{
	int is = my_status;
	int glbl_is=0;

	// all nodes have to agree on the next steps...
	MPI_Allreduce(&is, &glbl_is, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
	// cout << "after reduce operation: global status is " << glbl_is << endl;
	if(glbl_is==0) {
		// something went wrong
		// -> report some details (only on rank 0)
		vector<int> isbuf(comm_size_,0);
		stringstream ssfailed;
		char *errormsg;

		errormsg = dlerror(); // call 'dlerror' on ALL nodes to clear error flags

		MPI_Gather(&is, 1, MPI_INT, isbuf.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(my_rank_==0) {
			bool all_fail=true;
			for(int k=0; k<comm_size_; k++) {
				all_fail &= (isbuf.at(k)==0);
				if(isbuf.at(k)==0)
					ssfailed << k << " ";
			}

			cout << "DiagFieldHook: Diagnostic message \"" << msg << "\" ";
			if(all_fail) {
				cout << "on ALL ranks" << endl;
				cout << "Status on rank 0: dlerror=" << errormsg << endl;
			} else {
				cout << "on ranks " << ssfailed.str() << endl;
			}
		}
		return(false);
	}

	return(true);
}

// after calling this function, always call the _diag function (it calls 'dlerror' function)
bool DiagFieldHook::get_shared_lib(h_dynamic_lib *hout)
{
#if 0
	// test code: simulate issue with loading of shared library on one of the nodes
	if(my_rank_==1) {
		return(false);
	}
#endif

	// Remark on flags:
	// . RTLD_NOW is used because it is easier to debug any issues with missing symbols (they are directly triggered by this command)
	//   Possible reasons are missing C++ (base) classes, etc.
	h_dynamic_lib h = dlopen(libfile_.c_str(), RTLD_NOW);
	bool load_ok = (nullptr==h) ? false : true;

	// need to call dlerror() here?

	if(!load_ok) {
		// if loading the library failed, don't update the handle
		return(false);
	}

	if(my_rank_==0) {
		cout << "Rank 0: Loaded the library, handle is " << h << endl;
	}

	*hout = h;
	return(true);
}

void DiagFieldHook::report_infos(DiagFieldHookInfos *pi)
{
	cout << "Rank " << my_rank_ << ": Got class instance" << endl;
	cout << "   info_txt=\"" << pi->info_txt << "\"" << endl;
	cout << "   do_multi=" << pi->do_multi << endl;
	if(pi->obj_names == nullptr) {
		cout << "   Note: obj_names=NULL (could be a software issue)" << endl;
	} else {
		bool got_objs=false;
		const std::vector<const char *> *p = pi->obj_names;

		for(auto const& the_obj: *p) {
			cout << "   provides object named " << the_obj << endl;
			got_objs=true;
		}
		if(!got_objs) {
			cout << "   does not provide objects" << endl;
		}
	}
}
void DiagFieldHook::clone_obj_names(DiagFieldHookInfos *pi)
{
	if(pi->obj_names == nullptr)
		return;


	const std::vector<const char *> *p = pi->obj_names;

#if 1
	obj_names_.clear();
	for(auto const& the_obj: *p) {
		obj_names_.push_back(the_obj);
	}
#endif
}

bool DiagFieldHook::get_shared_lib_objs(h_dynamic_lib *h)
{
	void *p_f = dlsym(*h, "factory");
	void *p_d = dlsym(*h, "destroy");
	if(p_f==nullptr) {
		return(false);
	}
	if(my_rank_==0) {
		cout << "Rank 0: Factory location is " << p_f << endl;
		cout << "Rank 0: Destroyer location is " << p_d << endl;
	}
	factory_   = reinterpret_cast<fptr_f>(p_f);
	destroyer_ = reinterpret_cast<fptr_d>(p_d);


	/*
	 * Obtain instance of diagnostic class in shared library from the factory
	 * !!! ALWAYS use destructor function provided by library to free the object, do NOT just use 'delete' !!!
	 */
	if(my_rank_==0) {
		cout << "Rank 0: Calling factory" << endl;
	}
	pdiag_ = factory_();

	DiagFieldHookInfos infos;
	infos.version = DIAGFIELD_DATA_STRUCTVERSION; // version number that is incremented when structure layout changes
	infos.parameter = parameter_;
	infos.mpi_rank = my_rank_;
	infos.mpi_size = comm_size_;
	if(my_rank_==0) {
		cout << "Rank 0: Calling get_infos" << endl;
	}
	pdiag_->get_infos(&infos);
	multimode_ = infos.do_multi;
	clone_obj_names(&infos);
	if(my_rank_==0) {
		report_infos(&infos);
	}

	return(true);
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

	for(auto const& the_obj: obj_names_) {
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
		vector<double> dataout(obj_names_.size(),-1.);

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
		pdiag_->doit(&hd);
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
		for(int kk=0; kk<obj_names_.size(); kk++) {
			stringstream ss;
			ss << obj_prefix_;
			ss << "/";
			ss << obj_names_[kk];
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
		multi_dataout.at(k).resize(obj_names_.size(),-1.);


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
	pdiag_->doit(&hd);
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
		for(int kk=0; kk<obj_names_.size(); kk++) {
			stringstream ss;
			ss << obj_prefix_;
			ss << "/";
			ss << obj_names_[kk];
			update_data(val, ss.str(), idx, multi_dataout.at(is).at(kk));
		}
	}
}

void DiagFieldHook::getValues(Field *field, std::map<std::string,std::vector<double> >&val, int iz)
{
	if(interface_verbose_ && (my_rank_==0)) {
		cout << "DiagFieldHook::getValues for iz=" << iz << ", h=" << field->harm << " (obj_prefix=" << obj_prefix_ << ")" << endl;
	}
	
	if(multimode_) {
		if(interface_verbose_ && (my_rank_==0)) {
			cout << "multimode enabled" << endl;
		}
		getValues_multiworker(field, val, iz);
	} else {
		getValues_worker(field, val, iz);
	}
}

