#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <mpi.h>
#include <dlfcn.h>

#include "DiagnosticHook.h"
#include "DiagnosticPlugin.h"
#include "DiagnosticHookS.h"


// sources:
// https://stackoverflow.com/questions/496664/c-dynamic-shared-library-on-linux/497158#497158
// https://tldp.org/HOWTO/C++-dlopen/thesolution.html
// https://0x00sec.org/t/c-dynamic-loading-of-shared-objects-at-runtime/1498


LibraryInterface::LibraryInterface()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size_);
}

bool LibraryInterface::init_lib(string libfile_in)
{
	libfile_ = libfile_in;

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

	libok_ = true;
	return(true);
}

bool LibraryInterface::get_shared_lib_diag(bool my_status, const char *msg)
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
bool LibraryInterface::get_shared_lib(h_dynamic_lib *hout)
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

void LibraryInterface::report_infos(DiagFieldHookInfos *pi)
{
	cout << "Rank " << my_rank_ << ": Got class instance" << endl;
	cout << "   info_txt=\"" << pi->info_txt << "\"" << endl;
	cout << "   do_multi=" << pi->do_multi << endl;
	if(pi->obj_names == nullptr) {
		cout << "   Note: obj_names=NULL (could be a software issue)" << endl;
	} else {
		bool has_objs=false;
		const std::vector<const char *> *p = pi->obj_names;

		for(auto const& the_obj: *p) {
			cout << "   provides object named " << the_obj << endl;
			has_objs=true;
		}
		if(!has_objs) {
			cout << "   does not provide objects" << endl;
		}
	}
}

//void LibraryInterface::clone_obj_names(DiagFieldHookInfos *pi)
void LibraryInterface::clone_obj_names(const vector<const char *> *p)
{
#if 0
	if(pi->obj_names == nullptr)
		return;


	const std::vector<const char *> *p = pi->obj_names;
#endif
	if(p==nullptr)
		return;

	obj_names_.clear();
	for(auto const& the_obj: *p) {
		obj_names_.push_back(the_obj);
	}
}

bool LibraryInterface::get_shared_lib_objs(h_dynamic_lib *h)
{
	void *p_f = dlsym(*h, "factory");
	void *p_d = dlsym(*h, "destroy");
	void *p_fbeam = dlsym(*h, "factory_beamdiag");
	void *p_dbeam = dlsym(*h, "destroy_beamdiag");

	if(my_rank_==0) {
		cout << "Rank 0: FieldDiag Factory location is " << p_f << endl;
		cout << "Rank 0: FieldDiag Destroyer location is " << p_d << endl;
		cout << "Rank 0: BeamDiag Factory location is " << p_fbeam << endl;
		cout << "Rank 0: BeamDiag Destroyer location is " << p_dbeam << endl;
	}

	/* The found factory indicates if this is a field or a beam
	 * diag plugin. Ensure that exactly one is present.
	 */
	if((p_f==nullptr) && (p_fbeam==nullptr)) {
		return(false);
	}
	if((p_f!=nullptr) && (p_fbeam!=nullptr)) {
		return(false);
	}

	is_field_ = false;
	if(p_f!=nullptr) {
		is_field_ = true;
	}

	factory_            = reinterpret_cast<fptr_ff>(p_f);
	destroyer_          = reinterpret_cast<fptr_df>(p_d);
	beamdiag_factory_   = reinterpret_cast<fptr_fb>(p_fbeam);
	beamdiag_destroyer_ = reinterpret_cast<fptr_db>(p_dbeam);

	/*
	 * Obtain instance of diagnostic class in shared library from the factory
	 * !!! ALWAYS use destructor function provided by library to free the object, do NOT just use 'delete' !!!
	 */
	if(my_rank_==0) {
		cout << "Rank 0: Calling factory" << endl;
	}
	if(is_field_) {
		pdiagfield_ = factory_();
		pdiagbeam_ = nullptr;
	} else {
		pdiagbeam_ = beamdiag_factory_();
		pdiagfield_ = nullptr;
	}


	if(is_field_) {
		DiagFieldHookInfos infos;
		infos.version = DIAGFIELD_DATA_STRUCTVERSION; // version number that is incremented when structure layout changes
		infos.parameter = parameter_;
		infos.mpi_rank = my_rank_;
		infos.mpi_size = comm_size_;
		if(my_rank_==0) {
			cout << "Rank 0: Calling get_infos" << endl;
		}
		pdiagfield_->get_infos(&infos);
		multimode_ = infos.do_multi;
		clone_obj_names(infos.obj_names);
		if(my_rank_==0) {
			report_infos(&infos);
		}
	} else {
		DiagBeamHookInfos infos;
		infos.version = DIAGFIELD_DATA_STRUCTVERSION; // version number that is incremented when structure layout changes
		infos.parameter = parameter_;
		infos.mpi_rank = my_rank_;
		infos.mpi_size = comm_size_;
		if(my_rank_==0) {
			cout << "Rank 0: Calling get_infos" << endl;
		}
		pdiagbeam_->get_infos(&infos);
		multimode_ = infos.do_multi;
		clone_obj_names(infos.obj_names);
		if(my_rank_==0) {
			// report_infos(&infos);
			cout << "Info: Reporting for beam plugins to be implemented" << endl;
		}
	}

	return(true);
}

void LibraryInterface::close_lib(void)
{
	if(!libok_)
		return;

	if(destroyer_!=nullptr) {
		if(my_rank_==0) {
			cout << "calling destroyer function" << endl;
		}
		destroyer_(pdiagfield_);
	}
	if(beamdiag_destroyer_!=nullptr) {
		if(my_rank_==0) {
			cout << "calling destroyer function" << endl;
		}
		beamdiag_destroyer_(pdiagbeam_);
	}
	pdiagfield_ = nullptr;
	pdiagbeam_ = nullptr;
	libok_=false;
}

bool LibraryInterface::supports_multimode(void)
{
	return(multimode_);
}

bool LibraryInterface::is_libok(void)
{
	return(libok_);
}
