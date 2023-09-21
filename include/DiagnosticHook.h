#ifndef GENESIS_1_3_VERSION4__DIAGHOOK_H
#define GENESIS_1_3_VERSION4__DIAGHOOK_H

#include <map>
#include <string>
#include <vector>

#include "DiagnosticBase.h"

class DiagFieldHookedBase;
class DiagFieldHookInfos;
class DiagFieldPluginCfg;

// data type for factory and destructor (not C++ destrcutor)
typedef DiagFieldHookedBase* (*fptr_f)();
typedef void (*fptr_d)(DiagFieldHookedBase *);

typedef void* h_dynamic_lib;

class DiagFieldHook: public DiagFieldBase {
public:
	DiagFieldHook();
	~DiagFieldHook();

	std::map<std::string,OutputInfo> getTags(FilterDiagnostics &);
	void getValues(Field *, std::map<std::string,std::vector<double> > &, int);

	bool init(DiagFieldPluginCfg *);
	void set_runid(int);

private:
	bool get_shared_lib_diag(bool, const char *);
	bool get_shared_lib(h_dynamic_lib *);
	bool get_shared_lib_objs(h_dynamic_lib *);
	void report_infos(DiagFieldHookInfos *);
	void clone_obj_names(DiagFieldHookInfos *);
	bool update_data(std::map<std::string,std::vector<double> > &, std::string, size_t, double);
	
	void getValues_worker(Field *, std::map<std::string,std::vector<double> >&, int);
	void getValues_multiworker(Field *, std::map<std::string,std::vector<double> >&, int);

	int my_rank_, comm_size_;

	int runid_ {1};

	bool libok {false};
	std::string libfile_;
	h_dynamic_lib h_lib_;
	fptr_f factory_ {nullptr};
	fptr_d destroyer_ {nullptr};
	DiagFieldHookedBase *pdiag_ {nullptr};
	std::vector<std::string> obj_names_;
	bool multimode_ {false};


	std::string obj_prefix_;
	std::string parameter_;
	bool lib_verbose_ {false};
	bool interface_verbose_ {false};
};


#endif
