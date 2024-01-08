#ifndef GENESIS_1_3_VERSION4__DIAGHOOK_H
#define GENESIS_1_3_VERSION4__DIAGHOOK_H

#include <map>
#include <string>
#include <vector>

#include "DiagnosticBase.h"

class DiagCommonHookInfos;
class DiagBeamHookedBase;
class DiagBeamHookInfos;
class DiagBeamPluginCfg;
class DiagFieldHookedBase;
class DiagFieldHookInfos;
class DiagFieldPluginCfg;

// data type for factory and destructor (not C++ destructor)
typedef DiagFieldHookedBase* (*fptr_ff)();
typedef void (*fptr_df)(DiagFieldHookedBase *);
typedef DiagBeamHookedBase* (*fptr_fb)();
typedef void (*fptr_db)(DiagBeamHookedBase *);

typedef void* h_dynamic_lib;

class LibraryInterface {
public:
	LibraryInterface(void);
	
	bool init_lib(std::string);
	void close_lib(void);
	bool is_libok(void);
	bool is_plugintype_field(void);
	bool is_plugintype_beam(void);
	bool supports_multimode(void);
	const std::string& get_info_txt() const;

	DiagFieldHookedBase *pdiagfield_ {nullptr};
	DiagBeamHookedBase  *pdiagbeam_ {nullptr};
	std::vector<std::string> obj_names_;

	std::string parameter_;

private:
	bool get_shared_lib_diag(bool, const char *);
	bool get_shared_lib(h_dynamic_lib *);
	bool get_shared_lib_objs(h_dynamic_lib *);

	void report_common_infos(DiagCommonHookInfos *);
	void report_infos(DiagFieldHookInfos *);
	void report_infos(DiagBeamHookInfos *);
	void clone_obj_names(const std::vector<const char *> *);

	
	bool libok_ {false};
	std::string libfile_;
	h_dynamic_lib h_lib_;
	fptr_ff factory_ {nullptr};
	fptr_df destroyer_ {nullptr};
	fptr_fb beamdiag_factory_ {nullptr};
	fptr_db beamdiag_destroyer_ {nullptr};
	bool is_field_ {true};
	bool multimode_ {false};
	std::string info_txt_;
	
	int my_rank_, comm_size_;
};



class DiagBeamHook: public DiagBeamBase {
public:
	DiagBeamHook();
	~DiagBeamHook();

	std::map<std::string,OutputInfo> getTags(FilterDiagnostics &);
	void getValues(Beam *, std::map<std::string,std::vector<double> >&, int);

	bool init(DiagBeamPluginCfg *);
	void set_runid(int);
	const std::string& get_info_txt() const;

private:
	bool update_data(std::map<std::string,std::vector<double> > &, std::string, size_t, double);
	
	void getValues_worker(Beam *, std::map<std::string,std::vector<double> >&, int);
	void getValues_multiworker(Beam *, std::map<std::string,std::vector<double> >&, int);

	int my_rank_, comm_size_;

	int runid_ {1};

	std::string obj_prefix_;
	bool lib_verbose_ {false};
	bool interface_verbose_ {false};
	
	LibraryInterface li_;
};

class DiagFieldHook: public DiagFieldBase {
public:
	DiagFieldHook();
	~DiagFieldHook();

	std::map<std::string,OutputInfo> getTags(FilterDiagnostics &);
	void getValues(Field *, std::map<std::string,std::vector<double> > &, int);

	bool init(DiagFieldPluginCfg *);
	void set_runid(int);
	const std::string& get_info_txt() const;

private:
	bool update_data(std::map<std::string,std::vector<double> > &, std::string, size_t, double);
	
	void getValues_worker(Field *, std::map<std::string,std::vector<double> >&, int);
	void getValues_multiworker(Field *, std::map<std::string,std::vector<double> >&, int);

	int my_rank_, comm_size_;

	int runid_ {1};

	std::string obj_prefix_;
	bool lib_verbose_ {false};
	bool interface_verbose_ {false};
	
	LibraryInterface li_;
};


#endif
