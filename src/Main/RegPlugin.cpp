#include <cmath>
#include "RegPlugin.h"
#include "DiagnosticHookS.h"

using namespace std;

void AddPluginCommon::usage(){
	string s = nameliststr();
	cout << "List of keywords for " << s << endl;
	cout << "&" << s << endl;
	cout << " string libfile = <empty>" << endl;
	cout << " string obj_prefix = plugin" << endl;
	cout << " bool verbose = false" << endl;
	cout << " bool interface_verbose = false" << endl;
	cout << "&end" << endl << endl;
	return;
}

bool AddPluginCommon::init(int rank, int size, map<string,string> *arg, Setup *setup)
{
	map<string,string>::iterator end=arg->end();
	string libfile = "";
	string obj_prefix = "plugin";
	string parameter = "";
	bool verbose = false;
	bool interface_verbose = false;

	// extract parameters
	if (arg->find("libfile")!=end)      {libfile = arg->at("libfile");        arg->erase(arg->find("libfile"));}
	if (arg->find("obj_prefix")!=end)   {obj_prefix = arg->at("obj_prefix");  arg->erase(arg->find("obj_prefix"));}
	if (arg->find("parameter")!=end)    {parameter = arg->at("parameter");    arg->erase(arg->find("parameter"));}
	if (arg->find("verbose")!=end)      {verbose = atob(arg->at("verbose"));  arg->erase(arg->find("verbose"));}
	if (arg->find("interface_verbose")!=end)  {interface_verbose = atob(arg->at("interface_verbose"));  arg->erase(arg->find("interface_verbose"));}

	if (arg->size()!=0){
		if (rank==0) {
			cout << "*** Error: Unknown elements in &" << nameliststr() << endl;
			usage();
		}
		return(false);
	}


	// check parameters
	if(libfile.size()<1) {
		if (rank==0) {
			cout << "*** Error in " << nameliststr() << ": name of shared library to load not specified" << endl;
		}
		return(false);
	}

#if 0
	if (rank==0) {
		cout << "parameter=" << parameter << endl;
	}
#endif

	/* register diagnostic plugin for activation */
	if(id() == PLUGIN_FIELD) {
		DiagFieldPluginCfg cfg;
		cfg.obj_prefix = obj_prefix;
		cfg.libfile = libfile;
		cfg.parameter = parameter;
		cfg.lib_verbose = verbose; // controls verbosity in the loaded .so lib
		cfg.interface_verbose = interface_verbose;
		setup->diagpluginfield_.push_back(cfg);
	} else {
		DiagBeamPluginCfg cfg;
		cfg.obj_prefix = obj_prefix;
		cfg.libfile = libfile;
		cfg.parameter = parameter;
		cfg.lib_verbose = verbose; // controls verbosity in the loaded .so lib
		cfg.interface_verbose = interface_verbose;
		setup->diagpluginbeam_.push_back(cfg);
	}

	return(true);
}

PluginType AddPluginFieldDiag::id(void)
{
	return(PLUGIN_FIELD);
}
string AddPluginFieldDiag::nameliststr(void)
{
	return("add_plugin_fielddiag");
}
PluginType AddPluginBeamDiag::id(void)
{
	return(PLUGIN_BEAM);
}
string AddPluginBeamDiag::nameliststr(void)
{
	return("add_plugin_beamdiag");
}
