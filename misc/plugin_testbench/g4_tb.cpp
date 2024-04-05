#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <mpi.h>

#include "DiagnosticHook.h"
#include "DiagnosticHookS.h"
#include "GaussHermite.h"
#include "g4_tb_util.h"

using namespace std;

void prepare_field(TB_Cfg *cfg, vector<Field *> *fieldin)
{
	/*** based on LoadField::init in LoadField.cpp ***/
	Field *field = new Field;
	const int ngrid = cfg->ngrid;
	field->init(cfg->nslice,ngrid,cfg->dgrid,cfg->lambda,cfg->sample*cfg->lambda,cfg->s0,cfg->harm);
	fieldin->push_back(field);
	int idx=fieldin->size()-1;

	complex< double >  *fieldslice = new complex<double> [ngrid*ngrid];
	FieldSlice slice;
	GaussHermite gh;

	for (int j=0; j<cfg->nslice; j++)
	{
		// int i=j+time->getNodeOffset(); // FIXME
		int i=j; // +time->getNodeOffset(); // FIXME
		slice.lambda = cfg->lambda;
		slice.power = 1e6;
		slice.phase = 0;
		slice.z0 = 0;
		slice.w0 = 20e-6;
		slice.xcen = 0;
		slice.ycen = 0;
		slice.xangle = 0;
		slice.yangle = 0;
		slice.nx = 0;
		slice.ny = 0;
		slice.harm = cfg->harm;
		gh.loadGauss(fieldslice,&slice,cfg->dgrid,ngrid);
		for (int k=0; k<ngrid*ngrid;k++) {
			fieldin->at(idx)->field[j].at(k)=fieldslice[k];
		}
	}
	delete[] fieldslice;
}

void dump_result_mtx(TB_Cfg *cfg, const vector<double>& d)
{
	// Fast index is slice id
	size_t idx=0;
	for(int iz=0; iz<cfg->nz; iz++) {
		cout << "      [";
		for(int idslice=0; ; ) {
			cout << d.at(idx);
			idx++;
			idslice++;
			if(idslice>=cfg->nslice) {
				break;
			}
			cout << ",";
		}
		cout << "]";
		if(iz<(cfg->nz-1))
			cout << ",";
		cout << endl;
	}
}
void dump_results(TB_Cfg *pcfg, const map< string,vector<double> >& r)
{
	for(auto const &obj: r) {
		cout << "   data in \"" << obj.first << "\"" << endl;
		dump_result_mtx(pcfg, obj.second);
	}
}

int main(int argc, char **argv)
{
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/*** copied from RegPlugin.cpp (commit id 4681421, date: 2024-03-20) ***/
	string libfile = "./libdemo.so";
	string obj_prefix = "plugin";
	string parameter = "";
	bool verbose = false;
	bool interface_verbose = false;

	DiagFieldPluginCfg cfg;
	cfg.obj_prefix = obj_prefix;
	cfg.libfile = libfile;
	cfg.parameter = parameter;
	cfg.lib_verbose = verbose; // controls verbosity in the loaded .so lib
	cfg.interface_verbose = interface_verbose;

   
	/*** copied from Gencore.cpp (commit id 4681421, date: 2024-03-20) ***/
	if(rank==0) {
		cout << "Setting up DiagFieldHook for libfile=\"" << cfg.libfile << "\", obj_prefix=\"" << cfg.obj_prefix << "\"" << endl;
	}
	DiagFieldHook *pdfh = new DiagFieldHook(); /* !do not delete this instance, it will be destroyed when DiagFieldHook instance is deleted! */
	bool diaghook_ok = pdfh->init(&cfg);
	if(diaghook_ok)
	{
		// pdfh->set_runid(setup->getCount()); // propagate run id so that it can be used in the plugins, for instance for filename generation
		pdfh->set_runid(1); // propagate run id so that it can be used in the plugins, for instance for filename generation

		string tmp_infotxt = pdfh->get_info_txt();
		// und->plugin_info_txt.push_back(tmp_infotxt);
		stringstream tmp_prefix;
		tmp_prefix << "/Field/" << cfg.obj_prefix;
		// und->plugin_hdf5_prefix.push_back(tmp_prefix.str());

		// diag.add_field_diag(pdfh);
		if(rank==0) {
			cout << "DONE: Registered DiagFieldHook" << endl;
		}
	} else {
		delete pdfh;
		if(rank==0) {
			cout << "failed to set up DiagFieldHook, exiting" << endl;
		}
		MPI_Finalize();
		exit(1);
	}

	TB_Cfg *pcfg = new TB_Cfg;
	vector<Field *> *fields = new vector<Field *>;
	prepare_field(pcfg, fields);
	/*** Remark, CL, 2024-03-28: Not using field import code from src/Loading/ImportField.cpp yet, because one would have to prepare an instance of class 'Setup'. ***/

	/*** setup data structures for diagnostics data, code from Diagnostic::init in Diagnostic.cpp ***/
    
	// Remark on data organization:
	// in Diagnostic.h:
	// vector< map< string,vector<double> > > val;
	// here we use (just keeping data for a single harmonic field):
	map< string,vector<double> > results;

	FilterDiagnostics filter;
	for(auto const &tag: pdfh->getTags(filter)) {
		int size = pcfg->nslice*pcfg->nz;
		results[tag.first].resize(size);
	}

	/*** NOW, run the diagnostics code ***/
	for(int iz=0; iz<pcfg->nz; iz++)
	{
		for(int ifld=0; ifld<fields->size(); ifld++) {
			/* this calls the code provided by the plugin, so put breakpoints there if needed */
			pdfh->getValues(fields->at(ifld), results, iz);
		}

		// scale the field for the test
		for (int j=0; j<pcfg->nslice; j++) {
			for (int k=0; k<pcfg->ngrid*pcfg->ngrid;k++){
				const int idx=0;
				fields->at(idx)->field[j].at(k)*=2.0; 
			}
		}
	}

	cout << "**********************************" << endl;
	cout << "Dumping generated diagnostics data" << endl;
	dump_results(pcfg, results);
    
	/*** TODO: clean up all resources, unload the plugin module etc. ***/

	for(int kk=0; kk<fields->size(); kk++) {
		delete fields->at(kk);
	}
	delete fields;
	delete pcfg;

	MPI_Finalize();
	return(0);
}
