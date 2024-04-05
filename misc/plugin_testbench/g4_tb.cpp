/*
 * Testbench for plugins
 *
 * Christoph Lechner, European XFEL, 2024-Apr
 */
#include <fstream>
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

void prepare_field(TB_Cfg *ptbcfg, vector<Field *> *fieldin)
{
	/* based on LoadField::init in LoadField.cpp */
	Field *field = new Field;
	const int ngrid = ptbcfg->ngrid;
	field->init(ptbcfg->nslice,ngrid,ptbcfg->dgrid,ptbcfg->lambda,ptbcfg->sample*ptbcfg->lambda,ptbcfg->s0,ptbcfg->harm);
	fieldin->push_back(field);
	int idx=fieldin->size()-1;

	complex< double >  *fieldslice = new complex<double> [ngrid*ngrid];
	FieldSlice slice;
	GaussHermite gh;

	for (int j=0; j<ptbcfg->nslice; j++)
	{
		// int i=j+time->getNodeOffset(); // FIXME
		int i=j; // +time->getNodeOffset(); // FIXME
		slice.lambda = ptbcfg->lambda;
		slice.power =  ptbcfg->power;
		slice.phase = 0;
		slice.z0 = 0;
		slice.w0 = ptbcfg->w0;
		slice.xcen = 0;
		slice.ycen = 0;
		slice.xangle = 0;
		slice.yangle = 0;
		slice.nx = 0;
		slice.ny = 0;
		slice.harm = ptbcfg->harm;
		gh.loadGauss(fieldslice,&slice,ptbcfg->dgrid,ngrid);
		for (int k=0; k<ngrid*ngrid;k++) {
			fieldin->at(idx)->field[j].at(k)=fieldslice[k];
		}
	}
	delete[] fieldslice;
}

void dump_result_mtx(TB_Cfg *ptbcfg, const vector<double>& d)
{
	// Fast index is slice id
	size_t idx=0;
	for(int iz=0; iz<ptbcfg->nz; iz++) {
		cout << "      [";
		for(int idslice=0; ; ) {
			cout << d.at(idx);
			idx++;
			idslice++;
			if(idslice>=ptbcfg->nslice) {
				break;
			}
			cout << ",";
		}
		cout << "]";
		if(iz<(ptbcfg->nz-1))
			cout << ",";
		cout << endl;
	}
}
void dump_results(TB_Cfg *ptbcfg, const map< string,vector<double> >& r)
{
	for(auto const &obj: r) {
		cout << "   data in \"" << obj.first << "\"" << endl;
		dump_result_mtx(ptbcfg, obj.second);
	}
}

int main(int argc, char **argv)
{
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(argc!=2) {
		if(rank==0) {
			cout << "This program requires one parameter: The name of the configuration file." << endl << "Exiting." << endl;
		}
		MPI_Finalize();
		exit(1);
	}
	string fncfg(argv[1]);

	/* prepare the configuration parameters */
	TB_Cfg *ptbcfg = new TB_Cfg;

	ifstream ifs;
	ifs.open(fncfg, ifstream::in);
	if(!ifs.good()) {
		if(rank==0) {
			cout << "Issue reading configuration file '" << fncfg << "'." << endl << "Exiting." << endl;
		}
		MPI_Finalize();
		exit(1);
	}
	if(rank==0) {
		cout << "opened config file '" << fncfg << "." << endl;
	}
	bool parse_ok=ptbcfg->update_from_stream(ifs);
	ifs.close();
	if(rank==0) {
		if(!parse_ok)
			cout << "processed config file (there were errors)." << endl << endl;
		else
			cout << "processed config file." << endl << endl;
	}


	/***********************************************************/
	/*** set up library interface and attempt to load plugin ***/
	/***********************************************************/

	/* based on RegPlugin.cpp (commit id 4681421, date: 2024-03-20) */
	DiagFieldPluginCfg cfg;
	cfg.libfile = ptbcfg->libfile;
	cfg.parameter = ptbcfg->parameter;
	cfg.obj_prefix = "plugin";
	cfg.lib_verbose = false; // controls verbosity in the loaded .so lib
	cfg.interface_verbose = false;

   
	/* copied from Gencore.cpp (commit id 4681421, date: 2024-03-20) */
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

	/**********************/
	/*** prepare fields ***/
	/**********************/
	const size_t cfgout_size = ptbcfg->nslice*ptbcfg->nz;
	vector<Field *> *fields = new vector<Field *>;
	prepare_field(ptbcfg, fields);
	// Remark, CL, 2024-03-28: Not using field import code from src/Loading/ImportField.cpp yet, because one would have to prepare an instance of class 'Setup'.
	/* setup data structures for diagnostics data, code from Diagnostic::init in Diagnostic.cpp */
    
	// Remark on data organization:
	// in Diagnostic.h:
	// vector< map< string,vector<double> > > val;
	// here we use (just keeping data for a single harmonic field):
	map< string,vector<double> > results;

	FilterDiagnostics filter;
	for(auto const &tag: pdfh->getTags(filter)) {
		results[tag.first].resize(cfgout_size);
	}

	/**********************************************/
	/*** run the diagnostics code in the plugin ***/
	/**********************************************/
	for(int iz=0; iz<ptbcfg->nz; iz++)
	{
		for(int ifld=0; ifld<fields->size(); ifld++) {
			/* this calls the code provided by the plugin, so put breakpoints there if needed */
			pdfh->getValues(fields->at(ifld), results, iz);
		}

#if 0
		// scale the field for test of data organization
		for (int j=0; j<ptbcfg->nslice; j++) {
			for (int k=0; k<ptbcfg->ngrid*ptbcfg->ngrid;k++){
				const int idx=0;
				fields->at(idx)->field[j].at(k)*=2.0; 
			}
		}
#endif
	}

	if(0==rank) {
		cout << endl
		     << "**********************************" << endl
		     << "Dumping generated diagnostics data" << endl;
		dump_results(ptbcfg, results);
	}
    
	/*** TODO: clean up all resources, unload the plugin module etc. ***/

	for(int kk=0; kk<fields->size(); kk++) {
		delete fields->at(kk);
	}
	delete fields;
	delete ptbcfg;

	MPI_Finalize();
	return(0);
}
