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
#include <algorithm>
#include <cstring>
#include <mpi.h>

#include "DiagnosticHook.h"
#include "DiagnosticHookS.h"
#include "GaussHermite.h"
#include "g4_tb_util.h"
#include "version.h"

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

void dump_result_mtx(ostream& os, const vector<double>& d, const int nz, const int nslice)
{
	// Fast index is slice id
	size_t idx=0;
	for(int iz=0; iz<nz; iz++) {
		os << "      [";
		for(int idslice=0; ; ) {
			os << d.at(idx);
			idx++;
			idslice++;
			if(idslice>=nslice) {
				break;
			}
			os << ",";
		}
		os << "]";
		if(iz<(nz-1))
			os << ",";
		os << endl;
	}
}
void dump_results_core(ostream &os, const map< string,vector<double> >& r, const int nz, const int nslice)
{
	os << "# create empty Python dictionary" << endl;
	os << "data = {};" << endl;
	for(auto const &[name,data]: r) {
		// formatting of data matrix optimized for numpy.array
		os << "# data in \"" << name << "\"" << endl;
		os << "data['" << name << "'] = np.array([" << endl;
		dump_result_mtx(os, data, nz, nslice);
		os << "])" << endl << endl;
	}
}
void dump_results(ostream &os, TB_Cfg *ptbcfg, const map< string,vector<double> >& r)
{
	dump_results_core(os, r, ptbcfg->nz, ptbcfg->nslice);
}

void res_collect(TB_Cfg *ptbcfg, map< string,vector<double> > & globaldata, const map< string,vector<double> >& localdata)
{
	// remark: globaldata is only used on rank 0
	int rank, mpisize;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

	/* 1) enforce that all processes on MPI communicator have same map size (=number of elements)  */
	int map_nele = localdata.size();
	MPI_Bcast(&map_nele, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(map_nele != localdata.size()) {
		abort();
	}


	/* 2) compare keystrings (we know already that all processes have an identical number of map elements -> prevents hang-ups during MPI collective ops) */
	vector<string> keys;
	for (auto const &[k,v]: localdata) {
		keys.push_back(k);
	}
	stable_sort(keys.begin(), keys.end());
	for(auto const &k: keys) {
		const char *loc_s = k.c_str();
		char *buf = nullptr;
		int buf_slen = k.length();

		MPI_Bcast(&buf_slen, 1, MPI_INT, 0, MPI_COMM_WORLD);
		buf = new char[buf_slen + 1 /* space for trailing '\0' */]();
		if (0==rank)
			strncpy(buf, loc_s, buf_slen+1); /* source process */
		MPI_Bcast(buf, buf_slen, MPI_CHAR, 0, MPI_COMM_WORLD);

#if 0
		// dbg: inject deviation into received key-string, will it be detected?
		if(1==rank)
			buf[0] = '_';
#endif

		if(strcmp(buf, loc_s)) {
			abort(); // local key string is different from root process with rank 0
		}

		delete [] buf;
	}

	const size_t local_cfgout_size = ptbcfg->nslice*ptbcfg->nz;
	vector<double> rbuf, rbuf_sorted;
	double *pdest = nullptr;
	if(rank==0) {
		rbuf.resize(mpisize*local_cfgout_size);
		rbuf_sorted.resize(mpisize*local_cfgout_size);
		pdest = rbuf.data();
	}
	for (auto const &[k,v]: localdata) {
		// remark: not using the vector with keys compiled for the tests above (combination of "const map<... , ...>" and the [] operation does not compile, probably because [] operation could also insert a new object into the map if the key does not exist)
		MPI_Gather(v.data(), local_cfgout_size, MPI_DOUBLE, pdest, local_cfgout_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/* rearrange collected data
		 * - Every process we has a the local buffer, a 1-D vector with size nslice*nz
		 * - After the collection all these are stored one after another in one big 1-D array on rank 0
		 * -> there are nz interleaved sets of data blocks (block length is nslice, the stride is length of local data buffer, i.e. nslice*nz)
		 */
		if(rank==0) {
			const int nslice = ptbcfg->nslice;
			size_t idx_dest=0;
			for(size_t iz=0; iz<ptbcfg->nz; iz++) {
				size_t idx_src0 = iz*nslice;
				for (size_t irank=0; irank<mpisize; irank++) {
					size_t idx_src = idx_src0 + irank*local_cfgout_size;
					for(int kk=0; kk<nslice; kk++)
						rbuf_sorted.at(idx_dest++) = rbuf.at(idx_src++);
				}
			}
			globaldata[k] = rbuf_sorted;
		}
	}
}

int main(int argc, char **argv)
{
	int rank, mpisize;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

	if(argc!=2) {
		if(rank==0) {
			cout << "This program requires one parameter: The name of the configuration file." << endl << "Exiting." << endl;
		}
		MPI_Finalize();
		exit(1);
	}
	string fncfg(argv[1]);

	if(0==rank) {
		VersionInfo vi;
		cout << "g4_tb, build info: " << vi.Build() << endl
		     << "name of configuration file: " << fncfg << endl
		     << "mpisize = " << mpisize << endl << endl;
	}

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
		cout << "opened config file '" << fncfg << "'." << endl;
	}
	bool parse_ok=ptbcfg->update_from_stream(ifs);
	ifs.close();
	if(rank==0) {
		if(!parse_ok)
			cout << "   processed config file (there were errors)." << endl << endl;
		else
			cout << "   processed config file." << endl << endl;
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
	}

#if 0
	// test of data organization
	for(int iz=0; iz<ptbcfg->nz; iz++) {
		for (int is=0; is<ptbcfg->nslice; is++) {
			unsigned long long x = is + rank*ptbcfg->nslice; // 'global' slice id
			x = 10*iz;
			results["plugin/abc"].at(iz*ptbcfg->nslice + is) = x;
		}
	}
#endif

	if(0==rank) {
		cout << endl
		     << "**********************************" << endl
		     << "Dumping generated diagnostics data" << endl;
		dump_results(cout, ptbcfg, results);
	}

	/***************************************************************/
	/*** collect data from all processes on the MPI communicator ***/
	/***************************************************************/
	map< string,vector<double> > glbl_results; // !only populated with collected data on rank 0!
	res_collect(ptbcfg, glbl_results, results);
	if(0==rank) {
		const char *fnout = "plugin_data.txt";
		ofstream ofs;
		ofs.open(fnout, ofstream::out);
		dump_results_core(ofs, glbl_results, ptbcfg->nz, mpisize*ptbcfg->nslice);
		ofs.close();

		cout << "Data collected from all processes written to file '" << fnout << "'" << endl;
	}
    
	/*** TODO: clean up all resources ***/

	for(int kk=0; kk<fields->size(); kk++) {
		delete fields->at(kk);
	}
	delete fields;
	delete pdfh; // unloads library
	delete ptbcfg;

	MPI_Finalize();
	return(0);
}
