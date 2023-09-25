#ifndef GENESIS_1_3_VERSION4__DIAGHOOKS_H
#define GENESIS_1_3_VERSION4__DIAGHOOKS_H

// Interface for field diagnostics plugin.

#include <complex>
#include <vector>

/*
 * version number included in data structured used for information exchange.
 * Increment it if there is a change to the definitions and then recompile GENESIS and the shared libraries.
 */
#define DIAGFIELD_DATA_STRUCTVERSION 3

/* header file with the data structures/classes */

/* class holding configuration parameters extracted from input file */
class DiagFieldPluginCfg {
public:
	std::string obj_prefix;
	std::string libfile;
	std::string parameter;
	bool lib_verbose;
	bool interface_verbose;
};

class DiagFieldHookInfos {
public:
	// version number: increment if there is change to this definition and recompile GENESIS and the shared libraries
	int version;

	const char *info_txt;
	const std::vector<const char *> *obj_names;
	std::string parameter;
	
	// If true, the plugin processes all slices at once:
	// useful for plugins that need all slices at the same time,
	// such as spectral analysis.
	bool do_multi {false};

	// not really configuration parameters, but handy for debug msgs
	int mpi_rank, mpi_size;
};

/* data structure for exchange of information with shared library */
class DiagFieldHookData {
public:
	// version number: increment if there is change to this definition and recompile GENESIS and the shared libraries
	int version;

	int mpi_rank, mpi_size;

	// int runid;
	int iz;
	int is, ns;
	int harm;
	int ngrid;
	double dgrid;
	double xlambda;

	bool do_multi {false};
	const std::vector<std::complex <double> > *datain {nullptr};
	std::vector<double> *dataout {nullptr};
	std::vector<const std::vector<std::complex <double> > *> *multi_datain {nullptr};
	std::vector<std::vector<double> > *multi_dataout {nullptr};

	/* add new fields here */
	int verbose;

	int runid; // to be moved up once all plugins check the version of this structure
};

#endif //  GENESIS_1_3_VERSION4__DIAGHOOKS_H
