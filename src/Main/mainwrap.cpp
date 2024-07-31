#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <string>

#include <getopt.h>
#include <mpi.h>
#include <hdf5.h>
#ifdef FFTW
 #include <fftw3.h> // for 'fftw_version'
#endif

#include "genesis.h"
#include "version.h"

using namespace std;



void G4_usage(void)
{
    cout << endl
         << "-------------------------------------------------------------------------------"<< endl;
    cout << "Usage: genesis4 [-o Output] [-l Lattice] [-b Beamline] [-s Seed] [-h] Input" << endl;
    cout << "        -o Output:  equivalent to 'rootname' in the setup namelist " << endl;
    cout << "        -l Lattice:  equivalent to 'lattice' in the setup namelist " << endl;
    cout << "        -b Beamline: equivalent to 'beamline' in the setup namelist " << endl;
    cout << "        -s Seed: equivalent to 'seed' in the setup namelist " << endl;
    cout << "        -h: prints help-info " << endl;
    cout << "  Any specified command line argument will overwrite" << endl;
    cout << "  the corresponding entry in the setup namelist" << endl;
    cout << "-------------------------------------------------------------------------------" << endl;
}

void G4_usage_and_stop(void)
{
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(0==rank)
        G4_usage();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); // node turned off

    // non-zero exit code signals error, for instance to SLURM
    exit(1);
}

void G4_report_lib_versions(void)
{
    int rank;

    // only rank 0 prints infos...
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(0!=rank)
        return;


    cout << endl; // some whitespace

    char buf[MPI_MAX_LIBRARY_VERSION_STRING];
    int bufused=-1;
    MPI_Get_library_version(buf, &bufused);
    cout << "MPI library version string: " << buf << endl;

    // HDF5 library
    // Remark: there are macros H5_VERS_MAJOR/MINOR/RELEASE that define
    // the version at compile time, but we report version at runtime.
    // Reporting both if there is a difference...
    unsigned hdf5_maj=-1, hdf5_min=-1, hdf5_relnum=-1;
    H5get_libversion(&hdf5_maj, &hdf5_min, &hdf5_relnum);
    cout << "HDF5 library reports version: "
         << hdf5_maj << "." << hdf5_min << "." << hdf5_relnum << endl;
    // hdf5_relnum=42; // to test the following
    if((H5_VERS_MAJOR!=hdf5_maj) || (H5_VERS_MINOR!=hdf5_min) || (H5_VERS_RELEASE!=hdf5_relnum)) {
        cout << "   Note: HDF5 library from header file (at build time): "
             << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << endl;
    }


    // FFTW3: is there a function that gives the version at runtime? I didn't find one ...
    // But there is char fftw_version[] (checked for FFTW v3.3.8)
#ifdef FFTW
    cout << "FFTW version string: " << fftw_version << endl;
#else
    cout << "FFTW version string: not compiled in" << endl;
#endif

    cout << endl;
}

int main (int argc, char *argv[])
{
    //-------------------------------------------------------
    // init MPI and get rank
    //
    int rank;
    MPI_Init(&argc, &argv); //MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node

    if (rank == 0) {
        VersionInfo vi;
        cout << "---------------------------------------------" << endl;
        cout << "GENESIS - Version " << vi.Major() << "." << vi.Minor() << "." << vi.Rev();
        if (vi.isBeta()) { cout << " (beta)"; }
        cout << " has started..." << endl;
        cout << "Compile info: " << vi.Build() << endl;
    }


	// parse the command line arguments
	struct option opts[] = {
		{"output",	required_argument,	NULL, 'o'},
		{"lattice", 	required_argument,	NULL, 'l'},
		{"beamline",	required_argument,	NULL, 'b'},
		{"seed",	required_argument,	NULL, 's'},
		{"semaphore-file", required_argument,	NULL, 'S'},
		{"help",	no_argument,		NULL, 'h'},
		{"dbg-versions",no_argument,		NULL, 'd'},
		{NULL, 0, NULL, 0}	/* marks end of list => do not remove */
	};

	// Inhibit error messages by getopt from all MPI processes except rank 0
	// Note that these error messages go to stderr, so depending on your simulation setup, they can end up in a different file
	if(rank!=0) {
		// Inhibits error message: "unrecognized option '--xyz'"
		// Also inhibits other messages, for instance msg informing about missing required argument ("option '--xyz' requires an argument")
		opterr = 0;
	}

	bool got_filename=false;
	string filename;
	map<string,string> arguments; /* keys of this map correspond to the parameter names in &setup namelist */

	int curr_opt;
	while((curr_opt=getopt_long(
		argc, argv,
		"o:l:b:s:h",	/* options that can be specified by single letter (see docs for effect of ':') */
		opts,
		NULL)) != -1)
	{
		// Keep in mind that the value of 'optarg' changes with every iteration of this loop
		// => generate copy if needed
		switch(curr_opt)
		{
			case 'o':
				arguments["rootname"] = optarg;
				break;

			case 'l':
				arguments["lattice"] = optarg;
				break;

			case 'b':
				arguments["beamline"] = optarg;
				break;

			case 's':
				/* FIXME: add check to see if user-provided string can be converted to int */
				arguments["seed"] = optarg;
				break;

			case 'S':
				arguments["semaphore_file_name"] = optarg;
				break;


			case 'h':
                                G4_usage_and_stop();
				break;

			case 'd': // debug functions (at start-up phase)
				G4_report_lib_versions();
				break;


			case '?': /* fall-through */
			default:
				if(rank==0) {
					cout << "*** error: unknown command line option specified ***" << endl;
				}
                                G4_usage_and_stop();
				break;
		}
	}

	/* NOTE: GNU getopt functions permute argv, unless requested to not do so (POSIXLY_CORRECT environment variable). 'optind' points to the first non-option argument. */
	if(optind<argc)
	{
		got_filename=true;
		filename = argv[optind];
		// cout << "input file name: " << filename << endl;
	}

#if 0
    /* debug code */
    if(rank==0) {
        cout << endl << "### reporting options ###" << endl;

        if(got_filename)
            cout << "filename: " << filename << endl;
        else
            cout << "no input file specified" << endl;

        map<string,string>::iterator it;
        for(it = arguments.begin(); it != arguments.end(); ++it) {
            cout << it->first << ": " << it->second << endl;
        }
        cout << "### DONE ###" << endl;
    }
#endif

    if(!got_filename) {
        if(rank==0)
            cout << "*** error: no input file provided -- execution of GENESIS will end ***" << endl;

        G4_usage_and_stop();
    }

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); // node turned off
    exit(0);
#endif

    // call the core routine for genesis (returns 0 when successful)
    int sim_core_result = genmain(filename,arguments,false);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); // node turned off

    return(sim_core_result);
}
