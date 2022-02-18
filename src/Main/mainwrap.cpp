#include "genesis.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>

#include <getopt.h>

//#include "VersionInfo.h"
#include "version.h"
using namespace std;

// very basic wrapper for genesis. Most of the genesis stuff is moved into genmain.
#include "hdf5.h"

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
    exit(0);
}

int main (int argc, char *argv[])
{
    //-------------------------------------------------------
    // init MPI and get rank
    //
    int rank;
    MPI_Init(&argc, &argv); //MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node

//    string filename (argv[argc-1]);  // input file is always last element
//    map<string,string> arguments;
//    bool ok=true;

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
		{NULL, 0, NULL, 0}	/* marks end of list => do not remove */
	};

	// inhibit error messages by getopt from all MPI processes except rank 0
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

    // return 0;
    return(sim_core_result);

#if 0
	if (argc == 1){
        if (rank==0) {
            cout << "*** Error: No input file specified - Execution of Genesis will end" << endl;
        }
	} else {
	  for (int i = 1; i < argc; i++){
          if ((strcmp(argv[i],"-h") ==0) || (strcmp(argv[argc-1],"-h")==0)) {
                ok = false;
                if (rank == 0) {
                    cout << endl << "-------------------------------------------------------------------------------"<< endl;
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
                break;
          }
          // checking for the individual commandline argumen
          if (strcmp(argv[i],"-o")==0) {     // output file
              arguments["rootname"] = argv[i + 1];
              i++;
              continue;
          }
          if (strcmp(argv[i],"-l")==0) {     // output file
              arguments["lattice"] = argv[i + 1];
              i++;
              continue;
          }
          if (strcmp(argv[i],"-b")==0) {     // output file
              arguments["beamline"] = argv[i + 1];
              i++;
              continue;
          }
          if (strcmp(argv[i],"-s")==0) {     // output file
              arguments["seed"] = argv[i + 1];
              i++;
              continue;
          }
          // an unrecognized input is encountered
          if (i < (argc-1)) {
              if (rank == 0) {
                  cout << "*** Error: Invalid command line argument"<< endl;
                  cout << "*** Error: Argument " << argv[i] << " is not recognized" << endl;
                  cout << "*** Error: Execution of Genesis will end." << endl;
              }
              ok = false;
              break;
          }
	  }

      //   call the core routine for genesis
	  if (ok) {genmain(filename,arguments,false); }   //
	}

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); // node turned off

    return 0;
#endif
}
