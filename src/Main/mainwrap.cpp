#include "genesis.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

// very basic wrapper for genesis. Most of the genesis stuff is moved into genmain.


int main (int argc, char *argv[]) {


        //-------------------------------------------------------
        // init MPI and get size etc.
        //

        MPI_Init(&argc, &argv); //MPI

	int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node

      	string filename (argv[argc-1]);  // input file is always last element
	string latname  ("");
	string outname  ("");
	int    seed = 123456789;
	
	bool ok=true;

	// parse the command line arguments

	if (argc == 1){
	  if (rank == 0){
	    cout << "*** Error: No input file specified - Execution of Genesis will end" << endl;
	  }
	} else {
	  for (int i = 1; i < argc-1; i++){
	    if (strstr(argv[i],"-o")){     // output file
		if (i < argc-2){
		  outname=argv[i+1];
		  i++;
		  continue;
		}
	    }
	    if (strstr(argv[i],"-l")){
	      if (i < argc-2){         // lattice file
		  latname=argv[i+1];
		  i++;
		  continue;
		}
	    }

	    if (strstr(argv[i],"-s")){
	      if (i < argc-2){        // seed
		seed=atoi(argv[i+1]);
		  i++;
		  continue;
		}
	    }

	    if (rank == 0){
		cout << "*** Error: Invalid command line argument - Execution of Genesis will end" << endl;
	    }
	    ok=false;
	    break;
	  }       
	  if (ok) {genmain(filename,latname,outname,seed,false);}
	}
        MPI_Finalize(); // node turned off

        return 0;

}
