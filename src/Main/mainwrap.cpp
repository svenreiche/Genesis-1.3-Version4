

#include "genesis.h"

// very basic wrapper for genesis. Most of the genesis stuff is moved into genmain.


int main (int argc, char *argv[]) {


        //-------------------------------------------------------
        // init MPI and get size etc.
        //

        MPI_Init(&argc, &argv); //MPI
      	string filename (argv[argc-1]);
	string latname;
	genmain(filename,latname,false,false,false);

        MPI_Finalize(); // node turned off

        return 0;

}
