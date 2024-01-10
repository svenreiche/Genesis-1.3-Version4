/*
 * Simple handshaking with files.
 * 
 * C. Lechner, EuXFEL, 2023-Nov
 * 
 * Was implemented to trigger processing of dumped field distribution
 * before subsequent loading. Of course, the processing program needs
 * to be already running...
 * 
 * FIXME: Currently uses MPI_Barrier, which does busy waiting.
 */

#include <fstream>
#include <iostream>
#include <unistd.h>
#include <mpi.h>

#include "SimpleHandshake.h"

using namespace std;

SimpleHandshake::SimpleHandshake() {
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);

	// these are the filenames relative to the simulation output directory
	// -> not the complete filenames
	fn_wait_ = "g4_hs.wait";
	fn_resume_ = "g4_hs.resume";
}

void SimpleHandshake::drop_file(const string fn)
{
	// code from SemaFile class
	if (0==my_rank_) {
		ofstream ofs;

		// Note: already existing files are truncated
		ofs.open(fn, ofstream::out);
		if(!ofs) {
			cout << "error generating semaphore file " << fn << endl;
		}
		ofs.close();
	}
}

void SimpleHandshake::wait_for_file(const string fname)
{
	const string msg_indent = "   ";

	if(0==my_rank_) {
		cout << msg_indent << "Waiting for file \"" << fname << "\" to be created."<< endl;
		while(1) {
			sleep(1);

			// try to remove the flag file (if this operation fails with ENOENT the file was there yet)
			int r = unlink(fname.c_str());
			if(0==r) {
				// Success: this means the file was there and could therefore be deleted
				cout << msg_indent << "Detected (and removed) file \"" << fname << "\"" << endl;
				break;
			}
		}
	}

	/*
	 * MPI_Barrier: All processes wait until rank 0 arrives as well.
	 * FIXME: Busy waiting on all CPUs, replace the MPI_Barrier by something better that puts processes to sleep.
	 */
	MPI_Barrier(MPI_COMM_WORLD);
	if(0==my_rank_) {
		cout << msg_indent << "Execution of GENESIS continues" << endl;
	}
}

bool SimpleHandshake::doit(const string prefix)
{
	if(0==my_rank_) {
		cout << endl << "SimpleHandshake" << endl;
	}

	// if output directory is defined, prefix has a trailing "/"
	// if no output dir defined, prefix is empty string
	const string path_wait   = prefix+fn_wait_;
	const string path_resume = prefix+fn_resume_;

	// Delete any "resume" file that may still be there from a previous
	// run before dropping the "wait" file (not handling the error
	// that we get if there is no such file)
	unlink(path_resume.c_str());
	drop_file(path_wait);
	wait_for_file(path_resume);

	if(0==my_rank_) {
		cout << endl;
	}

	return(true);
}
