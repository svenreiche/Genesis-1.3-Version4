/*
 * Simple handshaking with files.
 * 
 * C. Lechner, EuXFEL, 2023-Nov/2024-Jan
 * 
 * Was implemented to trigger processing of dumped field distribution
 * before subsequent loading. Of course, the processing program needs
 * to be already running...
 */

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unistd.h>
#include <mpi.h>

#include "SimpleHandshake.h"

using namespace std;

SimpleHandshake::SimpleHandshake() {
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
	
	// controls MPI synchronization method
	busy_wait_ = false;

	// these are the file extensions
	ext_wait_ = ".wait";
	ext_resume_ = ".resume";
}

void SimpleHandshake::usage(void)
{
  cout << "List of keywords for simple_handshake" << endl;
  cout << "&simple_handshake" << endl;
  cout << " string file = <g4_hs>" << endl;
  cout << "&end" << endl << endl;
}

void SimpleHandshake::join(void)
{
	if(busy_wait_) {
		/*
		 * MPI_Barrier: All processes wait until rank 0 arrives as well.
		 * Note: With OpenMPI, this results in busy waiting on all CPUs,
		 */
		MPI_Barrier(MPI_COMM_WORLD);
		return;
	}
	
	/*
	 * Improved sync that puts waiting processes to sleep
	 * (so that CPUs can be be used for other tasks)
	 */
	MPI_Request req;
	MPI_Ibarrier(MPI_COMM_WORLD, &req);
	
	while(1) {
		int flag=-1;
		MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
		if(flag) {
			break;
		}
		
		sleep(1);
		// cout << "waiting on rank " << my_rank_ << endl;
	}
	return;
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

	/* Here, all processes wait until rank 0 arrives as well. */
	join();
	if(0==my_rank_) {
		cout << msg_indent << "Execution of GENESIS continues" << endl;
	}
}

bool SimpleHandshake::doit(map<string,string> *arg, const string prefix)
{
	if(0==my_rank_) {
		cout << endl << "SimpleHandshake" << endl;
	}

	auto end=arg->end();
	string file = "g4_hs";
	if (arg->find("file")!=end){file = arg->at("file"); arg->erase(arg->find("file"));}
	if (!arg->empty()){
	  if (my_rank_==0){ cout << "*** Error: Unknown elements in &simple_handshake" << endl; usage();}
	  return(false);
	}

	// If output directory is defined, prefix has a trailing "/",
	// if no output dir defined, prefix is empty string.
	// All names relative to the current working directory
	const string path_wait   = prefix+file+ext_wait_;
	const string path_resume = prefix+file+ext_resume_;

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
