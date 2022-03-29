#include <cstdio>
#include <iostream>
#include <fstream>
#include <mpi.h>

#include "SemaFile.h"

using namespace std;

SemaFile::SemaFile() {
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
}
SemaFile::~SemaFile() { }

void SemaFile::remove(string fn) {
	if (0==my_rank_) {
		// to be implemented: use std::remove to delete semaphore file
	}
}

void SemaFile::put(string fn) {
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


