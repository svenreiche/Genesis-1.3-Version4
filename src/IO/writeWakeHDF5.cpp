
#include "writeWakeHDF5.h"

extern bool MPISingle;

// constructor destructor
WriteWakeHDF5::WriteWakeHDF5() : fid(0) {}

WriteWakeHDF5::~WriteWakeHDF5() {}


bool WriteWakeHDF5::write(string fileroot, double *s, double *wakeres, double *wakegeo, double *wakerou, unsigned long ns)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node
  MPI_Comm_size(MPI_COMM_WORLD, &size); // assign rank to node
  if (MPISingle) {
    size = 1;
    rank = 0;
  }

  s0 = rank;  // selects the process which is writing the actual data
  string filename;
  filename = fileroot + ".wake.h5";
  if (rank == 0) { cout << "Writing single particle wakes to file: " << filename << " ..." << endl; }
  if (!create_outfile(&fid, filename)) {
    return (false);
  }

  this->writeSingleNodePointer(fid, "s", "m", s, ns);
  this->writeSingleNodePointer(fid, "resistive_wall", "m", wakeres, ns);
  this->writeSingleNodePointer(fid,"geometric","m",wakegeo,ns);
  this->writeSingleNodePointer(fid,"surface_roughness","m",wakerou,ns);
  H5Fclose(fid);
  return(true);
}

