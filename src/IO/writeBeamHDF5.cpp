
#include "writeBeamHDF5.h"
#include "MPEProfiling.h"

// constructor destructor
WriteBeamHDF5::WriteBeamHDF5()
{
}

WriteBeamHDF5::~WriteBeamHDF5()
{
}


void WriteBeamHDF5::write(string fileroot, Beam *beam){

  MPI::Status status;

  int size=MPI::COMM_WORLD.Get_size(); // get size of cluster
  int rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node

  int tag=1;
  int slicecount=0;
  int ntotal=0;
  int nslice=beam->beam.size();

  MPI::COMM_WORLD.Reduce(&nslice,&ntotal,1,MPI::INT,MPI::SUM,0);

  mpe.logIO(false,true,"Write Beam Dump to File");  


  char filename[100];
  sprintf(filename,"%s.par.h5",fileroot.c_str());

  if (rank==0){

    hid_t fid=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT); 
    this->writeGlobal(fid,beam->nbins,beam->one4one,beam->reflength,beam->slicelength,beam->s0,ntotal);
    slicecount+=this->writeSlice(fid, beam, slicecount);
    H5Fclose(fid);

    if (size > 1){
      	MPI::COMM_WORLD.Send(&slicecount, 1, MPI::INT, rank+1, tag);
    }

  } else {

    MPI::COMM_WORLD.Recv(&slicecount, 1, MPI::INT, rank-1, tag,status);

    hid_t fid=H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);   
    slicecount+=this->writeSlice(fid, beam, slicecount);
    H5Fclose(fid);

    if (rank < (size-1)){
      	MPI::COMM_WORLD.Send(&slicecount, 1, MPI::INT, rank+1, tag);
    }
  }

  mpe.logIO(true,true,"Write Beam Dump to File");  

  return;
}

void WriteBeamHDF5::writeGlobal(hid_t fid,int nbins,bool one4one, double reflen, double slicelen, double s0, int count)
{

  int tmp=0;
  if (one4one) { tmp=1; }

  this->writeDataDouble(fid, "slicelength", &reflen, 1);
  this->writeDataDouble(fid, "slicespacing", &slicelen, 1);
  this->writeDataDouble(fid, "refposition", &s0, 1);

  this->writeDataInt(fid, "slicecount", &count, 1);
  this->writeDataInt(fid, "beamletsize", &nbins, 1);
  this->writeDataInt(fid, "one4one", &tmp, 1);

  return;
}


int WriteBeamHDF5::writeSlice(hid_t fid, Beam *beam, int off)
{
  int nslice=beam->beam.size();
  
  int npart=0;
  vector<double> work;

  for (int islice=0;islice<nslice;islice++){
 
      int npart=beam->beam.at(islice).size();
      if (npart>work.size()){ work.resize(npart); }
      char slicename[20];
      sprintf(slicename,"slice%6.6d",islice+off+1);
      hid_t gid=H5Gcreate(fid,slicename,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

      
      this->writeDataDouble(gid, "current", &beam->current.at(islice), 1);

      for (int i=0; i<npart;i++){ work[i]=beam->beam.at(islice).at(i).gamma;}  
      this->writeDataDouble(gid, "gamma", &work[0], npart);

      for (int i=0; i<npart;i++){ work[i]=beam->beam.at(islice).at(i).theta;}  
      this->writeDataDouble(gid, "theta", &work[0], npart);

      for (int i=0; i<npart;i++){ work[i]=beam->beam.at(islice).at(i).x;}  
      this->writeDataDouble(gid, "x", &work[0], npart);

      for (int i=0; i<npart;i++){ work[i]=beam->beam.at(islice).at(i).y;}  
      this->writeDataDouble(gid, "y", &work[0], npart);

      for (int i=0; i<npart;i++){ work[i]=beam->beam.at(islice).at(i).px;}  
      this->writeDataDouble(gid, "px", &work[0], npart);

      for (int i=0; i<npart;i++){ work[i]=beam->beam.at(islice).at(i).py;}  
      this->writeDataDouble(gid, "py", &work[0], npart);


      H5Gclose(gid);
       

  }
  return nslice;
}



