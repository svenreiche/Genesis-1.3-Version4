
#include "writeBeamHDF5.h"

extern bool MPISingle;


// constructor destructor
WriteBeamHDF5::WriteBeamHDF5()
{
}

WriteBeamHDF5::~WriteBeamHDF5()
{
}


void WriteBeamHDF5::write(string fileroot, Beam *beam){


  size=MPI::COMM_WORLD.Get_size(); // get size of cluster
  rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node
  if (MPISingle){
    size=1;
    rank=0;
  }

  char filename[100];
  sprintf(filename,"%s.par.h5",fileroot.c_str()); 
  if (rank == 0) { cout << "Writing particle distribution to file: " <<filename << " ..." << endl;} 

  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  if (size>1){
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  fid=H5Fcreate(filename,H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
  H5Pclose(pid);

  s0=rank;
  int ntotal=size*beam->beam.size();

  // write global data

  this->writeGlobal(beam->nbins,beam->one4one,beam->reflength,beam->slicelength,beam->s0,ntotal);


  // write slices

  // loop through slices
  
  int smin=rank*beam->beam.size();
  int smax=smin+beam->beam.size();

  vector<double> work,cur;
  cur.resize(1);
  int nwork=0;
  int npart=0;

  for (int i=0; i<(ntotal);i++){
    s0=-1;
    char name[16];
    sprintf(name,"slice%6.6d",i+1);
    hid_t gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    int islice= i % beam->beam.size() ;   // count inside the given slice range

    if ((i>=smin) && (i<smax)){
      s0=0;    // select the slice which is writing
      npart=beam->beam.at(islice).size();
    }

    int root = i /beam->beam.size();  // the current rank which sends the informationof a slice
    if (size>1){
       MPI::COMM_WORLD.Bcast(&npart,1,MPI::INT, root);
    }

    if (npart != nwork){   // all cores do need to have the same length -> otherwise one4one crashes
	nwork=npart;
	work.resize(nwork);
    }

    if (s0==0) {
       for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip).gamma;}
    }
    this->writeSingleNode(gid,"gamma"  ,&work);

    if (s0==0) {
       for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip).theta;}
    }
    this->writeSingleNode(gid,"theta"  ,&work);

    if (s0==0) {
       for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip).x;}
    }
    this->writeSingleNode(gid,"x"  ,&work);

    if (s0==0) {
       for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip).y;}
    }
    this->writeSingleNode(gid,"y"  ,&work);

    if (s0==0) {
       for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip).px;}
    }
    this->writeSingleNode(gid,"px"  ,&work);

    if (s0==0) {
       for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip).py;}
    }
    this->writeSingleNode(gid,"py"  ,&work);

    if (s0==0){
      cur[0]=beam->current.at(islice);
    }
    this->writeSingleNode(gid,"current",&cur);

          
    H5Gclose(gid);
  }

  H5Fclose(fid);

 

  return;
}

void WriteBeamHDF5::writeGlobal(int nbins,bool one4one, double reflen, double slicelen, double s0, int count)
{



  vector<double> tmp;
  tmp.resize(1);

  tmp[0]=reflen;
  this->writeSingleNode(fid,"slicelength",&tmp);
  tmp[0]=slicelen;
  this->writeSingleNode(fid,"slicespacing",&tmp);
  tmp[0]=s0;
  this->writeSingleNode(fid,"refposition",&tmp);
  
  vector<int> itmp;
  itmp.resize(1);

  itmp[0]=nbins;
  this->writeSingleNodeInt(fid,"beamletsize",&itmp);
  itmp[0]=count;
  this->writeSingleNodeInt(fid,"slicecount",&itmp);
  itmp[0]=static_cast<int>(one4one);
  this->writeSingleNodeInt(fid,"one4one",&itmp);
  
  return;
}





