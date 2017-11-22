#include "Output.h"

Output::Output(){}
Output::~Output(){}

void Output::close(){


  H5Fclose(fid);

  
}


void Output::open(const char * file, int s0_in, int ds_in)
{
  // ns = total number of slices
  // nz = total number of integration steps
  // s0 = slice number of first slice for a given node
  // ds = number of slices, which are kept by a given node

  // set record range and allocate working memory
  s0=s0_in;
  ds=ds_in;
  
  // create the file for parallel access
  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  fid=H5Fopen(file,H5F_ACC_RDWR, pid); 
  H5Pclose(pid);

}



void Output::writeBeamBuffer(Beam *beam)
{


  // step 1 - create the group
  hid_t gid;
  gid=H5Gcreate(fid,"Beam",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  // step 2 - write individual datasets
  this->writeBuffer(gid, "energy",&beam->gavg);
  this->writeBuffer(gid, "energyspread",&beam->gsig);
  this->writeBuffer(gid, "xposition",&beam->xavg);
  this->writeBuffer(gid, "yposition",&beam->yavg);
  this->writeBuffer(gid, "pxposition",&beam->pxavg);
  this->writeBuffer(gid, "pyposition",&beam->pyavg);
  this->writeBuffer(gid, "xsize",&beam->xsig);
  this->writeBuffer(gid, "ysize",&beam->ysig);
  this->writeBuffer(gid, "bunching",&beam->bunch);
  this->writeBuffer(gid, "bunchingphase",&beam->bphi);
  
  this->writeBuffer(gid, "betax",&beam->bx);
  this->writeBuffer(gid, "betay",&beam->by);
  this->writeBuffer(gid, "alphax",&beam->ax);
  this->writeBuffer(gid, "alphay",&beam->ay);
  this->writeBuffer(gid, "emitx",&beam->ex);
  this->writeBuffer(gid, "emity",&beam->ey);
  this->writeBuffer(gid, "current",&beam->cu);
  // step 3 - close group and done

  H5Gclose(gid);

  // step 4 - write z-posiiton for plotting
  gid=H5Gopen(fid,"Global",H5P_DEFAULT);
  this->writeSingleNode(gid,"zplot",&beam->zpos);
  H5Gclose(gid);

  return;
}


void Output::writeFieldBuffer(Field *field)
{


  // step 1 - create the group
  hid_t gid;
  char name[10];

  int harm=field->harm;
  if (harm==1){
     sprintf(name,"Field");
  } else {
     sprintf(name,"Field%d",harm);
  }
  gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  // step 2 - write individual datasets
  this->writeBuffer(gid, "power",&field->power);
  this->writeBuffer(gid, "xsize",&field->xsig);
  this->writeBuffer(gid, "ysize",&field->ysig);
  this->writeBuffer(gid, "xposition",&field->xavg);
  this->writeBuffer(gid, "yposition",&field->yavg);
  this->writeBuffer(gid, "intensity-nearfield",&field->nf_intensity);
  this->writeBuffer(gid, "phase-nearfield",&field->nf_phi);
  this->writeBuffer(gid, "intensity-farfield",&field->ff_intensity);
  this->writeBuffer(gid, "phase-farfield",&field->ff_phi);


  // step 3 - close group and done

  H5Gclose(gid);
  return;

 

}


void Output::writeBuffer(hid_t gid, string dataset,vector<double> *data){


  // step 1 - calculate the file space and create dataset
  int size=MPI::COMM_WORLD.Get_size(); // get size of cluster


  hsize_t dz=data->size()/ds;
  hsize_t fblock[2]={dz,size*ds};
  hid_t filespace=H5Screate_simple(2,fblock,NULL);
  hid_t did=H5Dcreate(gid,dataset.c_str(),H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);   
  H5Sclose(filespace);



  // step 2 - file space
  hsize_t count[2]={dz,ds};
  hsize_t offset[2] = {0,s0};   // offset of record entry
  hid_t memspace=H5Screate_simple(2,count,NULL);


  // step 3 - set up hyperslab for file transfer.
  filespace=H5Dget_space(did);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);



  // step 4 - set up transfer and write
  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pid,H5FD_MPIO_COLLECTIVE);    
  H5Dwrite(did,H5T_NATIVE_DOUBLE,memspace,filespace,pid,&data->at(0));

  
  // close all HDF5 stuff 
  H5Dclose(did);
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(pid);

}

void Output::writeSingleNode(hid_t gid, string dataset,vector<double> *data){


  int nd = data->size();

  hsize_t fblock[1]={nd};
  hid_t filespace=H5Screate_simple(1,fblock,NULL);
  hid_t did=H5Dcreate(gid,dataset.c_str(),H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);   
  H5Sclose(filespace);
  
  hid_t memspace=H5Screate_simple(1,fblock,NULL);
  filespace=H5Dget_space(did);

  if (s0==0){
    H5Dwrite(did,H5T_NATIVE_DOUBLE,memspace,filespace,H5P_DEFAULT,&data->at(0));

  }
  /*
  hsize_t block[1],count[1],offset[1],stride[1];
  int dataset_rank=1;
  if (s0==0) {
    block[0] =nd; // chunck size
    count[0] = 1;     // repeatition of block 
    offset[0]= 0;   // offset of record entry
    stride[0]= 1;
  } else{
    block[0] = 0; // chunck size
    count[0] = 0;     // repeatition of block 
    offset[0]= 0;   // offset of record entry
    stride[0]= 1;
  }

  //  cout << offset[1] << " " << stride[1] << " "  << count[1] << " " << block[1] << endl; 

  
  // set up memory and file space - filespac eneeds to by a hyperslab
  hid_t memspace=H5Screate_simple(dataset_rank,block,NULL);
  filespace=H5Dget_space(did);

 
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,stride,count,block);
  //H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  
  // set up transfer and write
  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pid,H5FD_MPIO_INDEPENDENT);    
  H5Dwrite(did,H5T_NATIVE_DOUBLE,memspace,filespace,pid,data);
 
  // close all HDF5 stuff except for the file id fid
  H5Pclose(pid);
  */
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(did);
  

}





