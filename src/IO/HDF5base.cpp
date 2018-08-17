#include "HDF5base.h"

extern bool MPISingle;

HDF5Base::HDF5Base(){}
HDF5Base::~HDF5Base(){}


//----------------------------
// routines taken from Output and writeHDF5



void HDF5Base::writeBuffer(hid_t gid, string dataset,vector<double> *data){


  // step 1 - calculate the file space and create dataset
  int size=1;
  if (!MPISingle){
     size=MPI::COMM_WORLD.Get_size(); // get size of cluster
  }

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
  if (!MPISingle){
     H5Pset_dxpl_mpio(pid,H5FD_MPIO_COLLECTIVE);    
  }
  H5Dwrite(did,H5T_NATIVE_DOUBLE,memspace,filespace,pid,&data->at(0));

  
  // close all HDF5 stuff 
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(pid);
  H5Dclose(did);

}

void HDF5Base::writeSingleNode(hid_t gid, string dataset,vector<double> *data){


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
 
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(did);
  

}

void HDF5Base::writeSingleNodeInt(hid_t gid, string dataset,vector<int> *data){


  int nd = data->size();

  hsize_t fblock[1]={nd};
  hid_t filespace=H5Screate_simple(1,fblock,NULL);
  hid_t did=H5Dcreate(gid,dataset.c_str(),H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);   
  H5Sclose(filespace);
  
  hid_t memspace=H5Screate_simple(1,fblock,NULL);
  filespace=H5Dget_space(did);

  if (s0==0){
    H5Dwrite(did,H5T_NATIVE_INT,memspace,filespace,H5P_DEFAULT,&data->at(0));

  }
 
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(did);
  

}




void HDF5Base::writeSingleNodeString(hid_t gid, string dataset, string *data){


 
  int nd = data->size();

  hsize_t fblock[1]={1};
  hid_t filespace=H5Screate_simple(1,fblock,NULL);


   hid_t dtype = H5Tcopy (H5T_C_S1);
   herr_t status = H5Tset_size (dtype, nd);



  hid_t did=H5Dcreate(gid,dataset.c_str(),dtype,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);   
  H5Sclose(filespace);


  hid_t memspace=H5Screate_simple(1,fblock,NULL);
  filespace=H5Dget_space(did);

  if (s0==0){
    H5Dwrite(did,dtype,memspace,filespace,H5P_DEFAULT,data->c_str());

  }
 
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(did);
  

}









//--------------------- 
// reading procedures


void HDF5Base::readDouble1D(hid_t fid, const char *name, double *data, hsize_t datsize, hsize_t start)
{



  int dataset_rank=1;

  hsize_t count[1] = {datsize};     // length of record entry 
  hsize_t offset[1] = {start};   // offset of record entry


  hid_t did=H5Dopen(fid,name,H5P_DEFAULT);

  hid_t filespace=H5Dget_space(did);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  hid_t memspace=H5Screate_simple(dataset_rank,count,NULL);


  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pid,H5FD_MPIO_COLLECTIVE);    


  H5Dread(did,H5T_NATIVE_DOUBLE,memspace,filespace,pid,data);

 
  // close all HDF5 stuff except for the file id fid
  H5Dclose(did);
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(pid);

  return;
}


void HDF5Base::readDataDouble(hid_t fid, char *name, double *data, int size)
{

  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dopen(fid,name,H5P_DEFAULT);
  hid_t plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);

  H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,plist_id,data);
  H5Dclose(dataset_id);     
  H5Sclose(dataspace_id);
  H5Pclose(plist_id);
  return;
}

void HDF5Base::readDataChar(hid_t fid, char *name, char *data, int size)
{

  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dopen(fid,name,H5P_DEFAULT);
  hid_t plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);

  H5Dread(dataset_id,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,plist_id,data);
  H5Dclose(dataset_id);     
  H5Sclose(dataspace_id);
  H5Pclose(plist_id);

  return;
}

void HDF5Base::readDataInt(hid_t fid, char *name, int *data, int size)
{
  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dopen(fid,name,H5P_DEFAULT);
  hid_t plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);

  H5Dread(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist_id,data);
  H5Dclose(dataset_id);     
  H5Sclose(dataspace_id);
  H5Pclose(plist_id);
  return;
}


int HDF5Base::getDatasetSize(hid_t fid, char *name)
{

  hsize_t dims[1],maxdims[1];
  hid_t  dsid=H5Dopen(fid,name,H5P_DEFAULT);
  hid_t spaceid=H5Dget_space(dsid);
  H5Sget_simple_extent_dims(spaceid,dims,maxdims);
  H5Dclose(dsid);
  return dims[0];

}


//------------------------
// simple read functions
bool HDF5Base::simpleReadDouble1D(const string &path, vector<double> *data){
  vector<string> ele;
  stringstream ss(path);
  string file;
  string group;
  char delim='/';             // does not compile it I use double quotation marks
  if (getline(ss,file,delim)){
     if (!getline(ss,group)){
       return false;
     }
  } else {
    return false;
  }
  

  hid_t fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  int nsize=this->getDatasetSize(fid,(char *)group.c_str());
  
  data->resize(nsize);
  this->readDataDouble(fid, (char *)group.c_str(), &data->at(0), nsize);
  H5Fclose(fid);
 
  return true;
}


//------------------------------
// utility functions

bool HDF5Base::checkForLink(hid_t fid, string name){

  htri_t result=H5Lexists(fid,name.c_str(),H5P_DEFAULT);
  if (result==true){
    return true;
  } 
  return false;
}

 
