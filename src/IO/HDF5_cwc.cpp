#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>
#include <hdf5.h>

#include "HDF5_cwc.h"

using namespace std;

extern bool MPISingle;

/*
 * Write code with full flexibility (datatype is double)
 * . ndim is the dimensionality of the data
 * . totalsize is the size of the object in the file
 * . offset and count describe what this process writes
 */
HDF5_CollWriteCore::HDF5_CollWriteCore()
{
	ndim_=0;
	did_=0;
	did_valid_=false;
}
HDF5_CollWriteCore::~HDF5_CollWriteCore()
{
	if(did_valid_) {
		cout << "warning: HDF5 obj still open" << endl;
		close();
	}
}

void HDF5_CollWriteCore::create_and_prepare(hid_t gid, string dataset, string unit, vector<hsize_t> *totalsize, unsigned long ndim)
{
  if (totalsize->size()!=ndim)
  {
     cout << "obj size: incorrect element count" << endl;
     return;
  }

  // step 1 - create dataset
  hid_t filespace=H5Screate_simple(ndim,&(*totalsize)[0],NULL);
  did_=H5Dcreate(gid,dataset.c_str(),H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);   
  H5Sclose(filespace);

  // write attribute
  hid_t aid = H5Screate(H5S_SCALAR);
  hid_t atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, unit.size());
  H5Tset_strpad(atype,H5T_STR_NULLTERM);
  hid_t attr = H5Acreate2(did_,"unit",atype,aid,H5P_DEFAULT,H5P_DEFAULT);
  H5Awrite(attr,atype,unit.c_str());
  H5Sclose(aid);
  H5Tclose(atype);
  H5Aclose(attr);

  ndim_=ndim;
  did_valid_=true;
}

void HDF5_CollWriteCore::write(vector<double> *data, vector<hsize_t> *count, vector<hsize_t> *offset)
{
  bool do_write = true;
  double *pdata = NULL;
  hsize_t ntot_xfer;

  if(!did_valid_)
  {
    cout << "error: no valid HDF5 object is assigned to this instance" << endl;
    return;
  }
  if ((count->size()!=ndim_) && (offset->size()!=ndim_))
  {
    cout << "offset/count vectors: incorrect element count" << endl;
    return;
  }


  do_write=true;
  if(data==NULL)
    do_write=false;

  for(int k=0; k<count->size(); k++)
    if(count->at(k)==0)
      do_write=false;
  

  // step 2 - file space
  hid_t memspace=H5Screate_simple(ndim_,&(*count)[0],NULL);

  // step 3 - set up hyperslab for file transfer.
  hid_t filespace=H5Dget_space(did_);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,&(*offset)[0],NULL,&(*count)[0],NULL);

  pdata=NULL;
  if(do_write) {
    pdata = &data->at(0);
  } else {
    // No data is transferred by the local process ...
    // ... but it still has to call H5Dwrite
    H5Sselect_none(memspace);
    H5Sselect_none(filespace);
  }


  // step 4 - set up transfer and write
  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  if (!MPISingle){
    H5Pset_dxpl_mpio(pid,H5FD_MPIO_COLLECTIVE);    
  }
  H5Dwrite(did_,H5T_NATIVE_DOUBLE,memspace,filespace,pid,pdata);
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(pid);
}

void HDF5_CollWriteCore::close(void)
{
  if(!did_valid_)
    return;

  // close all HDF5 stuff 
  H5Dclose(did_);
  did_=0;
  did_valid_=false;
}

