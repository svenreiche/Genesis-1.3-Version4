#include "readFieldHDF5.h"

ReadFieldHDF5::ReadFieldHDF5(){
  isOpen=false;
  nwork=0;
}

ReadFieldHDF5::~ReadFieldHDF5()
{
  if (nwork>0){
   delete [] work;
  }
  if (isOpen){
   H5Fclose(fid);
  }  
}

bool ReadFieldHDF5::open(char *filename, int *ngrid, double *dgrid, double *lambda)
{

  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  fid=H5Fopen(filename,H5F_ACC_RDONLY,pid);
  H5Pclose(pid);


  isOpen=true;

  nwork=getDatasetSize(fid,(char *)"/slice000001/field-real");
  *ngrid=sqrt(nwork);
  readDataInt(fid,(char *)"slicecount",&count,1);
  readDataDouble(fid,(char *)"gridsize",dgrid,1);
  readDataDouble(fid,(char *)"wavelength",lambda,1);  
  readDataDouble(fid,(char *)"refposition",&s0,1);
  readDataDouble(fid,(char *)"slicespacing",&slen,1);



  // allocate some work arrays
  work=new double [nwork];

  return true; 
}


bool ReadFieldHDF5::readfield(double s, vector< complex<double> > *field){

  for(int i=0;i<nwork;i++){
      field->at(i)=complex<double> (0,0);
  }

  if(!isOpen){ return false; } // skip if partfile option is not selected
  
  double rslice=(s-s0)/slen;
  if (fabs(rslice-round(rslice))>1e-3){ return false; }
  int islice=static_cast<int> (round(rslice))+1;


  if ((islice<1)||(islice>count)){
    return true;
  }
   

  char name[20];


  sprintf(name,"slice%6.6d/field-real",islice);
  readDataDouble(fid,name,work,nwork);
  for (int i=0;i<nwork;i++){
    field->at(i)=work[i];
  }

  sprintf(name,"slice%6.6d/field-imag",islice);
  readDataDouble(fid,name,work,nwork);
  for (int i=0;i<nwork;i++){
    field->at(i)+=complex<double>(0,work[i]);
  }


  return true;
}


