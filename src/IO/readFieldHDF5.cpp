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
}


void ReadFieldHDF5::close(){
  if (isOpen){ H5Fclose(fid); }
}
  

bool ReadFieldHDF5::readGlobal(int rank, int size,string file, Setup *setup, Time *time, double offset, bool dotime)
{


  isOpen=false;
  // read global data

  double reflen,slen,gridsize;
  int count,ngrid;


  fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);  
  readDataDouble(fid,(char *)"refposition",&s0,1);
  readDataDouble(fid,(char *)"gridsize",&gridsize,1);
  readDataDouble(fid,(char *)"wavelength",&reflen,1);
  readDataDouble(fid,(char *)"slicespacing",&slicelen,1);
  readDataInt(fid,(char *)"slicecount",&count,1);
  readDataInt(fid,(char *)"gridpoints",&ngrid,1);
  isOpen=true;

  nwork=ngrid*ngrid;
  work = new double [nwork];

  

  double lambda=setup->getReferenceLength();                        // reference length for theta
  
  s0=s0+offset;  // add offset from input deck
  slen=slicelen*count;
  return true;
} 



/*
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
*/

