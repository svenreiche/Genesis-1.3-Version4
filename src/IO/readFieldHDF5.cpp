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
  

bool ReadFieldHDF5::readGlobal(int rank, int size,string file, Setup *setup, Time *time, int harm, bool dotime)
{


  isOpen=false;
  // read global data

  double reflen;


  fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);  
  readDataDouble(fid,(char *)"refposition",&s0,1);
  readDataDouble(fid,(char *)"gridsize",&dgrid,1);
  readDataDouble(fid,(char *)"wavelength",&reflen,1);
  readDataDouble(fid,(char *)"slicespacing",&slicelen,1);
  readDataInt(fid,(char *)"slicecount",&count,1);
  readDataInt(fid,(char *)"gridpoints",&ngrid,1);
  isOpen=true;

  nwork=ngrid*ngrid;
  work = new double [nwork];

  

  double lambda=setup->getReferenceLength();                        // reference length for theta
  
  s0=s0;  // add offset from input deck
  slen=slicelen*count;



  double ks=4.*asin(1)/reflen;
  scl=dgrid*eev/ks/sqrt(vacimp);
  scl=1./scl;
  dgrid=0.5*dgrid*static_cast<double>(ngrid-1);


  if (fabs((lambda-reflen)/lambda)>1e-6){
      if (rank==0){ cout << "*** Error: Mismatch between reference length of run and of input file" << endl;}
      return false;    
  }



  // set-up time window if not set.
  if ((!time->isInitialized())&& (count >1)){
      time->setSampleRate(slicelen/reflen); 
      time->setTimeWindowStart(s0);
      time->setTimeWindowLength(slen);

      map<string,string> arg;
      if (dotime) { arg["time"]="true"; } else { arg["time"]="false"; }
      bool check=time->init(rank, size, &arg, setup);
      if (!check){ return check; }
  }



  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  if (fabs(slicelen-lambda*sample)/slicelen>1e-6){
      if (rank==0){ cout << "*** Error: Mismatch in sample rate of run and of input file" << endl;}
      return false;
  }








  return true;
} 




bool ReadFieldHDF5::readSlice(double s, vector<complex<double> >*field){

  for(int i=0;i<nwork;i++){
      field->at(i)=complex<double> (0,0);
  }

  if(!isOpen){ return false; } // skip if partfile option is not selected
  
  double rslice=(s-s0)/slicelen;
  if (fabs(rslice-round(rslice))>1e-3){ return false; }
  int islice=static_cast<int> (round(rslice))+1;


  if ((islice<1)||(islice>count)){
    return true;
  }
   

  char name[30];


  sprintf(name,"slice%6.6d/field-real",islice);
  readDataDouble(fid,name,work,nwork);
  for (int i=0;i<nwork;i++){
    field->at(i)=scl*work[i];
  }

  sprintf(name,"slice%6.6d/field-imag",islice);
  readDataDouble(fid,name,work,nwork);
  for (int i=0;i<nwork;i++){
    field->at(i)+=complex<double>(0,scl*work[i]);
  }


  return true;
}


