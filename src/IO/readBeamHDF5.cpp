#include "readBeamHDF5.h"

ReadBeamHDF5::ReadBeamHDF5()
{
  isOpen=false;
  nwork=-1;
}

ReadBeamHDF5::~ReadBeamHDF5()
{
  if (nwork>0){ delete [] work; }
}

void ReadBeamHDF5::close(){
  if (isOpen){ H5Fclose(fid); }
}
  

bool ReadBeamHDF5::readGlobal(int rank, int size,string file, Setup *setup, Time *time, bool dotime)
{


  isOpen=false;
  // read global data

  double reflen;
  int nbins,one4one;


  fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);  
  readDataDouble(fid,(char *)"refposition",&s0,1);
  readDataDouble(fid,(char *)"slicelength",&reflen,1);
  readDataDouble(fid,(char *)"slicespacing",&slicelen,1);
  readDataInt(fid,(char *)"slicecount",&count,1);
  readDataInt(fid,(char *)"beamletsize",&nbins,1);
  readDataInt(fid,(char *)"one4one",&one4one,1);
  isOpen=true;


  double lambda=setup->getReferenceLength();                        // reference length for theta
  
  slen=slicelen*count;       // window of the imported distribution                 



  if (fabs((lambda-reflen)/lambda)>1e-6){
      if (rank==0){ cout << "*** Error: Mismatch between reference length of run and of input file" << endl;}
      return false;    
  }


  if (nbins!=setup->getNbins()){
      if (rank==0){ cout << "*** Error: Mismatch between beamlet size (NBINS) of run and of input file" << endl;}
      return false;    
  }

  
  if (((setup->getOne4One())&&(one4one==0)) || ((!setup->getOne4One())&&(one4one!=0))){
      if (rank==0){ cout << "*** Error: Mismatch between ONE4ONE flag of run and of input file" << endl;}
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




bool ReadBeamHDF5::readSlice(double s, vector<Particle> *slice, double *current, bool one4one){

  
  slice->resize(0);
  *current=0;

 
  if(!isOpen){ return false; } // skip if partfile option is not selected



  double rslice=(s-s0)/slicelen;

  if (fabs(rslice-round(rslice))>1e-3){ return false; }
  int islice=static_cast<int> (round(rslice))+1;

  double wei=1.;
  if ((islice<1) || (islice> count)){
    if (one4one){  return false; }   // one4one do nothing
    wei=0;                      // fill slice form 1st slice but put all values to zero
    islice=1;
  }


  // get size of data set from slice
  char name[20];
  sprintf(name,"slice%6.6d/gamma",islice);

  int nsize=getDatasetSize(fid, name);

  if (nsize>nwork){ // allocate extra work array to hold field
    if (nwork>0) {delete [] work;}
    nwork=nsize;
    work=new double [nwork];
  }


  // allocate data in the beam record
  slice->resize(nsize);

  // get current   
  sprintf(name,"slice%6.6d/current",islice);
  readDataDouble(fid,name,current,1);
  *current*=wei;

  sprintf(name,"slice%6.6d/gamma",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).gamma=work[i];
  }
  
  sprintf(name,"slice%6.6d/theta",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).theta=work[i];
  }

  sprintf(name,"slice%6.6d/x",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).x=work[i]*wei;
  }

  sprintf(name,"slice%6.6d/y",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).y=work[i]*wei;
  }

  sprintf(name,"slice%6.6d/px",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).px=work[i]*wei;
  }

  sprintf(name,"slice%6.6d/py",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).py=work[i]*wei;
  }
  
  return true;
}


