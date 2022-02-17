
#include "Field.h"
//#include "genesis_fortran_common.h"


#include <fstream>

Field::~Field(){}

Field::Field(){
  harm=1;
  polarization=false;
  disabled=false;
  out_global=true;

  /* controlled by calls to Field::setOutput */
  doFFT=true;
  doSpatial=true;
  doIntensity=true;
  doDumpField=true;
}


void Field::init(int nsize, int ngrid_in, double dgrid_in, double xlambda0, double slicelen_in, double s0_in, int harm_in)
{  

  slicelength=slicelen_in;
  s0=s0_in;
  harm=harm_in;   // harmonics
  xlambda=xlambda0/static_cast<double>(harm);   // wavelength : xlambda0 is the reference wavelength of the fundamental
  ngrid=ngrid_in;
  gridmax=dgrid_in;
  dgrid=2*gridmax/static_cast<double>(ngrid-1); // grid point  separation
  xks=4.*asin(1)/xlambda;
  first=0;                                // pointer to slice which correspond to first in the time window
  dz_save=0;

  if (field.size()!=nsize){                  // allocate the memory in advance
     field.resize(nsize);
  }
  
  if (field[0].size()!=ngrid*ngrid){
    for (int i=0;i<nsize;i++){
        field[i].resize(ngrid*ngrid); 
    }
    solver.init(ngrid);
  } 

  return;
}


void Field::setStepsize(double delz)
{

  if (delz!=dz_save){
    dz_save=delz;
    //getdiag_(&delz,&dgrid,&xks,&ngrid);
  }

}

bool Field::getLLGridpoint(double x, double y, double *wx, double *wy, int *idx){

      if ((x>-gridmax) && (x < gridmax) && (y>-gridmax) && (y < gridmax)){
           *wx = (x+gridmax)/dgrid;
           *wy = (y+gridmax)/dgrid;
           int ix= static_cast<int> (floor(*wx));
           int iy= static_cast<int> (floor(*wy));
           *wx=1+floor(*wx)-*wx;
           *wy=1+floor(*wy)-*wy;
           *idx=ix+iy*ngrid;
           return true;
      } else {
           return false;
      }
}


void Field::track(double delz, Beam *beam, Undulator *und)
{
  solver.getDiag(delz,dgrid,xks,ngrid);  // check whether step size has changed and recalculate aux arrays
  solver.advance(delz,this,beam,und);
  return;
}


bool Field::subharmonicConversion(int harm_in, bool resample)
{
  // note that only an existing harmonic will be preserved
  harm=harm * harm_in;
  if ((!resample)|| (harm_in==1)){
    return true;
  }

  slicelength*=static_cast<double>(harm_in);
  int nsize=field.size();

  // step zero - calculate the radiation power, per slice

  vector<double> pow;
  pow.resize(nsize);
  for (int i=0;i<nsize;i++){
    pow[i]=0;
    for (int k=0; k<ngrid*ngrid;k++){
      complex<double> loc =field[i].at(k);  
      pow[i]+=loc.real()*loc.real()+loc.imag()*loc.imag();
    }
  }

 
  // step one - merge adjacent files into the first. 
  for (int i=0; i<nsize;i=i+harm_in){
    int is0= (i+first) % nsize; // loop over slices in the right order
    for (int j=1; j<harm_in; j++){
       int is1= (i+j+first) % nsize; // loop over slices to be merged
       for (int k=0; k<ngrid*ngrid;k++){
	 field[is0].at(k)+=field[is1].at(k);  // add field to slice
       }
       pow[is0]+=pow[is1];                  // add up power

    }
  }

  // step two - copy the merged one in the right order
  double scl=1./static_cast<double>(harm_in);

  int di=first % harm_in;
  for (int i=0; i<nsize/harm_in;i++){
     for (int k=0; k<ngrid*ngrid;k++){
	field[i].at(k)=field[i*harm_in+di].at(k);
     }
     pow[i]=pow[i*harm_in+di]*scl;    // calculate the normalized power
  }
  first=(first-di)/harm_in;

  // step three - average field by the subharmonic number
  field.resize(nsize/harm_in);
  for (int i=0; i<field.size();i++){
    double powloc=0;
    for (int k=0; k<ngrid*ngrid;k++){
       complex<double> loc =field[i].at(k);  
       powloc+=loc.real()*loc.real()+loc.imag()*loc.imag();
    }
    if (pow[i]>0){
     for (int k=0; k<ngrid*ngrid;k++){
       field[i].at(k)*=sqrt(pow[i]/powloc);  // do mean average because field were added up.
     }
    }

  }


  return true;
}


bool Field::harmonicConversion(int harm_in, bool resample)
{
  // note that only an existing harmonic will be preserved
  harm=1;
  if ((!resample)|| (harm_in==1)){
    return true;
  }

  slicelength/=static_cast<double>(harm_in);

  int nloc=field.size();
  field.resize(nloc*harm_in);
  for (int i=nloc;i<nloc*harm_in;i++){
    field.at(i).resize(ngrid*ngrid); // allocate space
  }

  // copy old slices into full field record
  for (int i=nloc;i>0;i--){
    for (int j=0;j<harm_in;j++){
      int idx=i*harm_in-1-j;
      for (int k=0; k<ngrid*ngrid;k++){
	field[idx].at(k)=field[(i-1)].at(k);
      }
    }
  }
  first*=harm_in; // adjust pointer of first slice
  return true;
}


