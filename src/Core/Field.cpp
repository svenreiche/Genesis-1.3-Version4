
#include "Field.h"
//#include "genesis_fortran_common.h"


#include <fstream>



Field::~Field(){
#ifdef FFTW                // release the FFTW plan
  fftw_destroy_plan(p);
#endif

  delete[] in;
  delete[] out;
}

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
  dgrid=2*gridmax/static_cast<double>(ngrid-1); // grid pointe separation

  if (field.size()!=nsize){                  // allocate the memory in advance
     field.resize(nsize);
  }
  
  if (field[0].size()!=ngrid*ngrid){
    for (int i=0;i<nsize;i++){
        field[i].resize(ngrid*ngrid); 
    } 
  } 

  in = new complex<double> [ngrid*ngrid];
  out= new complex<double> [ngrid*ngrid];

#ifdef FFTW
  p  = fftw_plan_dft_2d(ngrid,ngrid,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out),FFTW_FORWARD,FFTW_MEASURE);
#endif

  xks=4.*asin(1)/xlambda;
  first=0;                                // pointer to slice which correspond to first in the time window
  dz_save=0;

  
  return;
}


// at each run the buffer should be cleared.
void Field::initDiagnostics(int nz)
{
  idx=0;
  accuslip=0;    
  int ns=field.size();


  power.resize(ns*nz);
  if (doSpatial){
    xavg.resize(ns*nz);
    xsig.resize(ns*nz);
    yavg.resize(ns*nz);
    ysig.resize(ns*nz);
  }
#ifdef FFTW
  if (doFFT){  // output of FFT is requested
    txavg.resize(ns*nz);
    txsig.resize(ns*nz);
    tyavg.resize(ns*nz);
    tysig.resize(ns*nz);
  } else {
    txavg.resize(0);
    txsig.resize(0);
    tyavg.resize(0);
    tysig.resize(0);
  }
#endif
  nf_intensity.resize(ns*nz);
  nf_phi.resize(ns*nz);
  ff_intensity.resize(ns*nz);
  ff_phi.resize(ns*nz);
  if (out_global){
    energy.resize(nz);
    if (doSpatial){
      gl_xsig.resize(nz);
      gl_xavg.resize(nz);
      gl_ysig.resize(nz);
      gl_yavg.resize(nz);
    } else {
      gl_xsig.resize(0);
      gl_xavg.resize(0);
      gl_ysig.resize(0);
      gl_yavg.resize(0);
    }
    gl_nf_intensity.resize(nz);
    gl_ff_intensity.resize(nz);
#ifdef FFTW
    if (doFFT){
      gl_txsig.resize(nz);
      gl_txavg.resize(nz);
      gl_tysig.resize(nz);
      gl_tyavg.resize(nz);
    } else {
      gl_txsig.resize(0);
      gl_txavg.resize(0);
      gl_tysig.resize(0);
      gl_tyavg.resize(0);
    }
#endif
  }

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



//  diagnostic part


void Field::diagnostics(bool output)
{

  if (!output) { return; }

  double shift=-0.5*static_cast<double> (ngrid);
  complex<double> loc;
  double loc2[2];
  int ds=field.size();
  int ioff=idx*ds;

  double acc_power=0;
  double acc_x=0;
  double acc_x2=0;
  double acc_y=0;
  double acc_y2=0;
  double acc_nf=0;
  double acc_ff=0;
#ifdef FFTW
  double acc_tpower=0;
  double acc_tx=0;
  double acc_tx2=0;
  double acc_ty=0;
  double acc_ty2=0;
#endif

  double ks=4.*asin(1)/xlambda;
  double scl=dgrid*eev/ks;
  double scltheta=xlambda/ngrid/dgrid;

  for (int is=0; is < ds; is++){
    int islice= (is+first) % ds ;   // include the rotation due to slippage
 
    double bpower=0;
    double bxavg=0;
    double byavg=0;
    double bxsig=0;
    double bysig=0;
    double fpower=0;
    double fxavg=0;
    double fyavg=0;
    double fxsig=0;
    double fysig=0;
    double bintensity=0;   
    double bfarfield=0;
    double bphinf=0;
    double bphiff=0;
    complex<double> ff = complex<double> (0,0);


    for (int iy=0;iy<ngrid;iy++){
      double dy=static_cast<double>(iy)+shift;
      for (int ix=0;ix<ngrid;ix++){
        double dx=static_cast<double>(ix)+shift;
	int i=iy*ngrid+ix;
        loc=field.at(islice).at(i);
	in[i]=loc;   // field for the FFT
        double wei=loc.real()*loc.real()+loc.imag()*loc.imag();
        ff+=loc;
        bpower+=wei;
        bxavg+=dx*wei;
        bxsig+=dx*dx*wei;
        byavg+=dy*wei;
        bysig+=dy*dy*wei;
      }
    }


    if (out_global){
      acc_power+=bpower;
      acc_x +=bxavg;
      acc_x2+=bxsig;  
      acc_y +=byavg;
      acc_y2+=bysig;
    }
    
    if (bpower>0){
      bxavg/=bpower;
      bxsig/=bpower;
      byavg/=bpower;
      bysig/=bpower;
    }
    bxsig=sqrt(abs(bxsig-bxavg*bxavg));
    bysig=sqrt(abs(bysig-byavg*byavg));


    // calculate the divergence angle
#ifdef FFTW
    if (doFFT){
      fftw_execute(p);
      for (int iy=0;iy<ngrid;iy++){
	double dy=static_cast<double>(iy)+shift;
	for (int ix=0;ix<ngrid;ix++){
	  double dx=static_cast<double>(ix)+shift;
	  int iiy=(iy+(ngrid+1)/2) % ngrid;
	  int iix=(ix+(ngrid+1)/2) % ngrid;
	  int ii=iiy*ngrid+iix;
	  loc=out[ii];
	  double wei=loc.real()*loc.real()+loc.imag()*loc.imag();
	  fpower+=wei;
	  fxavg+=dx*wei;
	  fxsig+=dx*dx*wei;
	  fyavg+=dy*wei;
	  fysig+=dy*dy*wei;
	}
      }

      if (out_global){
	acc_tpower+=fpower;
	acc_tx +=fxavg;
	acc_tx2+=fxsig;   // bxsig and bysig are in this part still variances
	acc_ty +=fyavg;
	acc_ty2+=fysig;
      }
 

      if (fpower>0){
	  fxavg/=fpower;
	  fxsig/=fpower;
	  fyavg/=fpower;
	  fysig/=fpower;
      }

      fxsig=sqrt(abs(fxsig-fxavg*fxavg));
      fysig=sqrt(abs(fysig-fyavg*fyavg));
    }

#endif





    bfarfield=ff.real()*ff.real()+ff.imag()*ff.imag();
    if (bfarfield > 0){
      bphiff=atan2(ff.imag(),ff.real());
    }
    int i=(ngrid*ngrid-1)/2;
    loc=field.at(islice).at(i);
    bintensity=loc.real()*loc.real()+loc.imag()*loc.imag();
    if (bintensity > 0) {
      bphinf=atan2(loc.imag(),loc.real()); 
    }
    if (out_global){
      acc_ff+=bfarfield;
      acc_nf+=bintensity;
    }


    bpower*=scl*scl/vacimp; // scale to W
    bxavg*=dgrid;
    bxsig*=dgrid;
    byavg*=dgrid;
    bysig*=dgrid;
    bintensity*=eev*eev/ks/ks/vacimp;  // scale to W/m^2

    power[ioff+is]=bpower;
    if (doSpatial) {
      xavg[ioff+is] =bxavg;
      xsig[ioff+is] =bxsig;
      yavg[ioff+is] =byavg;
      ysig[ioff+is] =bysig;
    }
#ifdef FFTW
    if (doFFT){
      txavg[ioff+is]=fxavg*scltheta;
      txsig[ioff+is]=fxsig*scltheta;
      tyavg[ioff+is]=fyavg*scltheta;
      tysig[ioff+is]=fysig*scltheta;
    }
#endif
    nf_intensity[ioff+is]=bintensity;
    nf_phi[ioff+is]=bphinf;
    ff_intensity[ioff+is]=bfarfield;
    ff_phi[ioff+is]=bphiff;

  } 
  
  // accumulate all data fromt eh cores
  double temp=0;
  int size;      
  if(out_global) {
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      if (size>1){
	MPI_Allreduce(&acc_power, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_power=temp;
	MPI_Allreduce(&acc_x, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_x=temp;
	MPI_Allreduce(&acc_x2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_x2=temp;
	MPI_Allreduce(&acc_y, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_y=temp;
	MPI_Allreduce(&acc_y2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_y2=temp;
	MPI_Allreduce(&acc_nf, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_nf=temp;
	MPI_Allreduce(&acc_ff, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_ff=temp;
#ifdef FFTW
	if (doFFT){
	  MPI_Allreduce(&acc_tpower, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  acc_tpower=temp;
	  MPI_Allreduce(&acc_tx, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  acc_tx=temp;
	  MPI_Allreduce(&acc_tx2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  acc_tx2=temp;
	  MPI_Allreduce(&acc_ty, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  acc_ty=temp;
	  MPI_Allreduce(&acc_ty2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  acc_ty2=temp;
	}
#endif
      }

      energy[idx]=acc_power*scl*scl/vacimp*slicelength/3e8;  // energy 
      if (acc_power == 0){
	acc_power=1.;
      }
      if (doSpatial){
        gl_xavg[idx]=acc_x/acc_power;
        gl_xsig[idx]=acc_x2/acc_power;
        gl_yavg[idx]=acc_y/acc_power;
        gl_ysig[idx]=acc_y2/acc_power;
        gl_xsig[idx]=sqrt(abs(gl_xsig[idx]-gl_xavg[idx]*gl_xavg[idx]));
        gl_ysig[idx]=sqrt(abs(gl_ysig[idx]-gl_yavg[idx]*gl_yavg[idx]));
        gl_xavg[idx]*=dgrid;
        gl_xsig[idx]*=dgrid;
        gl_yavg[idx]*=dgrid;
        gl_ysig[idx]*=dgrid;
      }
#ifdef FFTW
      if (doFFT){
	if (acc_tpower == 0){
	  acc_tpower=1.;
	}
	gl_txavg[idx]=acc_tx/acc_tpower;
	gl_txsig[idx]=acc_tx2/acc_tpower;
	gl_tyavg[idx]=acc_ty/acc_tpower;
	gl_tysig[idx]=acc_ty2/acc_tpower;
	gl_txsig[idx]=sqrt(abs(gl_txsig[idx]-gl_txavg[idx]*gl_txavg[idx]));
	gl_tysig[idx]=sqrt(abs(gl_tysig[idx]-gl_tyavg[idx]*gl_tyavg[idx]));
	gl_txavg[idx]*=scltheta;
	gl_txsig[idx]*=scltheta;
	gl_tyavg[idx]*=scltheta;
	gl_tysig[idx]*=scltheta;
      }
#endif
      double norm=1./double(size*field.size());
      gl_nf_intensity[idx]=acc_nf*norm*eev*eev/ks/ks/vacimp;
      gl_ff_intensity[idx]=acc_ff*norm;
  }
  
  idx++;
  
}
