::::::::::::::
src/Core/Beam.cpp
::::::::::::::
#include "Beam.h"
#include "Field.h"
#include "Sorting.h"

#ifdef VTRACE
#include "vt_user.h"
#endif

Beam::~Beam(){}
Beam::Beam(){}

void Beam::init(int nsize, int nbins_in, double reflen_in, double slicelen_in, double s0_in, bool one4one_in )
{  

  nbins=nbins_in;
  reflength=reflen_in;  // the length corresponding to 2pi in ponderomotive phase.
  slicelength=slicelen_in;  // reflength times samplerate.
  s0=s0_in;
  one4one=one4one_in;

  current.resize(nsize);
  beam.resize(nsize);

  return;
}




void Beam::initDiagnostics(int nz)
{
  
  idx=0;
  int ns=current.size();
  zpos.resize(nz);
  xavg.resize(nz*ns);
  xsig.resize(nz*ns);
  yavg.resize(nz*ns);
  ysig.resize(nz*ns);
  gavg.resize(nz*ns);
  gsig.resize(nz*ns);
  pxavg.resize(nz*ns);
  pyavg.resize(nz*ns);
  bunch.resize(nz*ns); 
  bphi.resize(nz*ns);

  bx.resize(ns);
  by.resize(ns);
  ax.resize(ns);
  ay.resize(ns);
  ex.resize(ns);
  ey.resize(ns);
  cu.resize(ns);
  
  bh.clear();
  ph.clear();
  if (bharm>1){
    bh.resize(bharm-1);
    ph.resize(bharm-1);
    for (int i=0;i<(bharm-1);i++){
      bh[i].resize(nz*ns);
      ph[i].resize(nz*ns);
    }
  }

}

// initialize the sorting routine
// reference position is in ponderomotive phase. The valid bucket size is from 0 to 2 pi*sample
void Beam::initSorting(int rank,int size,bool doshift,bool dosort)
{
  int isz=beam.size();
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.init(rank,size,doshift,dosort);
  sorting.configure(0,sl,0,sl*isz,0,sl*isz,false); 
  return;
}

int Beam::sort()
{
  int shift=0;
  double dQ=ce/slicelength;
  if (one4one){
    shift= sorting.sort(&beam);
    for (int i=0; i<beam.size();i++){    // correct the local current
      int np=beam.at(i).size();
      current.at(i)=static_cast<double>(np)*ce/slicelength;
    }
  }  
  return shift;
}


void Beam::track(double delz,vector<Field *> *field, Undulator *und){

#ifdef VTRACE
  VT_TRACER("Beam_Tracking");
#endif  

  
  for (int i=0; i<field->size();i++){
    field->at(i)->setStepsize(delz);
  }

  solver.track(delz*0.5,this,und,false);   // track transverse coordinates first half of integration step

  solver.advance(delz,this,field,und);     // advance longitudinal variables 

  incoherent.apply(this,und,delz);         // apply effect of incoherent synchrotron

  solver.applyR56(this,und,reflength);    // apply the longitudinalphase shift due to R56 if a chicane is selected.

  solver.track(delz*0.5,this,und,true);      // apply corrector settings and track second half for transverse coordinate
  return;


}

bool Beam::harmonicConversion(int harmonic, bool resample)
{
  if ((resample==true) && (one4one==false)) { return false;}  // resampling requires one4one
 
  reflength=reflength/static_cast<double>(harmonic);
  if (resample){ 
    slicelength=slicelength/static_cast<double>(harmonic);
  }
  for (int i=0; i<beam.size(); i++){
    for (int j=0; j<beam[i].size();j++){
      beam[i].at(j).theta*=static_cast<double>(harmonic); 
    }
  }
  if (!resample) { return true; }

  // blowing up the slice number
  int nsize=beam.size();

  beam.resize(harmonic*nsize);
  current.resize(harmonic*nsize);
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.configure(0,sl,0,sl*nsize*harmonic,0,sl*nsize*harmonic,false);   


  // first make all new slices zero length
  for (int i=nsize; i < harmonic*nsize; i++){
    beam.at(i).resize(0);    
    current.at(i)=0;
  }
  // second copy the old slices 0,1,.. n-1 to h*0,h*1,... h*(n-1)

  Particle p;
  for (int i=nsize-1; i>0; i--){  // runs down to 1 only because there is no need to copy from 0 to h*0 = 0
    int nloc=beam[i].size();
    for (int j=0; j<nloc; j++){
      p.gamma=beam[i].at(j).gamma;
      p.theta=beam[i].at(j).theta;
      p.x    =beam[i].at(j).x;
      p.y    =beam[i].at(j).y;
      p.px   =beam[i].at(j).px;
      p.py   =beam[i].at(j).px;
      beam[i*harmonic].push_back(p);
    }
    beam[i].clear(); 
  }

  // updating the sorting algorithm

  int shift=this->sort();  // sort the particles and update current
  return true;
}

bool Beam::subharmonicConversion(int harmonic, bool resample)
{
  if ((resample==true) && (one4one==false)) { return false;}
   
 
  reflength=reflength*static_cast<double>(harmonic);
  if (resample){ // needs a lot of working here............
    slicelength=slicelength*static_cast<double>(harmonic);
  }


  for (int i=0; i<beam.size(); i++){
    for (int j=0; j<beam[i].size();j++){
      beam[i].at(j).theta/=static_cast<double>(harmonic);   // preparing to push everything into first slice
    }
  }
  if (!resample) { return true; }

  


// prepare to copy everyting into the first slice
  int nsize=beam.size();
  Particle p;

  if ((nsize % harmonic) !=0) { return false;}  // check whether the number of slices cannot merged into a smaller number

  //return true;

  double dtheta=4.*asin(1)*slicelength/reflength/static_cast<double>(harmonic);

  for (int i=1; i<nsize;i++){
      for (int k=0; k<beam.at(i).size();k++){
	p.gamma=beam[i].at(k).gamma;
	p.theta=beam[i].at(k).theta+i*dtheta;
	p.x    =beam[i].at(k).x;
	p.y    =beam[i].at(k).y;
	p.px   =beam[i].at(k).px;
	p.py   =beam[i].at(k).py;
	beam[0].push_back(p);
      }
      beam[i].clear(); 
  }
  beam.resize(nsize/harmonic);
  current.resize(nsize/harmonic);

  // updating the sorting algorithm
  int isz=beam.size();
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.configure(0,sl,0,sl*isz,0,sl*isz,false); 


  int shift=this->sort();  // sort the particles and update current
  return true;
}



void Beam::diagnostics(bool output, double z)
{
  if (!output) { return; }
  

  zpos[idx]=z;

  int ds=beam.size();
  int ioff=idx*ds; 

  for (int is=0; is < ds; is++){
    double bgavg=0;
    double bgsig=0;
    double bxavg=0;
    double bxsig=0;
    double byavg=0;
    double bysig=0;
    double bpxavg=0;
    double bpyavg=0;
    double br=0;
    double bi=0;
    double bbavg=0;
    double bbphi=0;

    unsigned int nsize=beam.at(is).size();
    for (int i=0;i < nsize;i++){
      double xtmp=beam.at(is).at(i).x;
      double ytmp=beam.at(is).at(i).y;
      double pxtmp=beam.at(is).at(i).px;
      double pytmp=beam.at(is).at(i).py;
      double gtmp=beam.at(is).at(i).gamma;
      double btmp=beam.at(is).at(i).theta;
      bxavg +=xtmp;
      byavg +=ytmp;
      bpxavg+=pxtmp;
      bpyavg+=pytmp;
      bgavg +=gtmp;
      bgsig +=gtmp*gtmp;
      bxsig +=xtmp*xtmp;
      bysig +=ytmp*ytmp;
      br+=cos(btmp);
      bi+=sin(btmp);
    }

    double scl=1;
    if (nsize>0){
      scl=1./static_cast<double>(nsize);
    }
    bgavg*=scl;
    bxavg*=scl;
    byavg*=scl;
    bgsig*=scl;
    bxsig*=scl;
    bysig*=scl;
    bpxavg*=scl;
    bpyavg*=scl;
    bbavg=sqrt(bi*bi+br*br)*scl;
    bbphi=atan2(bi,br);
    bgsig=sqrt(fabs(bgsig-bgavg*bgavg));
    bxsig=sqrt(fabs(bxsig-bxavg*bxavg));
    bysig=sqrt(fabs(bysig-byavg*byavg));
    
    gavg[ioff+is]=bgavg;
    gsig[ioff+is]=bgsig;
    xavg[ioff+is]=bxavg;
    xsig[ioff+is]=bxsig;
    yavg[ioff+is]=byavg;
    ysig[ioff+is]=bysig;
    pxavg[ioff+is]=bpxavg;
    pyavg[ioff+is]=bpyavg;
    bunch[ioff+is]=bbavg;
    bphi[ioff+is]=bbphi;
    for (int ih=1; ih<bharm;ih++){   // calculate the harmonics of the bunching
      br=0;
      bi=0;
      for (int i=0;i < nsize;i++){
        double btmp=static_cast<double>(ih+1)*beam.at(is).at(i).theta;
        br+=cos(btmp);
        bi+=sin(btmp);
      }
      bh[ih-1][ioff+is]=sqrt(bi*bi+br*br)*scl;
      ph[ih-1][ioff+is]=atan2(bi,br);
    }
  }
  idx++;
}


void Beam::diagnosticsStart()
{
  double gx,gy,gammax,gammay;
  double x1,y1,x2,y2,px1,py1,px2,py2,g1,xpx,ypy;

  int ds=beam.size();


  for (int is=0; is<ds;is++){
    cu[is]=current[is];
    x1=0;
    x2=0;
    px1=0;
    px2=0;
    y1=0;
    y2=0;
    py1=0;
    py2=0;
    xpx=0;
    ypy=0;
    g1=0;
    unsigned int nsize=beam.at(is).size();
    for (int i=0;i < nsize;i++){
      double xtmp=beam.at(is).at(i).x;
      double ytmp=beam.at(is).at(i).y;
      double pxtmp=beam.at(is).at(i).px;
      double pytmp=beam.at(is).at(i).py;
      double gtmp=beam.at(is).at(i).gamma;
      x1+=xtmp;
      y1+=ytmp;
      px1+=pxtmp;
      py1+=pytmp;
      g1+=gtmp; 
      x2+=xtmp*xtmp;
      y2+=ytmp*ytmp;
      px2+=pxtmp*pxtmp;
      py2+=pytmp*pytmp;
      xpx+=xtmp*pxtmp;
      ypy+=ytmp*pytmp;
    }
    double norm=1.;
    if (nsize>0){ norm=1./static_cast<double>(nsize);}
    x1*=norm;
    x2*=norm;
    y1*=norm;
    y2*=norm;
    px1*=norm;
    px2*=norm;
    py1*=norm;
    py2*=norm;
    g1*=norm;
    xpx*=norm;
    ypy*=norm;
    
    // because genesis works with momenta and not divergence, the emittance does not need energy
    ex[is]=sqrt(fabs((x2-x1*x1)*(px2-px1*px1)-(xpx-x1*px1)*(xpx-x1*px1)));
    ey[is]=sqrt(fabs((y2-y1*y1)*(py2-py1*py1)-(ypy-y1*py1)*(ypy-y1*py1)));
    bx[is]=(x2-x1*x1)/ex[is]*g1;
    by[is]=(y2-y1*y1)/ey[is]*g1;
    ax[is]=-(xpx-x1*px1)*g1/ex[is];
    ay[is]=-(ypy-y1*py1)*g1/ey[is];
     
    gx=(1+ax[is]*ax[is])/bx[is];
    gy=(1+ay[is]*ay[is])/by[is];

  }


  return;
}

::::::::::::::
src/Core/BeamSolver.cpp
::::::::::::::
#include "BeamSolver.h"
#include "Field.h"
#include "Beam.h"

BeamSolver::BeamSolver()
{
  onlyFundamental=false;
}

BeamSolver::~BeamSolver(){}


void BeamSolver::advance(double delz, Beam *beam, vector< Field *> *field, Undulator *und)
{
   
  // here the harmonics needs to be taken into account

  vector<int> nfld;
  vector<double> rtmp;
  rpart.clear();
  rharm.clear();
  xks=1;  // default value in the case that no field is defined

   for (int i=0; i < field->size(); i++){
    int harm=field->at(i)->getHarm();
    if ((harm==1) || !onlyFundamental){
      xks=field->at(i)->xks/static_cast<double>(harm);    // fundamental field wavenumber used in ODE below
      nfld.push_back(i);
      rtmp.push_back(und->fc(harm)/field->at(i)->xks);      // here the harmonics have to be taken care
      rpart.push_back(0);
      rharm.push_back(static_cast<double>(harm));
    }
  }  


  xku=und->getku();
  if (xku==0){   // in the case of drifts - the beam stays in phase if it has the reference energy // this requires that the phase slippage is not applied
    xku=xks*0.5/und->getGammaRef()/und->getGammaRef();
  }
	    
  double aw=und->getaw();

  double autophase=und->autophase();


  // Runge Kutta solver to advance particle

  

  for (int is=0; is<beam->beam.size(); is++){    


      for (int ip=0; ip<beam->beam.at(is).size();ip++){
        gamma=beam->beam.at(is).at(ip).gamma;
        theta=beam->beam.at(is).at(ip).theta+autophase; // add autophase here
        double x =beam->beam.at(is).at(ip).x;
        double y =beam->beam.at(is).at(ip).y;
        double px=beam->beam.at(is).at(ip).px;
        double py=beam->beam.at(is).at(ip).py;
	double awloc=und->faw(x,y);                 // get the transverse dependence of the undulator field
        btpar=1+px*px+py*py+aw*aw*awloc*awloc;	  

	ez=0;

	cpart=0;
	double wx,wy;
	int idx;
        for (int ifld=0;ifld<nfld.size();ifld++){

	  int islice=(is+field->at(nfld[ifld])->first) % field->at(nfld[ifld])->field.size(); 

	  if (field->at(nfld[ifld])->getLLGridpoint(x,y,&wx,&wy,&idx)){ // check whether particle is on grid
           cpart=field->at(nfld[ifld])->field[islice].at(idx)*wx*wy;
           idx++;
           cpart+=field->at(nfld[ifld])->field[islice].at(idx)*(1-wx)*wy;
           idx+=field->at(nfld[ifld])->ngrid-1;
           cpart+=field->at(nfld[ifld])->field[islice].at(idx)*wx*(1-wy);
           idx++;
           cpart+=field->at(nfld[ifld])->field[islice].at(idx)*(1-wx)*(1-wy);
           rpart[ifld]=rtmp[ifld]*awloc*conj(cpart);
	  }
	}
	this->RungeKutta(delz);

        beam->beam.at(is).at(ip).gamma=gamma;
        beam->beam.at(is).at(ip).theta=theta; 
      }
    
  }
  return;
}

void BeamSolver::RungeKutta(double delz)
{
  // Runge Kutta Solver 4th order - taken from pushp from the old Fortran source


  // first step
  k2gg=0;
  k2pp=0;

  this->ODE(gamma,theta);

  // second step
  double stpz=0.5*delz;

  gamma+=stpz*k2gg;
  theta+=stpz*k2pp;
  
  k3gg=k2gg;
  k3pp=k2pp;

  k2gg=0;
  k2pp=0;

  this->ODE(gamma,theta);

  // third step
  gamma+=stpz*(k2gg-k3gg);
  theta+=stpz*(k2pp-k3pp);

  k3gg/=6;
  k3pp/=6;

  k2gg*=-0.5;
  k2pp*=-0.5;

  this->ODE(gamma,theta);

  // fourth step
  stpz=delz;

  gamma+=stpz*k2gg;
  theta+=stpz*k2pp;

  k3gg-=k2gg;
  k3pp-=k2pp;

  k2gg*=2;
  k2pp*=2;

  this->ODE(gamma,theta);
  gamma+=stpz*(k3gg+k2gg/6.0);
  theta+=stpz*(k3pp+k2pp/6.0);

  return;
}


void BeamSolver::ODE(double tgam,double tthet)
{

  // differential equation for longitudinal motion
  double ztemp1=-2./xks;
  complex<double> ctmp=0;
  for (int i=0; i<rpart.size();i++){
    ctmp+=rpart[i]*complex<double> (cos(rharm[i]*tthet), -sin(rharm[i]*tthet));
  }
  double btper0=btpar+ztemp1*ctmp.real();   //perpendicular velocity
  double btpar0=sqrt(1.-btper0/(tgam*tgam));     //parallel velocity
  k2pp+=xks*(1.-1./btpar0)+xku;             //dtheta/dz
  k2gg+=ctmp.imag()/btpar0/tgam-ez;         //dgamma/dz

  return; 
}

::::::::::::::
src/Core/Collective.cpp
::::::::::::::
// include all routine which are arising from a collective contribution from the beam. This are:
// Wake
// Longrange Space Charge
// CSR

#include "Collective.h"

#include <iostream>
#include <fstream>

Collective::Collective(){
  doWakes=false;
  doSpaceCharge=false;
  doCSR=false;

// default values, should come from init
  r=2.5e-3;
  sigma=5.813e7;
  tau=8.1e-6;
  roundpipe=false;
  gap=0.5e-3;
  lgap=5;

  hrough=100e-9;
  lrough=100e-6; 

// setting up current (flat)
  ds=1e-10;
double Q=200e-12;
double currentmax=3000;
double smax=Q/currentmax*3e8;
ns=static_cast<int>(round(smax/ds));

current.resize(ns);
wakeres.resize(ns);
wakegeo.resize(ns);

for (int i=0;i<ns;i++){
current[i]=currentmax;
wakeres[i]=0;
wakegeo[i]=0;
}
}

Collective::~Collective(){}


void Collective::WakeRou()
{
  double pi=2*asin(1.);
  double kappa=2*pi/lrough;

  double ra=r*hrough*hrough*kappa*kappa*kappa/8;  
  double coef=2*ra/pi/r/r;

}

void Collective::WakeGeo(){

  double pi=2*asin(1.);
  double coef=-vacimp*ce/(pi*pi*r*lgap)*2*sqrt(0.5*gap); // scaling coefficient

  if (!roundpipe) { coef*=0.956; }     

  for (int is = 0;is<ns;is++){
     wakegeo[is]=coef*sqrt(ds*is);     
  }
}

void Collective::WakeRes(){

   // physical constants

   double c=3e8;
   double pi=2*asin(1.);



   // step 1 - clear wake;
   for (int i=0;i<ns;i++){
       wakeres[i]=0;
   }
   if (sigma<=0) { return; }


   // step 2 - calculate the impedance
   double s0=pow(2*r*r/vacimp/sigma,1./3.); // characteristic length in SI units
   double gamma=tau/s0;
   double coef = r/(s0*s0);

   int nk=1024;
   vector<double> Zre,Zim;
   Zre.resize(nk);
   Zim.resize(nk);
   Zre[0]=0;
   Zim[0]=0;

   double kappamax=6;  // kappa is k*s0!!! and chosen to resolve the resonance quite well.
   double coef2=-kappamax/nk/s0*c/pi*(vacimp*ce/4/pi); 
   
   if (roundpipe){ 
     for (int i=1; i<nk;i++){
       double kappa=(i)*kappamax/nk;  // value of kappa	 
       double t = kappa/sqrt(1+kappa*kappa*gamma*gamma);
       double lambdaRe=coef*sqrt(t)*sqrt(1.-t*gamma);
       double lambdaIm=coef*sqrt(t)*sqrt(1.+t*gamma)-kappa*kappa*r*0.5/s0/s0;
       double nomi=2.*kappa/(c*r*s0)/(lambdaRe*lambdaRe+lambdaIm*lambdaIm)*coef2;
       Zre[i]=lambdaRe*nomi;   
       Zim[i]=-lambdaIm*nomi;   
     }
   }else{
     int nq = 8*nk;
     vector<double> coh,sih;
     coh.resize(nq);
     sih.resize(nq);   //  this is actually sinh(q)/q
     double dq=15/static_cast<double>(nq-1);
      for (int i =1; i<nq;i++){
	   coh[i]=0.5*(exp(dq*i)+exp(-dq*i));
	   sih[i]=coh[i]-exp(-dq*i);
	   sih[i]/=dq*i;
     }
     coh[0]=1.;
     sih[0]=1.;   
     for (int i =1; i<nk; i++){       
        double kappa=(i+1.)*kappamax/nk;  // value of kappa	 
 	double t = kappa/sqrt(1+kappa*kappa*gamma*gamma);
        double scale=2.*15.*kappa/(c*r*s0*(2*nq-1));
	Zre[i]=0;
	Zim[i]=0;	 
        // integrate over q-> infty which is actually exp(30)
        for (int j=1;j<nq;j++){
               double lambdaRe=coef*sqrt(t)*sqrt(1.-t*gamma)*coh[j]*coh[j];
    	       double lambdaIm=coef*sqrt(t)*sqrt(1.+t*gamma)*coh[j]*coh[j]-kappa*kappa*r*0.5/s0/s0*sih[j]*coh[j];
	       double nomi=scale/(lambdaRe*lambdaRe+lambdaIm*lambdaIm)*coef2;
	       Zre[i]+=lambdaRe*nomi;
	       Zim[i]+=-lambdaIm*nomi;	 
	 }
     }
  }

  // step3 - construct the single particle wake 
  // using a recursive algorithm for sin(np + p) = sin(np)*cos(p)+cos(np)*sin(p)

  for (int is=0; is< ns; is++){
    double phi0=kappamax/nk/s0*ds*is;
    double cphi0=cos(phi0);
    double sphi0=sin(phi0);
    double cphin=1;
    double sphin=0;
     
    for (int ik=1 ; ik <nk; ik++){ // starts at ik=1 because Z[0]=0
      double sphin1=cphin*sphi0+sphin*cphi0;
      double cphin1=cphin*cphi0-sphin*sphi0;
      wakeres[is]+=Zre[ik]*cphin1+Zim[ik]*sphin1;
      cphin=cphin1;
      sphin=sphin1;
    }
  }
}


::::::::::::::
src/Core/Control.cpp
::::::::::::::
#include <sstream>
#include "Control.h"
#include "writeFieldHDF5.h"
#include "writeBeamHDF5.h"


#ifdef VTRACE
#include "vt_user.h"
#endif



Control::Control()
{
  nwork=0;
}


Control::~Control()
{
}


bool Control::applyMarker(Beam *beam, vector<Field*>*field, Undulator *und)
{

  bool sort=false;

  int marker=und->getMarker();
  // possible file names
  int istepz=und->getStep();
  stringstream sroot;
  sroot << "." << istepz;

  if ((marker & 1) != 0){
    WriteFieldHDF5 dump;
    dump.write(root+sroot.str(),field);
  }
  
  if ((marker & 2) != 0){
    WriteBeamHDF5 dump;
    dump.write(root+sroot.str(),beam);
  }
  
  if ((marker & 4) != 0){
    sort=true;   // sorting is deferred after the particles have been pushed by Runge-Kutta
  }

  // bit value 8 is checked in und->advance()
  

  return sort;
}


void Control::output(Beam *beam, vector<Field*> *field, Undulator *und)
{


  
  Output *out=new Output;

  string file=root.append(".out.h5");
  out->open(file,noffset,nslice);
  
  out->writeGlobal(und->getGammaRef(),reflen,sample,slen,one4one,timerun,scanrun);
  out->writeLattice(beam,und);



  for (int i=0; i<field->size();i++){
       out->writeFieldBuffer(field->at(i));
  }

  out->writeBeamBuffer(beam);
  out->close();
 
  delete out;
  return;



}


bool Control::init(int inrank, int insize, const char *file, Beam *beam, vector<Field*> *field, Undulator *und, bool inTime, bool inScan)
{

  rank=inrank;
  size=insize;
  accushift=0;

  stringstream sroot(file);
  root=sroot.str();
  root.resize(root.size()-7);  // remove the extension ".h5"

  one4one=beam->one4one;
  reflen=beam->reflength;
  sample=beam->slicelength/reflen;

  timerun=inTime;
  scanrun=inScan;
 
 

  // cross check simulation size

  nslice=beam->beam.size();
  noffset=rank*nslice;
  ntotal=size*nslice;  // all cores have the same amount of slices

  slen=ntotal*sample*reflen;


  if (rank==0){
    if(scanrun) { 
       cout << "Scan run with " << ntotal << " slices" << endl; 
    } else {
       if(timerun) { 
         cout << "Time-dependent run with " << ntotal << " slices" << " for a time window of " << slen*1e6 << " microns" << endl; 
       } else { 
         cout << "Steady-state run" << endl;
       }
    }
  }
  


  // initial diagnostic

  if (rank==0) { cout << "Initial analysis of electron beam and radiation field..."  << endl; }

  beam->initDiagnostics(und->outlength());
  beam->diagnostics(true,0);
  beam->diagnosticsStart();
  for (int i=0; i<field->size();i++){
      field->at(i)->initDiagnostics(und->outlength());
      field->at(i)->diagnostics(true);  // initial values
  }	

  return true;  
}



void Control::applySlippage(double slippage, Field *field)
{

#ifdef VTRACE
  VT_TRACER("Slippage");
#endif  

  if (timerun==false) { return; }

 
  // update accumulated slippage
  accushift+=slippage;

  // allocate working space

  if(nwork<field->ngrid*field->ngrid*2){
    nwork=field->ngrid*field->ngrid*2;
    work=new double [nwork];
  } 
  
  MPI::Status status;

  // following routine is applied if the required slippage is alrger than 80% of the sampling size

  int direction=1;


  while(abs(accushift)>(sample*0.8)){
      // check for anormal direction of slippage (backwards slippage)
      if (accushift<0) {direction=-1;} 

      accushift-=sample*direction; 

      // get adjacent node before and after in chain
      int rank_next=rank+1;
      int rank_prev=rank-1;
      if (rank_next >= size ) { rank_next=0; }
      if (rank_prev < 0 ) { rank_prev = size-1; }	

      // for inverse direction swap targets
      if (direction<0) {
	int tmp=rank_next;
        rank_next=rank_prev;
        rank_prev=tmp; 
      }

      int tag=1;
   
      // get slice which is transmitted
      int last=(field->first+field->field.size()-1)  %  field->field.size();
      // get first slice for inverse direction
      if (direction<0){
	last=(last+1) % field->field.size();  //  this actually first because it is sent backwards
      }

      if (size>1){
        if ( (rank % 2)==0 ){                   // even nodes are sending first and then receiving field
           for (int i=0; i<nwork/2; i++){
	     work[2*i]  =field->field[last].at(i).real();
	     work[2*i+1]=field->field[last].at(i).imag();
	   }
	   MPI::COMM_WORLD.Send(work, nwork, MPI::DOUBLE, rank_next, tag);
	   MPI::COMM_WORLD.Recv(work, nwork, MPI::DOUBLE, rank_prev, tag, status);
	   for (int i=0; i<nwork/2; i++){
	     complex <double> ctemp=complex<double> (work[2*i],work[2*i+1]);
	     field->field[last].at(i)=ctemp;
	   }
	} else {                               // odd nodes are receiving first and then sending

	  MPI::COMM_WORLD.Recv(work, nwork, MPI::DOUBLE, rank_prev, tag, status);
	  for (int i=0; i<nwork/2; i++){
	    complex <double> ctemp=complex<double> (work[2*i],work[2*i+1]);
	    work[2*i]  =field->field[last].at(i).real();
	    work[2*i+1]=field->field[last].at(i).imag();
	    field->field[last].at(i)=ctemp;
	  }
	  MPI::COMM_WORLD.Send(work, nwork, MPI::DOUBLE, rank_next, tag);
	}
      }

      // first node has emptz field slipped into the time window
      if ((rank==0) && (direction >0)){
        for (int i=0; i<nwork/2;i++){
	  field->field[last].at(i)=complex<double> (0,0);
        }
      }

      if ((rank==(size-1)) && (direction <0)){
        for (int i=0; i<nwork/2;i++){
	  field->field[last].at(i)=complex<double> (0,0);
        }
      }

      // last was the last slice to be transmitted to the succeding node and then filled with the 
      // the field from the preceeding node, making it now the start of the field record.
      field->first=last;
      if (direction<0){
	field->first=(last+1) % field->field.size();
      }
  }

}
::::::::::::::
src/Core/EFieldSolver.cpp
::::::::::::::
#include "EFieldSolver.h"

EFieldSolver::EFieldSolver(){
  nz=0;
  nphi=0;
  ngrid_ref=0;

}

EFieldSolver::~EFieldSolver(){}

void EFieldSolver::init(double rmax_in, int ngrid_in, int nz_in, int nphi_in, double lambda){

  rmax_ref=rmax_in;
  ngrid_ref=ngrid_in;
  nz=nz_in;
  nphi=nphi_in;
  ks = 4*asin(1)/lambda;
}

void EFieldSolver::shortRange(vector<Particle> *beam,vector<double> &ez, double current, double gammaz){

  double npart=beam->size();
  double rmax=rmax_ref;
  int ngrid=ngrid_ref;

  double pi=2.*asin(1);

  // gammaz = gamma/sqrt(1+aw^2)
  double gz2=gammaz*gammaz;
  double xcuren=current;



  if ((nz==0)||(npart==0)||(ngrid<2)) {return;}
  
  // adjust working arrays
  if (npart>idx.size()){
    idx.resize(npart);
    azi.resize(npart);
  }
  

  // calculate center of beam slice
  double xcen=0;
  double ycen=0;
  for (int i=0; i<npart; i++){
    xcen+=beam->at(i).x;
    ycen+=beam->at(i).y; 
  }
  xcen/=static_cast<double>(npart);
  ycen/=static_cast<double>(npart);

  // calculate radial grid point position

  double dr=rmax/(ngrid-1);
  double rbound=0;

  for (int i=0; i<npart; i++){
    double tx=beam->at(i).x-xcen;
    double ty=beam->at(i).y-ycen;
    double r=sqrt(tx*tx + ty*ty);
    if (r>rbound) { rbound=r; }
    idx[i]=static_cast<int> (floor(r/dr));
    azi[i]=atan2(ty,tx);
  }
  
  // check for beam sizes not fitting in the radial grid.
  if (2*rbound>rmax){
    ngrid=static_cast<int> (floor(2*rbound/dr))+1;
    rmax=(ngrid-1)*dr;
  }

  // adjust working arrays
  if (ngrid!=csrc.size()){
   vol.resize(ngrid);
   lmid.resize(ngrid);
   llow.resize(ngrid);
   lupp.resize(ngrid);
   celm.resize(ngrid);
   csrc.resize(ngrid);
   clow.resize(ngrid);
   cmid.resize(ngrid);
   cupp.resize(ngrid);
   gam.resize(ngrid);
  }



  // r_j are the center grid points. The boundaries are given +/-dr/2. Thus r_j=(j+1/2)*dr
  for (int j=0;j<ngrid;j++){
    vol[j]=dr*dr*(2*j+1);
    lupp[j]=2*(j+1)/vol[j];  // Eq 3.30 of my thesis
    llow[j]=2*j/vol[j];      // Eq 3.29
    if (j==0){
      rlog[j]=0;
    }else{
      rlog[j]=2*log((j+1)/j)/vol[j];  // Eq. 3.31
    }
  }
  lupp[ngrid-1]=0;



  for (int l=1;l<nz+1;l++){                 // loop over longitudinal modes
       for (int i=1; i<ngrid; i++){
           lmid[i]=-llow[i]-lupp[i]-l*l*ks*ks/gz2; // construct the diagonal terms for m=0
           csrc[i]=complex<double> (0,0);           // clear source term 
       }
       

       // normalize source term       
       double coef=vacimp/eev*xcuren/static_cast<double>(npart)/pi;

       // add the azimuthal dependence
       for (int m=-nphi;m<=nphi;m++){
	 // construct the source term

	 for (int i=0;i<ngrid;i++){
	   csrc[i]=complex<double> (0,0);
	   cmid[i]=lmid[i]-rlog[i]*m*m;
	 }

         for (int i=0; i<npart; i++){
	   double phi=l*beam->at(i).theta+m*azi[i];	   
	   csrc[idx[i]]+=complex<double>(cos(phi),-sin(phi));
         }

	 for (int i=0;i<ngrid;i++){
	   csrc[i]*=complex<double> (0,coef/vol[i]);
	 }

	 // solve tridiag(clow,cmid,cupp,csrc,celm,ngrid);

	 complex<double> bet=cmid[0];
	 celm[0]=csrc[0]/bet;
	 for (int i=1;i<ngrid;i++){
	   gam[i]=cupp[i-1]/bet;
	   bet=cmid[i]-clow[i]*gam[i];
	   celm[i]=(csrc[i]-clow[i]*celm[i-1])/bet;
	 }
	 for (int i=ngrid-2;i>-1;i--){
	   celm[i]=celm[i]-gam[i+1]*celm[i+1];
	 }


         // calculate the field at the electron bposition
	 for (int i=0;i<npart;i++){
	    double phi=l*beam->at(i).theta+m*azi[i];	   
	    complex<double> ctmp=celm[idx[i]]*complex<double>(cos(phi),sin(phi));
	    ez[i]+=2*ctmp.real();
	 }
       }
       
  }
  
  return;

}
::::::::::::::
src/Core/Field.cpp
::::::::::::::

#include "Field.h"
#include "genesis_fortran_common.h"

#ifdef VTRACE
#include "vt_user.h"
#endif


Field::~Field(){}
Field::Field(){
  harm=1;
  polarization=false;
  disabled=false;
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
  xavg.resize(ns*nz);
  xsig.resize(ns*nz);
  yavg.resize(ns*nz);
  ysig.resize(ns*nz);
  nf_intensity.resize(ns*nz);
  nf_phi.resize(ns*nz);
  ff_intensity.resize(ns*nz);
  ff_phi.resize(ns*nz);

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
  
#ifdef VTRACE
  VT_TRACER("Field_Tracking");
#endif  


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
  int ds=field.size();
  int ioff=idx*ds;

  for (int is=0; is < ds; is++){
    int islice= (is+first) % ds ;   // include the rotation due to slippage
 
    double bpower=0;
    double bxavg=0;
    double byavg=0;
    double bxsig=0;
    double bysig=0;
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
        double wei=loc.real()*loc.real()+loc.imag()*loc.imag();
        ff+=loc;
        bpower+=wei;
        bxavg+=dx*wei;
        bxsig+=dx*dx*wei;
        byavg+=dy*wei;
        bysig+=dy*dy*wei;
      }
    }
    
    if (bpower>0){
      bxavg/=bpower;
      bxsig=sqrt(abs(bxsig/bpower-bxavg*bxavg));
      byavg/=bpower;
      bysig=sqrt(abs(bysig/bpower-byavg*byavg));
    }

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

    double ks=4.*asin(1)/xlambda;
    double scl=dgrid*eev/ks;
    bpower*=scl*scl/vacimp; // scale to W
    bxavg*=dgrid;
    bxsig*=dgrid;
    byavg*=dgrid;
    bysig*=dgrid;
    bintensity*=eev*eev/ks/ks/vacimp;  // scale to W/m^2

    power[ioff+is]=bpower;
    xavg[ioff+is] =bxavg;
    xsig[ioff+is] =bxsig;
    yavg[ioff+is] =byavg;
    ysig[ioff+is] =bysig;
    nf_intensity[ioff+is]=bintensity;
    nf_phi[ioff+is]=bphinf;
    ff_intensity[ioff+is]=bfarfield;
    ff_phi[ioff+is]=bphiff;

  }
  
  idx++;
  
}
::::::::::::::
src/Core/FieldSolver.cpp
::::::::::::::
#include "FieldSolver.h"
#include "Field.h"
#include "Beam.h"

FieldSolver::FieldSolver()
{
  delz_save=0;
}

FieldSolver::~FieldSolver(){}



void FieldSolver::advance(double delz, Field *field, Beam *beam, Undulator *und)
{
 
  for (int ii=0; ii<field->field.size();ii++){  // ii is index for the beam
    int i= (ii+field->first) % field->field.size();           // index for the field
  
    // clear source term
    for (int i=0; i< ngrid*ngrid; i++){
      crsource[i]=0;
    }

    int harm=field->getHarm();
    // construc source term
    if (und->inUndulator()&& field->isEnabled()) { 
      double scl=und->fc(harm)*vacimp*beam->current[ii]*field->xks*delz;
       scl/=4*eev*beam->beam[ii].size()*field->dgrid*field->dgrid;
       complex<double> cpart;
       double part,weight,wx,wy;
       int idx;

       for (int ip=0;ip<beam->beam.at(ii).size();ip++){
	 double x    =beam->beam.at(ii).at(ip).x;
	 double y    =beam->beam.at(ii).at(ip).y;
	 double theta=static_cast<double>(harm)*beam->beam.at(ii).at(ip).theta;
	 double gamma=beam->beam.at(ii).at(ip).gamma;

         if (field->getLLGridpoint(x,y,&wx,&wy,&idx)){

           part=sqrt(und->faw2(x,y))*scl/gamma;
           // tmp  should be also normalized with beta parallel
           cpart=complex<double>(sin(theta),cos(theta))*part;

           weight=wx*wy;
	   crsource[idx]+=weight*cpart;
           weight=(1-wx)*wy;
           idx++;
	   crsource[idx]+=weight*cpart;
           weight=wx*(1-wy);
           idx+=ngrid-1;
           crsource[idx]+=weight*cpart;
           weight=(1-wx)*(1-wy);
           idx++;
	   crsource[idx]+=weight*cpart;

	 }
       } 
    }  // end of source term construction

    this->ADI(field->field[i]);
  } 	  

  return;
}


void FieldSolver::ADI(vector<complex<double> > &crfield)
{
  int ix,idx;
  // implicit direction in x
  for (idx=0;idx<ngrid;idx++){
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+ngrid]-2.0*crfield[idx]);
  }
  for (idx=ngrid;idx<ngrid*(ngrid-1);idx++){
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+ngrid]-2.0*crfield[idx]+crfield[idx-ngrid]);
  }
  for (idx=ngrid*(ngrid-1);idx<ngrid*ngrid;idx++){
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx-ngrid]-2.0*crfield[idx]);
  }

  // solve tridiagonal system in x
  this->tridagx(crfield);

  // implicit direction in y
  for(ix=0;ix<ngrid*ngrid;ix+=ngrid){
    idx=ix;
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+1]-2.0*crfield[idx]);
    for(idx=ix+1;idx<ix+ngrid-1;idx++){
      r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+1]-2.0*crfield[idx]+crfield[idx-1]);
    }
    idx=ix+ngrid-1;
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx-1]-2.0*crfield[idx]);
  }

  // solve tridiagonal system in y
  this->tridagy(crfield);


  return;
}


void FieldSolver::tridagx(vector<complex<double > > &u)
{ 
  for (int i=0;i<ngrid*(ngrid-1);i+=ngrid){
    u[i]=r[i]*cbet[0];
    for (int k=1; k<ngrid;k++){
      u[k+i]=(r[k+i]-c[k]*u[k+i-1])*cbet[k];
    }
    for (int k=ngrid-2;k>=0;k--){
      u[k+i]-=cwet[k+1]*u[k+i+1];
    }
  }
  return;
}

void FieldSolver::tridagy(vector<complex<double > > &u)
{
  for (int i=0; i<ngrid; i++){
    u[i]=r[i]*cbet[0];
  }
  for (int k=1;k<ngrid-1;k++){
    int n=k*ngrid;
    for (int i=0;i<ngrid;i++){
      u[n+i]=(r[n+i]-c[k]*u[n+i-ngrid])*cbet[k];
    }
  }
  for (int k=ngrid-2;k>=0;k--){
    int n=k*ngrid;
    for (int i=0; i<ngrid; i++){
      u[n+i]-=cwet[k+1]*u[n+i+ngrid];
    }
  }

  return;
}


void FieldSolver::getDiag(double delz,double dgrid, double xks, int ngrid_in)
{

  if (delz==delz_save){
    return;
  }
  delz_save=delz;
  ngrid=ngrid_in;


  double rtmp=0.25*delz/(xks*dgrid*dgrid); //factor dz/(4 ks dx^2)
  cstep = complex<double> ( 0, rtmp );

  double *mupp = new double[ngrid];
  double *mmid = new double[ngrid];
  double *mlow = new double[ngrid];
  complex<double> *cwrk1= new complex<double>[ngrid];
  complex<double> *cwrk2= new complex<double>[ngrid];
  if (c.size()!=ngrid){
    c.resize(ngrid);
    r.resize(ngrid*ngrid);
    cbet.resize(ngrid);
    cwet.resize(ngrid);
    crsource.resize(ngrid*ngrid);
  }

  mupp[0]=rtmp;
  mmid[0]=-2*rtmp;
  mlow[0]=0;
  for (int i=1;i<(ngrid-1);i++){
    mupp[i]=rtmp;
    mmid[i]=-2*rtmp;
    mlow[i]=rtmp;
  }
  mupp[ngrid-1]=0;
  mmid[ngrid-1]=-2*rtmp;
  mlow[ngrid-1]=rtmp;

  for (int i=0; i <ngrid; i++){
    cwrk1[i] =complex<double>(0,-mupp[i]);
    cwrk2[i] =complex<double>(1,-mmid[i]);
    c[i]     =complex<double>(0,-mlow[i]);
  }

 
  cbet[0]=1./cwrk2[0];
  cwet[0]=0.;
  for (int i=1; i<ngrid; i++){
    cwet[i]=cwrk1[i-1]*cbet[i-1];
    cbet[i]=1./(cwrk2[i]-c[i]*cwet[i]); 
    
  }

  delete[] mupp;
  delete[] mmid;
  delete[] mlow;
  delete[] cwrk1;
  delete[] cwrk2;
}
::::::::::::::
src/Core/Gencore.cpp
::::::::::::::
#include "Gencore.h"

extern bool MPISingle;

#ifdef VTRACE
#include "vt_user.h"
#endif


int Gencore::run(const char *file, Beam *beam, vector<Field*> *field, Undulator *und,bool isTime, bool isScan)
{


        //-------------------------------------------------------
        // init MPI and get size etc.
        //
#ifdef VTRACE
  VT_TRACER("Core");
#endif  
        int size=1;
        int rank=0;

	if (!MPISingle){
           size=MPI::COMM_WORLD.Get_size(); // get size of cluster
           rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node
        }

	if (rank==0) {
          cout << endl << "Running Core Simulation..." << endl;
        }

        //-----------------------------------------
	// init beam, field and undulator class

        Control   *control=new Control;

	control->init(rank,size,file,beam,field,und,isTime,isScan); 


        //------------------------------------------
        // main loop
	       	
	while(und->advance(rank)){
	  double delz=und->steplength();

	  // ----------------------------------------
	  // step 1 - apply most marker action  (always at beginning of a step)

	  bool sort=control->applyMarker(beam, field, und);


	  // ---------------------------------------
	  // step 2 - Advance electron beam

	  beam->track(delz,field,und);

	  // -----------------------------------------
	  // step 3 - Beam post processing, e.g. sorting


	  if (sort){
	    int shift=beam->sort();

	    if (shift!=0){
	      for (int i=0;i<field->size();i++){
		control->applySlippage(shift, field->at(i));  
	      }
	    }
	  }
  
	  // ---------------------------------------
	  // step 4 - Advance radiation field

	  for (int i=0; i<field->size();i++){
	    field->at(i)->track(delz,beam,und);
          }


	  //-----------------------------------------
	  // step 5 - Apply slippage

	  for (int i=0;i<field->size();i++){
	    control->applySlippage(und->slippage(), field->at(i));  
	  }

	  //-------------------------------
	  // step 6 - Calculate beam parameter stored into a buffer for output

	  beam->diagnostics(und->outstep(),und->getz());
	  for (int i=0;i<field->size();i++){
	    field->at(i)->diagnostics(und->outstep());
	  }

          
        }
     
        //---------------------------
        // end and clean-up 

	// perform last marker action

	bool sort=control->applyMarker(beam, field, und);
	if (sort){
	    int shift=beam->sort();

	    if (shift!=0){
	      for (int i=0;i<field->size();i++){
		control->applySlippage(shift, field->at(i));  
	      }
	    }
	}


	// write out diagnostic arrays

	if (rank==0){
	  cout << "Writing output file..." << endl;
	}

	control->output(beam,field,und);



	delete control;
      
        if (rank==0){
	  cout << endl << "Core Simulation done." << endl;
        }


        return 0;

}
::::::::::::::
src/Core/Incoherent.cpp
::::::::::::::
#include "Incoherent.h"
#include "Beam.h"

Incoherent::Incoherent(){
  sran=NULL;
  doLoss=false;
  doSpread=false;
}

Incoherent::~Incoherent(){}

void Incoherent::init(int base, int rank, bool doLoss_in,bool doSpread_in)
{

  doLoss=doLoss_in;
  doSpread=doSpread_in;


  RandomU rseed(base);
  double val;
  for (int i=0; i<=rank;i++){
    val=rseed.getElement();
  }
  val*=1e9;
  int locseed=static_cast<int> (round(val));
  if (sran !=NULL) { delete sran; }
  sran  = new RandomU (locseed);
  return;
}



void Incoherent::apply(Beam *beam, Undulator *und, double delz)
{  

  if (!und->inUndulator()) { return; }
  if ((!doLoss) && (!doSpread)) { return; }

  double gam0=und->getGammaRef();
  double awz=und->getaw();
  double xkw0=und->getku();


  double dgamsig=1.015e-27* xkw0 * xkw0 * awz * awz;

  if (und->isHelical()){
    dgamsig*= 1.42 *awz + 1./(1.+1.5*awz+0.95*awz*awz);
  } else {
    dgamsig*= 1.697*awz + 1./(1.+1.88*awz+0.8*awz*awz); 
  }

  if (!doSpread){ dgamsig=0;}

  dgamsig=sqrt(dgamsig*gam0*gam0*gam0*gam0*xkw0*delz)*sqrt(3.);


  double dgamavg=xkw0*gam0*awz;
  if(!doLoss) { dgamavg=0;}

  dgamavg=dgamavg*dgamavg*1.88e-15*delz;


  // apply energy change to electorn bunch
  int nbins=beam->nbins;
  if (beam->one4one){ nbins=1;}
  double dg=0;

  for (int islice=0;islice< beam->beam.size();islice++){
    int npart=beam->beam.at(islice).size();
    for (int ip=0; ip<npart; ip++){
      if ((ip % nbins) == 0){
         dg=-dgamavg+dgamsig*(2*sran->getElement()-1);
      }
      beam->beam.at(islice).at(ip).gamma+=dg;
    }
  }


  return;

}

::::::::::::::
src/Core/TrackBeam.cpp
::::::::::::::
#include "TrackBeam.h"
#include "Beam.h"

TrackBeam::TrackBeam(){}
TrackBeam::~TrackBeam(){}


void TrackBeam::track(double delz, Beam *beam,Undulator *und,bool lastStep=true)
{

  // get undulator parameter for the given step
  double aw,dax,day,ku,kx,ky;
  double qf,dqx,dqy;
  double cx,cy;
  double angle,lb,ld,lt;

  double gamma0=und->getGammaRef();
  und->getUndulatorParameters( &aw,&dax,&day,&ku,&kx,&ky);
  und->getQuadrupoleParameters(&qf,&dqx,&dqy);
  und->getCorrectorParameters(&cx,&cy);
  und->getChicaneParameters(&angle,&lb,&ld,&lt);

  double betpar0=sqrt(1-(1+aw*aw)/gamma0/gamma0);
  
  // effective focusing in x and y with the correct energy dependence
  double qquad=qf*gamma0;
  double qnatx=kx*aw*aw/gamma0/betpar0;  // kx has already the scaling with ku^2
  double qnaty=ky*aw*aw/gamma0/betpar0;  // same with ky

  double qx= qquad+qnatx;
  double qy=-qquad+qnaty;

  double xoff= qquad*dqx+qnatx*dax;
  double yoff=-qquad*dqx+qnatx*dax;
  // cout << "qnaty: " << qnaty/gamma0 << " Gamma0: " << gamma0 << endl;

  if (lastStep){
    if ((cx!=0) || (cy!=0)) { this->applyCorrector(beam,cx*gamma0,cy*gamma0); }
  } else {
    if (angle!=0) { this->applyChicane(beam,angle,lb,ld,lt,gamma0); }
  }
  // handle the different cases (drift, focusing and defocusing) with function pointers to member functions

  if (qx==0){ 
    this->ApplyX=&TrackBeam::applyDrift; 
  }else{
    xoff=xoff/qx;
    if (qx>0){
      this->ApplyX=&TrackBeam::applyFQuad;
    } else {
      this->ApplyX=&TrackBeam::applyDQuad;
    }
  }

  if (qy==0){ 
    this->ApplyY=&TrackBeam::applyDrift; 
  }else{
    yoff=yoff/qy;
    if (qy>0){
      this->ApplyY=&TrackBeam::applyFQuad;
    } else {
      this->ApplyY=&TrackBeam::applyDQuad;
    }
  }

  for (int i=0; i<beam->beam.size();i++){
    for (int j=0; j<beam->beam.at(i).size();j++){
      Particle *p=&beam->beam.at(i).at(j);
      double gammaz=sqrt(p->gamma*p->gamma-1- aw*aw - p->px*p->px - p->py*p->py); // = gamma*betaz=gamma*(1-(1+aw*aw)/gamma^2);
      (this->*ApplyX)(delz,qx,&(p->x),&(p->px),gammaz,xoff);
      (this->*ApplyY)(delz,qy,&(p->y),&(p->py),gammaz,yoff);
    }
  }



  return;
} 


void TrackBeam::applyDrift(double delz, double qf, double *x, double *px, double gammaz, double dx)
{
  *x+=(*px)*delz/gammaz;
  return;
}


void TrackBeam::applyFQuad(double delz, double qf, double *x, double *px, double gammaz, double dx)
{
  double foc=sqrt(qf/gammaz);
  double omg=foc*delz;
  double a1=cos(omg);
  double a2=sin(omg)/foc;
  double a3=-a2*foc*foc;
  double xtmp=*x-dx;
  *x =a1*xtmp+a2*(*px)/gammaz+dx;
  *px=a3*xtmp*gammaz+a1*(*px);
  return;
}


void TrackBeam::applyDQuad(double delz, double qf, double *x, double *px, double gammaz, double dx)
{
  double foc=sqrt(-qf/gammaz);
  double omg=foc*delz;
  double a1=cosh(omg);
  double a2=sinh(omg)/foc;
  double a3=a2*foc*foc;
  double xtmp=*x-dx;
  *x =a1*xtmp+a2*(*px)/gammaz+dx;
  *px=a3*xtmp*gammaz+a1*(*px);
  return;
}

void TrackBeam::applyCorrector(Beam *beam, double cx, double cy)
{ 

  for (int i=0; i<beam->beam.size();i++){
    for (int j=0; j<beam->beam.at(i).size();j++){
      beam->beam.at(i).at(j).px+=cx;
      beam->beam.at(i).at(j).py+=cy;
    }
  }
  return;
}


void TrackBeam::applyChicane(Beam *beam, double angle, double lb, double ld, double lt, double gamma0)
{ 
  // the tracking is done my applying the transfer matrix for the chicane and  backtracking for a drift over the length of the chicane
  // the effect of the R56 is applied here to the particle phase.  Then the normal tracking should do the momentum dependent change in the 
  // longitudinal position

  // the transfer matrix order is
  //  m -> bp -> ep -> d1 -> en -> bn -> d2 -> bn -> en-> d1 -> ep-> bp ->d3
  

  double m[4][4];
  double d1[4][4];
  double d2[4][4];
  double d3[4][4];
  double bp[4][4];
  double bn[4][4];
  double ep[4][4];
  double en[4][4];

  // construct the transfer matrix
  for (int i=0; i<4;i++){
    for (int j=0; j<4; j++){
      m[i][j]=0;
      d1[i][j]=0;
      d2[i][j]=0;
      d3[i][j]=0;
      bp[i][j]=0;
      bn[i][j]=0;
      ep[i][j]=0;
      en[i][j]=0;
    }
    m[i][i]=1;
    d1[i][i]=1;
    d2[i][i]=1;
    d3[i][i]=1;
    bp[i][i]=1;
    bn[i][i]=1;
    ep[i][i]=1;
    en[i][i]=1;
  }
  d1[0][1]=ld/cos(angle);  // drift between dipoles
  d1[2][3]=ld/cos(angle);   
  d2[0][1]=lt-4*lb-2*ld;   // drift in the middle
  d2[2][3]=lt-4*lb-2*ld;   
  d3[0][1]=-lt;            // negative drift over total chicane to get a zero element
  d3[2][3]=-lt;

  double R=lb/sin(angle);
  double Lpath=R*angle; 
  bp[2][3]=Lpath;  // positive deflection angle
  bp[0][0]=cos(angle);
  bp[0][1]=R*sin(angle);
  bp[1][0]=-sin(angle)/R;
  bp[1][1]=cos(angle);
 
  bn[2][3]=Lpath; // negative deflection angle
  bn[0][0]=cos(-angle);
  bn[0][1]=R*sin(-angle)*-1;
  bn[1][0]=-sin(-angle)/R*-1;
  bn[1][1]=cos(-angle);

  double efoc=tan(angle)/R;
  ep[1][0]=efoc;
  ep[3][2]=-efoc;
  en[1][0]=-efoc*-1;
  en[3][2]=efoc*-1;


  this->matmul(m,bp);
  this->matmul(m,ep);  
  this->matmul(m,d1);
  this->matmul(m,en);
  this->matmul(m,bn);
  this->matmul(m,d2);
  this->matmul(m,bn);
  this->matmul(m,en);
  this->matmul(m,d1);
  this->matmul(m,ep);
  this->matmul(m,bp);

  // transport matrix has been cross checked with Madx.  
  /*
  cout << "lt = " << lt << " angle = " << angle << " lb = " << lb << " ld = " << ld <<  endl;
  cout << m[0][0] << " " << m[0][1] << endl;
  cout << m[1][0] << " " << m[1][1] << endl;
  cout << m[2][2] << " " << m[2][3] << endl;
  cout << m[3][2] << " " << m[3][3] << endl;
  */

  this->matmul(m,d3);  // transport backwards because the main tracking still has to do the drift
  

  for (int i=0; i<beam->beam.size();i++){
    for (int j=0; j<beam->beam.at(i).size();j++){
      Particle *p=&beam->beam.at(i).at(j);
      double gammaz=sqrt(p->gamma*p->gamma-1- p->px*p->px - p->py*p->py); // = gamma*betaz=gamma*(1-(1+aw*aw)/gamma^2);

      double tmp=p->x;
      p->x =m[0][0]*tmp        +m[0][1]*p->px/gammaz;
      p->px=m[1][0]*tmp*gammaz +m[1][1]*p->px;
      tmp=p->y;
      p->y =m[2][2]*tmp        +m[2][3]*p->py/gammaz;
      p->py=m[3][2]*tmp*gammaz +m[3][3]*p->py;
      
    }
  }

  return;
}

void TrackBeam::matmul(double m[][4], double e[][4])
{
  double t[4][4];

  for (int i=0;i<4;i++){
    for (int j=0; j<4;j++){
      t[i][j]=0;
      for (int k=0;k<4;k++){
	t[i][j]+=e[i][k]*m[k][j];
      }
    }
  }
  for (int i=0;i<4;i++){
    for (int j=0; j<4;j++){
      m[i][j]=t[i][j];
    }
  }
}

void TrackBeam::applyR56(Beam *beam, Undulator *und, double lambda0)
{
  double angle,lb,ld,lt;

  double gamma0=und->getGammaRef();
  und->getChicaneParameters(&angle,&lb,&ld,&lt);
  if (angle==0) { return;}
  double R56=(4*lb/sin(angle)*(1-angle/tan(angle))+2*ld*tan(angle)/cos(angle))*angle;
  //    cout << "R56: " << R56 << endl;
  R56=R56*4*asin(1)/lambda0/gamma0;
  for (int i=0; i<beam->beam.size();i++){
    for (int j=0; j<beam->beam.at(i).size();j++){
      beam->beam.at(i).at(j).theta+=R56*(beam->beam.at(i).at(j).gamma-gamma0);
    }
  }
  return;

}
::::::::::::::
src/Core/Undulator.cpp
::::::::::::::
#include "Undulator.h"

Undulator::~Undulator(){}

Undulator::Undulator()
{
  istepz=-1;
  zstop=1e9;
}


void Undulator::updateOutput(double zstop_in,int nzout)
{


// calculate the size of the output record

  zstop=zstop_in;


  istepz=-1;
  nstepz=aw.size();
  nout=1;
  out.resize(nstepz+1);
  out[0]=true;  // first is always output
  for (int i=1;i<(nstepz+1);i++){   
    out[i]=false;
    if (((i % nzout)==0)&&(z[i-1]<zstop)){
      out[i]=true;
      nout++;
    } 
  }
  return;
}



void Undulator::updateMarker(int nfld, int npar, int nsort, double zstop)
{
  for (int i=0; i<marker.size();i++){
    if (nfld > 0){  // field dump
     if ((i % nfld) == 0) {
      marker[i]|=1;
     }
    }
    if (npar > 0){    // particle dump
     if ((i % npar) == 0) {
      marker[i]|=2;
     }
    }
    if (nsort > 0){    // sorting
     if ((i % nsort) == 0) {
      marker[i]|=4;
     }
    }
    if (z[i]>zstop) {  // stop calculation
      marker[i]|=8;
    }
  }  
  return;
} 


/*

bool Undulator::init(hid_t fid)
{

  nzout=1; 
  zstop=1e9;

  readDataDouble(fid,(char *)"/Global/gamma0",&gammaref,1);
  readDataDouble(fid,(char *)"/Global/zstop",&zstop,1);
  readDataInt(fid,(char *)"/Global/nzout",&nzout,1);
  

  int nwork=getDatasetSize(fid,(char *)"/Lattice/aw");

  aw.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/aw",&aw[0],nwork);
  helical.resize(nwork);
  readDataInt(fid,(char *)"/Lattice/helical",&helical[0],nwork);
  ax.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ax",&ax[0],nwork);
  ay.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ay",&ay[0],nwork);
  ku.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ku",&ku[0],nwork);
  kx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/kx",&kx[0],nwork);
  ky.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ky",&ky[0],nwork);
  gradx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/gradx",&gradx[0],nwork);
  grady.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/grady",&grady[0],nwork);
  qf.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/qf",&qf[0],nwork);
  qx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/qx",&qx[0],nwork);
  qy.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/qy",&qy[0],nwork);
  z.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/z",&z[0],nwork);
  dz.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/dz",&dz[0],nwork);
  slip.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/slippage",&slip[0],nwork);
  phaseshift.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/phaseshift",&phaseshift[0],nwork);
  cx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/cx",&cx[0],nwork);
  cy.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/cy",&cy[0],nwork);
  chic_lb.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/lb",&chic_lb[0],nwork);
  chic_ld.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ld",&chic_ld[0],nwork);
  chic_lt.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/lt",&chic_lt[0],nwork);
  chic_angle.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/delay",&chic_angle[0],nwork);
  marker.resize(nwork);
  readDataInt(fid,(char *)"/Lattice/marker",&marker[0],nwork);

  // convert delay into bending angle
  

  for (int i=0; i<aw.size(); i++){
    if (chic_angle[i]!=0){
       double delay=fabs(chic_angle[i]);
       double tmin=0;
       double tmax=asin(1)-0.001;
       bool converged=false;
       double theta,d;
       while (!converged){
         theta=0.5*(tmax+tmin);
         d=4*chic_lb[i]*(theta/sin(theta)-1)+2*chic_ld[i]*(1/cos(theta)-1);
         if (d>delay) {
           tmax=theta;
         } else {
          tmin=theta;
	 }
         if (fabs(delay-d)<1e-15) { converged=true; }
       }
       chic_angle[i]=theta;
    }
  }
 






// calculate the size of the output record

  istepz=-1;
  nstepz=aw.size();
  nout=1;
  out.resize(nstepz+1);
  out[0]=true;  // first is always output
  for (int i=1;i<(nstepz+1);i++){   
    out[i]=false;
    if (((i % nzout)==0)&&(z[i]<zstop)){
      out[i]=true;
      nout++;
    } 
  }

  return true;

}
*/


bool Undulator::advance(int rank)
{
  istepz++;
  
  if (istepz >= nstepz){  // end reached?
    return false;
  }

  if ((marker[istepz]&8)>0){ // check for 3rd bit set in marker value for stoping calculation
    if (rank==0){
      cout << "Calculation terminated due to requested stop. Missing output padded with zeros in output file" << endl;
    }
    return false; 
  }

  int dstepz=nstepz/10;
  if (dstepz<1){dstepz=1;}

  if (((istepz % dstepz) == 0) && (istepz >0) && (rank==0)){
    cout << "  Calculation: " <<10*istepz/dstepz << "% done" << endl;
  }
  return true;
}

double Undulator::fc(int h)
{
  BesselJ bessj;

  double coup=aw[istepz]; 
  if (this->isHelical()){
    if (h==1){
      return coup;
    } else {
      return 0;
    }
  } else {
    double xi = aw[istepz]*aw[istepz];
    xi=0.5 * xi /(1+xi)*static_cast<double>(h);
    int h0,h1;
    if ((h % 2) == 1){
      h0=(h-1)/2;
      h1 = h0+1;
      return coup*(bessj.value(h0,xi)-bessj.value(h1,xi))*pow(-1.,h0);
    } else {
      h0=(h-2)/2;
      h1 = h0+2;
      return coup*0.5*(bessj.value(h0,xi)-bessj.value(h1,xi))*pow(-1.,h0);
    }
  }
}


double Undulator::faw2(double x, double y){  // square of the transverse dependence of the undulator field.
  double dx=x-ax[istepz];
  double dy=y-ay[istepz]; 
  return (1+kx[istepz]*dx*dx+ky[istepz]*dy*dy+2*(gradx[istepz]*dx+grady[istepz]*dy)); // note kx is scaled as XKX*ku*ku in Lattice.cpp, gradx as ku*GRADX.
}



double Undulator::faw(double x, double y){  // transverse dependence of the undulator field.
  double dx=x-ax[istepz];
  double dy=y-ay[istepz]; 
  return (1+0.5*(kx[istepz]*dx*dx+ky[istepz]*dy*dy)+gradx[istepz]*dx+grady[istepz]*dy); // note kx is scaled as XKX*ku*ku in Lattice.cpp, gradx as ku*GRADX.
}
::::::::::::::
src/IO/HDF5base.cpp
::::::::::::::
#include "HDF5base.h"

HDF5Base::HDF5Base(){}
HDF5Base::~HDF5Base(){}


//----------------------------
// routines taken from Output and writeHDF5



void HDF5Base::writeBuffer(hid_t gid, string dataset,vector<double> *data){


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








//-------------------------------
// write of a single node in parallal HDF5




void HDF5Base::writeDouble1DExist(hsize_t datsize, double *data, hid_t gid, string dataset)
{

  int dataset_rank=1;
  hsize_t dims[1] = {datsize}; // dataset size must be the same for all nodes

  hid_t did=H5Dopen(gid,dataset.c_str(),H5P_DEFAULT);

  hid_t filespace=H5Dget_space(did);
  hid_t memspace=H5Screate_simple(dataset_rank,dims,NULL);


  // set up transfer and write
  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pid,H5FD_MPIO_INDEPENDENT);    
  H5Dwrite(did,H5T_NATIVE_DOUBLE,memspace,filespace,pid,data);
 
  // close all HDF5 stuff except for the file id fid
  H5Dclose(did);
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Pclose(pid);


  return;
}





void HDF5Base::writeChar1D(hsize_t reclen, hsize_t datsize, const char *data, hid_t gid, string name)
{

  int dataset_rank=1;
  hsize_t dims[1] = {datsize}; // dataset size must be the same for all nodes
  hid_t filespace = H5Screate_simple(dataset_rank,dims,NULL);

  hid_t did=H5Dcreate(gid,name.c_str(),H5T_NATIVE_CHAR,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);  
  H5Sclose(filespace);

 

  hsize_t count[1] = {reclen};     // length of record entry 
  hsize_t offset[1] = {0};   // offset of record entry
  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pid,H5FD_MPIO_COLLECTIVE);    

  hid_t memspace=H5Screate_simple(dataset_rank,count,NULL);
  if (reclen==0) {H5Sselect_none(memspace);}

  
  filespace= H5Dget_space(did);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  if (reclen==0){ H5Sselect_none(filespace);  }

  if (reclen==0) {
      H5Dwrite(did,H5T_NATIVE_CHAR,memspace,filespace,pid,NULL);
  } else { 
      H5Dwrite(did,H5T_NATIVE_CHAR,memspace,filespace,pid,data);
  }
  

  H5Dclose(did);
  H5Pclose(pid);
  H5Sclose(memspace);
  H5Sclose(filespace);

  return;
}


void HDF5Base::writeDouble1D(hsize_t reclen, hsize_t datsize, double *data, hid_t gid, string name)
{

  int dataset_rank=1;
  hsize_t dims[1] = {datsize}; // dataset size must be the same for all nodes
  hid_t filespace = H5Screate_simple(dataset_rank,dims,NULL);

  hid_t did=H5Dcreate(gid,name.c_str(),H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);  
  H5Sclose(filespace);

 

  hsize_t count[1] = {reclen};     // length of record entry 
  hsize_t offset[1] = {0};   // offset of record entry
  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pid,H5FD_MPIO_COLLECTIVE);    

  hid_t memspace=H5Screate_simple(dataset_rank,count,NULL);
  if (reclen==0) {H5Sselect_none(memspace);}

  
  filespace= H5Dget_space(did);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  if (reclen==0){ H5Sselect_none(filespace);  }

  if (reclen==0) {
      H5Dwrite(did,H5T_NATIVE_DOUBLE,memspace,filespace,pid,NULL);
  } else { 
      H5Dwrite(did,H5T_NATIVE_DOUBLE,memspace,filespace,pid,data);
  }
  

  H5Dclose(did);
  H5Pclose(pid);
  H5Sclose(memspace);
  H5Sclose(filespace);

  return;
}

void HDF5Base::writeInt1D(hsize_t reclen, hsize_t  datsize, int *data, hid_t gid, string name)
{

  int dataset_rank=1;
  hsize_t dims[1] = {datsize}; // dataset size must be the same for all nodes
  hid_t filespace = H5Screate_simple(dataset_rank,dims,NULL);

  hid_t did=H5Dcreate(gid,name.c_str(),H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);  
  H5Sclose(filespace);

 

  hsize_t count[1] = {reclen};     // length of record entry 
  hsize_t offset[1] = {0};   // offset of record entry
  hid_t pid =  H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(pid,H5FD_MPIO_COLLECTIVE);    

  hid_t memspace=H5Screate_simple(dataset_rank,count,NULL);
  if (reclen==0) {H5Sselect_none(memspace);}

  
  filespace= H5Dget_space(did);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  if (reclen==0){ H5Sselect_none(filespace);  }

  if (reclen==0) {
      H5Dwrite(did,H5T_NATIVE_INT,memspace,filespace,pid,NULL);
  } else { 
      H5Dwrite(did,H5T_NATIVE_INT,memspace,filespace,pid,data);
  }
  

  H5Dclose(did);
  H5Pclose(pid);
  H5Sclose(memspace);
  H5Sclose(filespace);

  return;
}



//----------------------
// generating expandable dataset

void HDF5Base::createExpDataset(hid_t fid, char *name, hsize_t nz, hsize_t ns)
{
  hsize_t dims[2] = {nz, ns};
  hsize_t maxdims[2]={H5S_UNLIMITED, H5S_UNLIMITED};
  hid_t dataspace = H5Screate_simple(2,dims,maxdims);
  hid_t cparms=H5Pcreate(H5P_DATASET_CREATE);
  hsize_t chunk_dims[2]={10,10};
  herr_t status=H5Pset_chunk( cparms, 2, chunk_dims);
  hid_t dataset=H5Dcreate(fid,name,H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT,cparms,H5P_DEFAULT);

  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  return;

}

void HDF5Base::expandDataset(hid_t fid, vector<double> *rec, int pos, hsize_t recsize, hsize_t slice, char *name)
{

  double *data= new double[recsize];
  
  int stride=rec->size()/recsize;
  for(int i=0;i<recsize;i++){
    data[i]=rec->at(i*stride+pos);
    //    cout << data[i] << endl;
  }
  //  return;

  hid_t did=H5Dopen(fid,name,H5P_DEFAULT);
  hsize_t size[2] = {recsize, slice};  
  herr_t status=H5Dset_extent(did,size); 

  hid_t filespace=H5Dget_space(did);  
  hsize_t offset[2] = {0,slice-1};
  hsize_t dims[2]={recsize,1};
  hid_t dataspace=H5Screate_simple(2,dims,dims);  

  status=H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,dims,NULL);
  status=H5Dwrite(did,H5T_NATIVE_DOUBLE,dataspace,filespace,H5P_DEFAULT,data);
  
  status=H5Sclose(dataspace);
  status=H5Sclose(filespace);
  status=H5Dclose(did);

  delete [] data;
  return;
}


//------------------------
// writing procedures

void HDF5Base::writeDataDouble(hid_t fid, const char *name, const double *data, int size)
{
  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dcreate(fid,name,H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

void HDF5Base::writeDataChar(hid_t fid, const char *name, const char *data, int size)
{
  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dcreate(fid,name,H5T_NATIVE_CHAR,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Dwrite(dataset_id,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}



void HDF5Base::writeDataInt(hid_t fid, const char *name, const int *data, int size)
{
  hsize_t dims[1];
  dims[0]=size;
  hid_t dataspace_id=H5Screate_simple(1,dims,NULL);
  hid_t dataset_id=H5Dcreate(fid,name,H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Dwrite(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
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

 
::::::::::::::
src/IO/HDF5Beam.cpp
::::::::::::::
#include "SDDSBeam.h"
#include "MPEProfiling.h"

SDDSBeam::SDDSBeam()
{
  file="";
  ds=0.01;
  output=false; 
  xcen=0;
  ycen=0;
  pxcen=0;
  pycen=0;
  gamma=12000;
  betax=15;
  alphax=0;
  betay=15;
  alphay=0;
  charge=0;
  center=false;
  match=false;
  matchs0=0;
  matchs1=1;
  align=0;
  aligns0=0;
  aligns1=0;
  
}

SDDSBeam::~SDDSBeam(){}

void SDDSBeam::usage(){

  cout << "List of keywords for sddsbeam" << endl;
  cout << "&sddsbeam" << endl;
  cout << " string file = <empty> " << endl;
  cout << " double charge   = 0" << endl;
  cout << " double slicewidth = 0.01" << endl;
  cout << " bool output = false " << endl;
  cout << " bool center = false " << endl;
  cout << " double gamma0 = gammaref " << endl;
  cout << " double x0 = 0 " << endl;
  cout << " double y0 = 0 " << endl;
  cout << " double px0 = 0 " << endl;
  cout << " double py0 = 0 " << endl;
  cout << " bool match = false " << endl;
  cout << " double betax  = 15 / matched" << endl;
  cout << " double alphax  = 0 / matched" << endl;
  cout << " double betay  = 15 / matched" << endl;
  cout << " double alphay  = 0 / matched" << endl;
  cout << " double match_start = 0 " << endl;
  cout << " double match_end = 1 " << endl;
  cout << " int align = 0 " << endl;
  cout << " double align_start = 0 " << endl;
  cout << " double align_end = 1 " << endl;
  cout << "&end" << endl << endl;
  return;
}

bool SDDSBeam::init(int inrank, int insize, map<string,string> *arg, Beam *beam, Setup *setup, Time *time, Lattice *lat)
{

  rank=inrank;
  size=insize;

  gamma=setup->getReferenceEnergy();           // get default energy from setup input deck
  lat->getMatchedOptics(&betax,&alphax,&betay,&alphay);  // use matched value if calculated

  double lambda=setup->getReferenceLength();   // reference length for theta
  double sample=static_cast<double>(time->getSampleRate());         // check slice length

  bool one4one=setup->getOne4One();
  bool shotnoise=setup->getShotNoise();
  int npart=setup->getNpart();
  int nbins=setup->getNbins();

  double theta0=4.*asin(1.);
  if (one4one) {
      nbins=1;
      theta0*=sample;
  }
  if ( (npart % nbins) != 0){
    if (rank==0) { cout << "*** Error: NPART is not a multiple of NBINS" << endl; } 
    return false;
  }
  
  theta0/=static_cast<double>(nbins);
 
  map<string,string>::iterator end=arg->end();

  if (arg->find("file")!=end)       {file   = arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("charge")!=end)     {charge = atof(arg->at("charge").c_str());     arg->erase(arg->find("charge"));}
  if (arg->find("slicewidth")!=end) {ds     = atof(arg->at("slicewidth").c_str()); arg->erase(arg->find("slicewidth"));}
  if (arg->find("match_start")!=end)    {matchs0    = atof(arg->at("match_start").c_str());    arg->erase(arg->find("match_start"));}
  if (arg->find("match_end")!=end)      {matchs1    = atof(arg->at("match_end").c_str());      arg->erase(arg->find("match_end"));}
  if (arg->find("align_start")!=end)    {aligns0    = atof(arg->at("align_start").c_str());    arg->erase(arg->find("align_start"));}
  if (arg->find("align_end")!=end)      {aligns1    = atof(arg->at("align_end").c_str());      arg->erase(arg->find("align_end"));}
  if (arg->find("betax")!=end)    {betax = atof(arg->at("betax").c_str()); arg->erase(arg->find("betax"));}
  if (arg->find("betay")!=end)    {betay = atof(arg->at("betay").c_str()); arg->erase(arg->find("betay"));}
  if (arg->find("alphax")!=end)   {alphax= atof(arg->at("alphax").c_str());arg->erase(arg->find("alphax"));}
  if (arg->find("alphay")!=end)   {alphay= atof(arg->at("alphay").c_str());arg->erase(arg->find("alphay"));}
  if (arg->find("x0")!=end)       {xcen  = atof(arg->at("x0").c_str());    arg->erase(arg->find("x0"));}
  if (arg->find("y0")!=end)       {ycen  = atof(arg->at("y0").c_str());    arg->erase(arg->find("y0"));}
  if (arg->find("px0")!=end)      {pxcen = atof(arg->at("px0").c_str());   arg->erase(arg->find("px0"));}
  if (arg->find("py0")!=end)      {pycen = atof(arg->at("py0").c_str());   arg->erase(arg->find("py0"));}
  if (arg->find("gamma0")!=end)   {gamma = atof(arg->at("gamma0").c_str());arg->erase(arg->find("gamma0"));}
  if (arg->find("align")!=end)    {align = atoi(arg->at("align").c_str()); arg->erase(arg->find("align"));}
  if (arg->find("match")!=end)    {match = atob(arg->at("match").c_str()); arg->erase(arg->find("match"));}
  if (arg->find("center")!=end)   {center= atob(arg->at("center").c_str());arg->erase(arg->find("center"));}
  if (arg->find("output")!=end)   {output= atob(arg->at("output").c_str());arg->erase(arg->find("output"));}



  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &sddsbeam" << endl; this->usage();}
    return false;
  }

  mpe.logCalc(false,true,"System call SDDS2HDF5");
  if (rank==0) { cout << "Converting SDDS distribution file " << file << " into HDF5 file... " << endl;}
  string command="sdds2hdf-dist.sh " + file;
  int status=0;
  if (rank==0){
    status=system(command.c_str());
    if (status!=0){ cout << "*** Error: SDDS Distribution file " << file << " does not exist" << endl;} 
  }

  mpe.logCalc(true,true,"System call SDDS2HDF5");


  MPI::COMM_WORLD.Bcast(&status,1,MPI::INT,0); // this statement keeps also all nodes on hold till the converison has been done
  if (status!=0) {return false;}



  if (rank==0) { cout << "Importing converted distribution file... " << endl;}


  mpe.logIO(false,true,"Read External Distribution"); 
  string anafile=file;
  string h5file=file.append(".h5");
  

  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  hid_t fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,pid);
  H5Pclose(pid);

  string dset="t";
  int ntotal=getDatasetSize(fid, (char *)dset.c_str());
  double dQ=charge/static_cast<double> (ntotal);

  int nchunk=ntotal/size;
  if ((ntotal % size) !=0) {nchunk++;}

  int nsize=nchunk;
  if ((rank*nchunk+nsize)>ntotal) { nsize=ntotal-rank*nchunk; }

  t.resize(nsize);
  g.resize(nsize);
  x.resize(nsize);
  y.resize(nsize);
  px.resize(nsize);
  py.resize(nsize);

  string dname="t";
  readDouble1D(fid,dname.c_str(),&t[0],nsize,rank*nchunk);
  dname="p"; 
  readDouble1D(fid,dname.c_str(),&g[0],nsize,rank*nchunk);
  dname="x"; 
  readDouble1D(fid,dname.c_str(),&x[0],nsize,rank*nchunk);
  dname="xp"; 
  readDouble1D(fid,dname.c_str(),&px[0],nsize,rank*nchunk);
  dname="y"; 
  readDouble1D(fid,dname.c_str(),&y[0],nsize,rank*nchunk);
  dname="yp"; 
  readDouble1D(fid,dname.c_str(),&py[0],nsize,rank*nchunk);

  H5Fclose(fid);

  mpe.logIO(true,true,"Read External Distribution"); 


  if (rank==0) { cout << "Analysing external distribution... " << endl;}


  for (int i=0; i<nsize; i++){
    t[i]*=-3e8;
    g[i]+=1.;
  }


  double tmin,tmax;

  double tmp=*min_element(t.begin(),t.end());
  MPI::COMM_WORLD.Allreduce(&tmp,&tmin,1,MPI::DOUBLE,MPI::MIN);
  tmp=*max_element(t.begin(),t.end());
  MPI::COMM_WORLD.Allreduce(&tmp,&tmax,1,MPI::DOUBLE,MPI::MAX);
  tmp=accumulate(g.begin(),g.end(),0);


  double ttotal=tmax-tmin;

  for (int i=0; i<nsize; i++){
    t[i]-=tmin;
  }


  if (rank==0) {
    cout << "Analysis of the imported distribution" << endl;  
    cout << "   Total Bunch Length  (microns): " << ttotal*1e6 << endl;
  }

  this->analyse(ttotal,nsize);

  if (center) {
    if (rank==0){cout << "Centering external distribution..." << endl; }
    double ratio=sqrt(gavg/gamma);
    gavg=gamma-gavg;
    xavg=xcen-xavg;
    yavg=ycen-yavg;
    pxavg=pxcen-pxavg;
    pyavg=pycen-pyavg;
    for (int i=0; i<nsize; i++){
      g[i]+=gavg;
      x[i]+=xavg;
      y[i]+=yavg;
      px[i]+=pxavg;
      py[i]+=pyavg;
      x[i]*=ratio;   // rescaling is needed to preserve emittance
      y[i]*=ratio;
      px[i]*=ratio;
      py[i]*=ratio;
    }
  }

  if (match) {
    if (rank==0){cout << "Matching external distribution..." << endl; }
    for (int i=0; i<nsize; i++){
      px[i]+=(ax/bx)*x[i];
      py[i]+=(ay/by)*y[i];
      x[i]*=sqrt(betax/bx);
      y[i]*=sqrt(betay/by);
      px[i]*=sqrt(bx/betax);
      py[i]*=sqrt(by/betay);
      px[i]-=(alphax/betax)*x[i];
      py[i]-=(alphay/betay)*y[i];
    }
  } 

  if ((match)||(center)){
    if (rank==0) {cout << "Reanalysing matched and aligned distribution..." << endl; }
    this->analyse(ttotal,nsize);
  }


  // slicing external distribution to fill given slizes

  vector<double> s;
  int nslice=time->getPosition(&s);
  int node_off=time->getNodeOffset();
  int node_len=time->getNodeNSlice();
  beam->beam.resize(node_len);

  double smin=s[node_off];
  double smax=s[node_off+node_len-1];
  
  if (rank==0) {cout << "Sorting external distribution..." << endl; }

  double dslen=ds*ttotal;  // ds is the relative width to extract the samples (equivalent to 1/NDCUT)

  vector<vector<Particle> > dist;
  dist.resize(1);

  // copying all particles into the dist vector to enable sorting
  Particle part;
  for (int i=0; i<nsize; i++){
      part.theta=t[i];
      part.gamma=g[i];
      part.x=x[i];
      part.y=y[i];
      part.px=px[i]*g[i];
      part.py=py[i]*g[i];
      dist[0].push_back(part);
  }
  t.clear();
  g.clear();
  x.clear();
  y.clear();
  px.clear();
  py.clear();


  Sorting sort;
  sort.init(rank,size,smin,smax,dslen,6,false,true); 
  sort.globalSort(&dist,0);  

  // now each node has all the particles, which is needed for the phase space reconstruction
  if (rank==0) {cout << "Generarting internal particle distribution..." << endl; }
  mpe.logLoading(false,"External Distribution");

  beam->current.resize(node_len);
  beam->beam.resize(node_len);

  this->initRandomSeq(setup->getSeed());
  ShotNoise sn;
  sn.init(setup->getSeed(),rank);

  int nwork=100;
  Particle *work;
  work=new Particle [nwork];
  for (int islice=0; islice<node_len;islice++){

    // step 1 - select all particles needed for the reconstruction
    double sloc=s[islice+node_off];
    for (int i=0; i<dist[0].size();i++){
      if ((dist[0].at(i).theta>(sloc-0.5*dslen))&&(dist[0].at(i).theta<(sloc+0.5*dslen))){
	beam->beam.at(islice).push_back(dist[0].at(i));
      }
    }
   
    // step 2 - calculate the current and number of particles.
    int ncount = beam->beam.at(islice).size();
    int mpart;
    beam->current[islice]=static_cast<double>(ncount)*dQ*3e8/dslen;
    if (one4one){
      npart=static_cast<int>(round(beam->current[islice]*lambda*sample/ce));
      mpart=npart;
    } else {
      mpart=npart/nbins;
    }

    // step 3 - bring initial distribution to the right size

    if (beam->beam.at(islice).size() >= mpart){
      this->removeParticles(&beam->beam.at(islice),mpart);
    } else {
      this->addParticles(&beam->beam.at(islice),mpart);
    }

    // step 4 - refill particle phase completely new
    for (int i=0;i<beam->beam.at(islice).size();i++){
      beam->beam.at(islice).at(i).theta=theta0*ran->getElement();
    }

    if (!one4one){
      mpart=beam->beam.at(islice).size();
      beam->beam.at(islice).resize(mpart*nbins);
      for (int i=mpart; i>0; i--){
        int i1=i-1;
        int i2=nbins*i1;
        for (int j=0;j<nbins;j++){
          beam->beam.at(islice).at(i2+j).gamma=beam->beam.at(islice).at(i1).gamma;
  	  beam->beam.at(islice).at(i2+j).x    =beam->beam.at(islice).at(i1).x;
          beam->beam.at(islice).at(i2+j).y    =beam->beam.at(islice).at(i1).y;
          beam->beam.at(islice).at(i2+j).px   =beam->beam.at(islice).at(i1).px;
          beam->beam.at(islice).at(i2+j).py   =beam->beam.at(islice).at(i1).py;
          beam->beam.at(islice).at(i2+j).theta=beam->beam.at(islice).at(i1).theta+j*theta0;     
	}
      }
      double ne=round(beam->current[islice]*lambda*sample/ce);
      if (mpart*nbins>nwork){
      	nwork=mpart*nbins;
        delete[] work;
        work=new Particle [nwork];
      }
      for (int i=0;i<mpart;i++){
      	work[i].theta=beam->beam.at(islice).at(i).theta;  
      }
      sn.applyShotNoise(work,mpart*nbins,nbins,ne); /////////////////////////////////////////////////////////////////////// shouldn-t it be copy back
      for (int i=0;i<mpart;i++){
      	beam->beam.at(islice).at(i).theta=work[i].theta;  
      }

    }
   
  }

  mpe.logLoading(true,"External Distribution");

  dist[0].clear();
  delete ran;
  delete [] work;
  
  if (output){
    mpe.logIO(false,true,"Write Distribution Analysis"); 
    if (rank==0) {cout << "Writing analysis to file..." << endl; }

    Output *out=new Output;
    string outfile=anafile.append(".ana.h5");

    hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
    hid_t fid=H5Fcreate(outfile.c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
    H5Pclose(pid);
    H5Fclose(fid);

    out->open(outfile.c_str(),size,NULL,nslice,1,node_off,node_len);
    out->writeBeam(0,beam);
    out->writeGlobalBeam(beam,bx,by,ax,ay);
    out->close();
    delete out;
    mpe.logIO(true,true,"Write Distribution Analysis"); 

  }


  mpe.logEvent("End: SDDSBeam::init");


  return true;


}


void SDDSBeam::addParticles(vector<Particle> *beam, int mpart){
   
  int  ndist=beam->size();
  if (ndist==0){
    return;
  }

  // step 1 - calculate the center and the rms beam size
  double g1,x1,y1,px1,py1,g2,x2,y2,px2,py2;
  g1=0;
  x1=0;
  y1=0;
  px1=0;
  py1=0;
  g2=0;
  x2=0;
  y2=0;
  px2=0;
  py2=0;
  for (int i=0; i<ndist; i++){
    g1 +=beam->at(i).gamma;
    g2 +=beam->at(i).gamma * beam->at(i).gamma;
    x1 +=beam->at(i).x;
    x2 +=beam->at(i).x * beam->at(i).x;
    px1+=beam->at(i).px;
    px2+=beam->at(i).px * beam->at(i).px;
    y1 +=beam->at(i).y;
    y2 +=beam->at(i).y * beam->at(i).y;
    py1+=beam->at(i).py;
    py2+=beam->at(i).py * beam->at(i).py;
  }
  double scl=1/static_cast<double>(ndist);

  g1*=scl;
  g2=sqrt(fabs(g2*scl-g1*g1));
  x1*=scl;
  x2=sqrt(fabs(x2*scl-x1*x1));
  px1*=scl;
  px2=sqrt(fabs(px2*scl-px1*px1));
  y1*=scl;
  y2=sqrt(fabs(y2*scl-y1*y1));
  py1*=scl;
  py2=sqrt(fabs(py2*scl-py1*py1));

  // step 2 - invert the beam size for normalization and check for "cold" dimensions, e.g. zero energy spread
  if (g2==0) { g2=1; } else { g2=1/g2; }
  if (x2==0) { x2=1; } else { x2=1/x2; }
  if (y2==0) { y2=1; } else { y2=1/y2; }
  if (px2==0) { px2=1; } else { px2=1/px2; }
  if (py2==0) { py2=1; } else { py2=1/py2; }

  // step 3 - normalize distribution so that it is aligned to the origin and has an rms size of unity in all dimensions
  for (int i=0; i<ndist; i++){
    beam->at(i).gamma=(beam->at(i).gamma - g1)*g2;
    beam->at(i).x    =(beam->at(i).x     - x1)*x2;
    beam->at(i).y    =(beam->at(i).y     - y1)*y2;
    beam->at(i).px   =(beam->at(i).px    - px1)*px2;
    beam->at(i).py   =(beam->at(i).py    - py1)*py2;
  }
  
  // step 4 - add particles
  int ndist0=ndist;
  Particle par;
  while (ndist<mpart){
    int n1=static_cast<int>(floor(static_cast<double>(ndist0)*ran->getElement()));
    double rmin=1e9;
    int n2=n1;
    for (int i=0; i<ndist0;i++){
       double r=this->distance(beam->at(n1),beam->at(i)); 
       if ((r<rmin) && ( i!=n1 )) {
           n2=i;
           rmin=r;
       }
    }
    par.gamma=0.5*(beam->at(n1).gamma+beam->at(n2).gamma)+(2*ran->getElement()-1)*(beam->at(n1).gamma-beam->at(n2).gamma);
    par.x =0.5*(beam->at(n1).x +beam->at(n2).x) +(2*ran->getElement()-1)*(beam->at(n1).x -beam->at(n2).x);
    par.px=0.5*(beam->at(n1).px+beam->at(n2).px)+(2*ran->getElement()-1)*(beam->at(n1).px-beam->at(n2).px);
    par.y =0.5*(beam->at(n1).y +beam->at(n2).y) +(2*ran->getElement()-1)*(beam->at(n1).y -beam->at(n2).y);
    par.py=0.5*(beam->at(n1).py+beam->at(n2).py)+(2*ran->getElement()-1)*(beam->at(n1).py-beam->at(n2).py);
    beam->push_back(par);    
    ndist++;
  }
 

  // step 5 - scale back

  for (int i=0; i<beam->size(); i++){
    beam->at(i).gamma=beam->at(i).gamma/g2 + g1;
    beam->at(i).x    =beam->at(i).x/x2 + x1;
    beam->at(i).y    =beam->at(i).y/y2 + y1;
    beam->at(i).px   =beam->at(i).px/px2 + px1; 
    beam->at(i).py   =beam->at(i).py/py2 + py1;
  }


  return;
}


double  SDDSBeam::distance(Particle p1, Particle p2){

  double tmp=p1.gamma-p2.gamma;  
  double r=tmp*tmp*ran->getElement();
  tmp=p1.x-p2.x;  
  r+=tmp*tmp*ran->getElement();
  tmp=p1.y-p2.y;  
  r+=tmp*tmp*ran->getElement();
  tmp=p1.px-p2.px;  
  r+=tmp*tmp*ran->getElement();
  tmp=p1.py-p2.py;  
  r+=tmp*tmp*ran->getElement();
  return r;
}


void SDDSBeam::removeParticles(vector<Particle> *beam,int mpart)
{
  int ndist=beam->size();
  while(ndist>mpart){
    int idx=static_cast<int>(floor(static_cast<double>(ndist)*ran->getElement()));
    beam->at(idx)=beam->at(ndist-1);
    ndist--;
  }
  beam->resize(mpart);
  return;
}

void SDDSBeam::analyse(double ttotal,int nsize)
{

  int ncount = 0;
  int nmean;
  double mt0=matchs0*ttotal;
  double mt1=matchs1*ttotal;

  double a2,a1,b2,b1,ab; 
  double c2,c1,d2,d1,cd; 
  double e1=0;
  a2=0;a1=0;b2=0;b1=0;ab=0;
  c2=0;c1=0;d2=0;d1=0;cd=0;
  
  //  cout << "Rank: " << rank << " Size: " << nsize << endl;

  for (int i=0; i <nsize;i++){
    if ((t[i]>mt0)&&(t[i]<mt1)){
      ncount++;
      a1+=x[i];
      a2+=x[i]*x[i];
      b1+=px[i];
      b2+=px[i]*px[i];
      ab+=x[i]*px[i];
      c1+=y[i];
      c2+=y[i]*y[i];
      d1+=py[i];
      d2+=py[i]*py[i];
      cd+=y[i]*py[i];
      e1+=g[i];
    }
  }


  MPI::COMM_WORLD.Allreduce(&ncount,&nmean,1,MPI::INT,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&e1,&gavg, 1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&a1,&xavg, 1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&b1,&pxavg,1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&c1,&yavg, 1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&d1,&pyavg,1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&a2,&xvar, 1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&b2,&pxvar,1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&ab,&xpx,  1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&c2,&yvar, 1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&d2,&pyvar,1,MPI::DOUBLE,MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cd,&ypy,  1,MPI::DOUBLE,MPI::SUM);

  if (nmean>0){
    double tmp=1./static_cast<double>(nmean);
    gavg*=tmp;
    xavg*=tmp;
    yavg*=tmp;
    pxavg*=tmp;
    pyavg*=tmp;
    xvar*=tmp;
    yvar*=tmp;
    pxvar*=tmp;
    pyvar*=tmp;
    xpx*=tmp;
    ypy*=tmp;
  }

  ex=sqrt(fabs((xvar-xavg*xavg)*(pxvar-pxavg*pxavg)-(xpx-xavg*pxavg)*(xpx-xavg*pxavg)))*gavg;
  ey=sqrt(fabs((yvar-yavg*yavg)*(pyvar-pyavg*pyavg)-(ypy-yavg*pyavg)*(ypy-yavg*pyavg)))*gavg;
  bx=(xvar-xavg*xavg)/ex*gavg;
  by=(yvar-yavg*yavg)/ey*gavg;
  ax=-(xpx-xavg*pxavg)*gavg/ex;
  ay=-(ypy-yavg*pyavg)*gavg/ey;

  if (rank==0) {
       cout << "   Length for Matching (microns): " << (mt1-mt0)*1e6 << endl; 
       cout << "   Energy                  (MeV): " << gavg*eev*1e-6 << endl;
       cout << "   Norm. Emittance in x (micron): " << ex*1e6 << endl;
       cout << "   Norm. Emittance in y (micron): " << ey*1e6 << endl;
       cout << "   Beta Function in x        (m): " << bx << endl;
       cout << "   Beta Function in y        (m): " << by << endl;
       cout << "   Alpha Function in x          : " << ax << endl;
       cout << "   Alpha Function in y          : " << ay << endl;
       cout << "   Beam center in x     (micron): " << xavg*1e6 << endl;
       cout << "   Beam center in y     (micron): " << yavg*1e6 << endl;
       cout << "   Beam center in px            : " << pxavg << endl;
       cout << "   Beam center in py            : " << pyavg << endl;
  }
  return;


}

void SDDSBeam::initRandomSeq(int base)
{
     RandomU rseed(base);
     double val;
     for (int i=0; i<=rank+10000;i++){
        val=rseed.getElement();
     }
     val*=1e9;
     int locseed=static_cast<int> (round(val));
     ran  =  new RandomU (locseed);
     return;
}
::::::::::::::
src/IO/Output.cpp
::::::::::::::
#include "Output.h"
#include <stdlib.h> 

// Meta info
#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#include <fstream>
#include <streambuf>

extern bool MPISingle;

Output::Output(){
  fid=-1;
}

Output::~Output(){}

void Output::close(){

  // output of unclosed elements

  
   int norphans = H5Fget_obj_count(fid, H5F_OBJ_ALL);
   if (norphans > 1) { /* expect 1 for the file we have not closed */
       int i;
       H5O_info_t info;
       char name[64];
       hid_t *objects = (hid_t *) calloc(norphans, sizeof(hid_t));
       H5Fget_obj_ids(fid, H5F_OBJ_ALL, -1, objects);
       for (i=0; i<norphans; i++) {
           H5Oget_info(objects[i], &info);
           H5Iget_name(objects[i], name, 64);
	   cout << "Slice " <<s0/ds << " : " << i+1 << " of " << norphans << " things still open: " << objects[i] << " with name " << name << " of type " << info.type << endl; 
	   //           printf("%d of %zd things still open: %d with name %s of type %d", i, norphans, objects[i], name, info.type);
       }
   }



  H5Fclose(fid);

  
}


void Output::open(string file, int s0_in, int ds_in)
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
  if (!MPISingle){
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  fid=H5Fcreate(file.c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
  H5Pclose(pid);

}


void Output::writeMeta()
{
  hid_t gid,gidsub;
  vector<double> tmp;
  tmp.resize(1);

  gid=H5Gcreate(fid,"Meta",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  gidsub=H5Gcreate(gid,"Version",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  tmp[0]=versionmajor;
  this->writeSingleNode(gidsub,"Major",&tmp);
  tmp[0]=versionminor;
  this->writeSingleNode(gidsub,"Minor",&tmp);
  tmp[0]=versionrevision;
  this->writeSingleNode(gidsub,"Revision",&tmp);
  tmp[0]=0;
  if (versionbeta) { tmp[0]=1;}
  this->writeSingleNode(gidsub,"Beta",&tmp);
  H5Gclose(gidsub);  
  
  time_t timer;
  time(&timer);
  string tim (ctime(&timer));
  this->writeSingleNodeString(gid,"TimeStamp", &tim);


  struct passwd *pws;
  pws = getpwuid(getuid());
  string user (pws->pw_name);
  this->writeSingleNodeString(gid,"User", &user);


  ifstream inFile (meta_inputfile->c_str());
  stringstream buffer;
  buffer << inFile.rdbuf();//read the file
  string str=buffer.str();
  this->writeSingleNodeString(gid,"InputFile", &str);
  inFile.close();

  ifstream inFile2 (meta_latfile->c_str());
  stringstream buffer2;
  buffer2 << inFile2.rdbuf();//read the file
  string str2=buffer2.str();
  this->writeSingleNodeString(gid,"LatticeFile", &str2);
  inFile2.close(); 



  H5Gclose(gid);


}

void Output::writeGlobal(double gamma, double lambda, double sample, double slen, bool one4one, bool time, bool scan)
{

  this->writeMeta();   
  vector<double> tmp;
  tmp.resize(1);
  hid_t gid;

  gid=H5Gcreate(fid,"Global",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  tmp[0]=gamma;
  this->writeSingleNode(gid,"gamma0",&tmp);
  tmp[0]=lambda;
  this->writeSingleNode(gid,"lambdaref",&tmp);
  tmp[0]=sample;
  this->writeSingleNode(gid,"sample",&tmp);
  tmp[0]=slen;
  this->writeSingleNode(gid,"slen",&tmp);
  tmp[0] = one4one ? 1. : 0 ;
  this->writeSingleNode(gid,"one4one",&tmp);
  tmp[0] = time ? 1. : 0 ;
  this->writeSingleNode(gid,"time",&tmp);
  tmp[0] = scan ? 1. : 0 ;
  this->writeSingleNode(gid,"scan",&tmp);

  H5Gclose(gid);

}



 
void Output::writeLattice(Beam * beam,Undulator *und)
{
  hid_t gid;
  gid=H5Gcreate(fid,"Lattice",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  this->writeSingleNode(gid,"zplot",&beam->zpos);
  this->writeSingleNode(gid,"z",&und->z);
  this->writeSingleNode(gid,"dz",&und->dz);
  this->writeSingleNode(gid,"aw",&und->aw);
  this->writeSingleNode(gid,"ax",&und->ax);
  this->writeSingleNode(gid,"ay",&und->ay);
  this->writeSingleNode(gid,"ku",&und->ku);
  this->writeSingleNode(gid,"kx",&und->kx);
  this->writeSingleNode(gid,"ky",&und->ky);
  this->writeSingleNode(gid,"qf",&und->qf);
  this->writeSingleNode(gid,"qx",&und->qx);
  this->writeSingleNode(gid,"qy",&und->qy);
  this->writeSingleNode(gid,"cx",&und->cx);
  this->writeSingleNode(gid,"cy",&und->cy);
  this->writeSingleNode(gid,"gradx",&und->gradx);
  this->writeSingleNode(gid,"grady",&und->grady);
  this->writeSingleNode(gid,"slippage",&und->slip);
  this->writeSingleNode(gid,"phaseshift",&und->phaseshift);
  this->writeSingleNode(gid,"chic_angle",&und->chic_angle);
  this->writeSingleNode(gid,"chic_lb",&und->chic_lb);
  this->writeSingleNode(gid,"chic_ld",&und->chic_ld);
  this->writeSingleNode(gid,"chic_lt",&und->chic_lt);

  H5Gclose(gid);


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

  int bh=beam->getBunchingHarmonics();
  char bgroup[20];
  for (int i=1; i<bh;i++){
    sprintf(bgroup,"bunching%d",(i+1));
    this->writeBuffer(gid, bgroup,  &beam->bh[i-1]);
    sprintf(bgroup,"bunchingphase%d",(i+1));
    this->writeBuffer(gid, bgroup,  &beam->ph[i-1]);
  }

  // step 3 - close group and done

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



  vector<double> tmp;
  tmp.resize(1);
  tmp[0]=field->dgrid;
  this->writeSingleNode(gid,"dgrid",&tmp);
  
  vector<int> itmp;
  itmp.resize(1);
  itmp[0]=field->ngrid;
  this->writeSingleNodeInt(gid,"ngrid",&itmp);



  // step 3 - close group and done






  H5Gclose(gid);
  return;

 

}


::::::::::::::
src/IO/readBeamHDF5.cpp
::::::::::::::
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
  

bool ReadBeamHDF5::readGlobal(int rank, int size,string file, Setup *setup, Time *time, double offset, bool dotime)
{


  isOpen=false;
  // read global data

  double reflen,slen;
  int count,nbins,one4one;


  fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);  
  readDataDouble(fid,(char *)"refposition",&s0,1);
  readDataDouble(fid,(char *)"slicelength",&reflen,1);
  readDataDouble(fid,(char *)"slicespacing",&slicelen,1);
  readDataInt(fid,(char *)"slicecount",&count,1);
  readDataInt(fid,(char *)"beamletsize",&nbins,1);
  readDataInt(fid,(char *)"one4one",&one4one,1);
  isOpen=true;


  double lambda=setup->getReferenceLength();                        // reference length for theta
  
  s0=s0+offset;  // add offset from input deck
  slen=slicelen*count;



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




  double ts0=time->getTimeWindowStart();
  double tslen=time->getTimeWindowLength(); 
  if ((ts0<s0) || (ts0+tslen)>(s0+slen)){
      if (rank==0){ cout << "*** Error: Defined time window is larger than given by input file" << endl;}
      return false;
  }

  // check if time window fits or whether it is defined by the external file


  return true; 


}




bool ReadBeamHDF5::readSlice(double s, vector<Particle> *slice, double *current){

  
  slice->resize(0);
  *current=0;

  if(!isOpen){ return false; } // skip if partfile option is not selected



  double rslice=(s-s0)/slicelen;

  if (fabs(rslice-round(rslice))>1e-3){ return false; }
  int islice=static_cast<int> (round(rslice))+1;


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
    slice->at(i).x=work[i];
  }

  sprintf(name,"slice%6.6d/y",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).y=work[i];
  }

  sprintf(name,"slice%6.6d/px",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).px=work[i];
  }

  sprintf(name,"slice%6.6d/py",islice);
  readDataDouble(fid,name,work,nsize);
  for (int i=0;i<nsize;i++){
    slice->at(i).py=work[i];
  }
  
  return true;
}


::::::::::::::
src/IO/readFieldHDF5.cpp
::::::::::::::
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

::::::::::::::
src/IO/SDDSBeam.cpp
::::::::::::::
#include "SDDSBeam.h"


SDDSBeam::SDDSBeam()
{
  file="";
  ds=0.01;
  xcen=0;
  ycen=0;
  pxcen=0;
  pycen=0;
  gamma=12000;
  betax=15;
  alphax=0;
  betay=15;
  alphay=0;
  charge=0;
  settime=false;
  center=false;
  match=false;
  matchs0=0;
  matchs1=1;
  align=0;
  aligns0=0;
  aligns1=0;
}

SDDSBeam::~SDDSBeam(){}

void SDDSBeam::usage(){

  cout << "List of keywords for importdistribution" << endl;
  cout << "&importdistribution" << endl;
  cout << " string file = <empty> " << endl;
  cout << " double charge   = 0 / <distribution file>" << endl;
  cout << " double slicewidth = 0.01" << endl;
  cout << " bool settimewindow = true" << endl;
  cout << " bool center = false " << endl;
  cout << " double gamma0 = gammaref " << endl;
  cout << " double x0 = 0 " << endl;
  cout << " double y0 = 0 " << endl;
  cout << " double px0 = 0 " << endl;
  cout << " double py0 = 0 " << endl;
  cout << " bool match = false " << endl;
  cout << " double betax  = 15 / matched" << endl;
  cout << " double alphax  = 0 / matched" << endl;
  cout << " double betay  = 15 / matched" << endl;
  cout << " double alphay  = 0 / matched" << endl;
  cout << " double match_start = 0 " << endl;
  cout << " double match_end = 1 " << endl;
  cout << " int align = 0 " << endl;
  cout << " double align_start = 0 " << endl;
  cout << " double align_end = 1 " << endl;
  cout << "&end" << endl << endl;
  return;
}

bool SDDSBeam::init(int inrank, int insize, map<string,string> *arg, Beam *beam, Setup *setup, Time *time, Lattice *lat)
{

  rank=inrank;
  size=insize;

  gamma=setup->getReferenceEnergy();           // get default energy from setup input deck
  lat->getMatchedOptics(&betax,&alphax,&betay,&alphay);  // use matched value if calculated

  double lambda=setup->getReferenceLength();   // reference length for theta
  double sample=static_cast<double>(time->getSampleRate());         // check slice length

  bool one4one=setup->getOne4One();
  bool shotnoise=setup->getShotNoise();
  int npart=setup->getNpart();
  int nbins=setup->getNbins();

  double theta0=4.*asin(1.);
  if (one4one) {
      nbins=1;
      theta0*=sample;
  }
  if ( (npart % nbins) != 0){
    if (rank==0) { cout << "*** Error: NPART is not a multiple of NBINS" << endl; } 
    return false;
  }
  
  theta0/=static_cast<double>(nbins);
 
  map<string,string>::iterator end=arg->end();

  if (arg->find("file")!=end)       {file   = arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("charge")!=end)     {charge = atof(arg->at("charge").c_str());     arg->erase(arg->find("charge"));}
  if (arg->find("slicewidth")!=end) {ds     = atof(arg->at("slicewidth").c_str()); arg->erase(arg->find("slicewidth"));}
  if (arg->find("match_start")!=end)    {matchs0    = atof(arg->at("match_start").c_str());    arg->erase(arg->find("match_start"));}
  if (arg->find("match_end")!=end)      {matchs1    = atof(arg->at("match_end").c_str());      arg->erase(arg->find("match_end"));}
  if (arg->find("align_start")!=end)    {aligns0    = atof(arg->at("align_start").c_str());    arg->erase(arg->find("align_start"));}
  if (arg->find("align_end")!=end)      {aligns1    = atof(arg->at("align_end").c_str());      arg->erase(arg->find("align_end"));}
  if (arg->find("betax")!=end)    {betax = atof(arg->at("betax").c_str()); arg->erase(arg->find("betax"));}
  if (arg->find("betay")!=end)    {betay = atof(arg->at("betay").c_str()); arg->erase(arg->find("betay"));}
  if (arg->find("alphax")!=end)   {alphax= atof(arg->at("alphax").c_str());arg->erase(arg->find("alphax"));}
  if (arg->find("alphay")!=end)   {alphay= atof(arg->at("alphay").c_str());arg->erase(arg->find("alphay"));}
  if (arg->find("x0")!=end)       {xcen  = atof(arg->at("x0").c_str());    arg->erase(arg->find("x0"));}
  if (arg->find("y0")!=end)       {ycen  = atof(arg->at("y0").c_str());    arg->erase(arg->find("y0"));}
  if (arg->find("px0")!=end)      {pxcen = atof(arg->at("px0").c_str());   arg->erase(arg->find("px0"));}
  if (arg->find("py0")!=end)      {pycen = atof(arg->at("py0").c_str());   arg->erase(arg->find("py0"));}
  if (arg->find("gamma0")!=end)   {gamma = atof(arg->at("gamma0").c_str());arg->erase(arg->find("gamma0"));}
  if (arg->find("align")!=end)    {align = atoi(arg->at("align").c_str()); arg->erase(arg->find("align"));}
  if (arg->find("match")!=end)    {match = atob(arg->at("match").c_str()); arg->erase(arg->find("match"));}
  if (arg->find("center")!=end)   {center= atob(arg->at("center").c_str());arg->erase(arg->find("center"));}
  if (arg->find("settimewindow")!=end)   {settime= atob(arg->at("settimewindow").c_str());arg->erase(arg->find("settimewindow"));}



  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &importdistribution" << endl; this->usage();}
    return false;
  }

  

  if (rank==0) { cout << "Importing distribution file... " << endl;}

 
  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  if (size>1){
     H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  hid_t fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,pid);
  H5Pclose(pid);

  // step 1 - check whether the field "charge" is present in the distribution

  int hascharge = H5Lexists(fid,"/charge",H5P_DEFAULT);
  if (hascharge>0){
    double extcharge;
    readDouble1D(fid,"charge" ,&extcharge,1,0);
    if (charge <=0) {
      charge=extcharge;
    }
  }
  if (rank==0) {cout << "Charge of external distribution: " << charge << endl; }


  // step 2 - get size of particle distribution

  string dset="t";
  int hasrecord = H5Lexists(fid,"/t",H5P_DEFAULT);
  if (hasrecord<=0){
    if (rank==0) { cout << "*** Error: Missing dataset in distribution file "<< file << endl;}
     H5Fclose(fid);
    return false;
  }

  int ntotal=getDatasetSize(fid, (char *)dset.c_str());
  double dQ=charge/static_cast<double> (ntotal);

  if (rank==0) {cout << "Particles in external distribution: " << ntotal << endl;}

 
  // step 3 - read datasets - this has to be checked that in parallel read the file can be not equivalent to the chunck size....

  int nchunk=ntotal/size;
  if ((ntotal % size) !=0) {nchunk++;}

  int nsize=nchunk;
  if ((rank*nchunk+nsize)>ntotal) { nsize=ntotal-rank*nchunk; }
  
  t.resize(nsize);
  g.resize(nsize);
  x.resize(nsize);
  y.resize(nsize);
  px.resize(nsize);
  py.resize(nsize);

  bool error=false;
  string dname="t";
  hasrecord=H5Lexists(fid,dname.c_str(),H5P_DEFAULT);
  if (hasrecord>0){    
     readDouble1D(fid,dname.c_str(),&t[0],nsize,rank*nchunk);
  } else {
    error=true;
  }
  dname="p";
  hasrecord=H5Lexists(fid,dname.c_str(),H5P_DEFAULT);
  if (hasrecord>0){    
     readDouble1D(fid,dname.c_str(),&g[0],nsize,rank*nchunk);
  } else {
    error=true;
  }
  dname="x";
  hasrecord=H5Lexists(fid,dname.c_str(),H5P_DEFAULT);
  if (hasrecord>0){    
     readDouble1D(fid,dname.c_str(),&x[0],nsize,rank*nchunk);
  } else {
    error=true;
  }
  dname="xp";
  hasrecord=H5Lexists(fid,dname.c_str(),H5P_DEFAULT);
  if (hasrecord>0){    
     readDouble1D(fid,dname.c_str(),&px[0],nsize,rank*nchunk);
  } else {
    error=true;
  }
  dname="y";
  hasrecord=H5Lexists(fid,dname.c_str(),H5P_DEFAULT);
  if (hasrecord>0){    
     readDouble1D(fid,dname.c_str(),&y[0],nsize,rank*nchunk);
  } else {
    error=true;
  }
  dname="yp";
  hasrecord=H5Lexists(fid,dname.c_str(),H5P_DEFAULT);
  if (hasrecord>0){    
     readDouble1D(fid,dname.c_str(),&py[0],nsize,rank*nchunk);
  } else {
    error=true;
  }
  H5Fclose(fid);

  if (error) {
    if (rank==0) { cout << "*** Error: Missing dataset in distribution file " << file << endl;} 
    return false;
  }


  // step 4 - analysing data file and setting up time window


  if (rank==0) { cout << "Analysing external distribution... " << endl;}


  for (int i=0; i<nsize; i++){
    t[i]*=-3e8;       // convert to positin in meters
    g[i]+=1.;         // convert from kinetic energy to total energy
  }


  double tmin,tmax;

  double tmp=*min_element(t.begin(),t.end());
  if (size==1){
    tmin=tmp;
  } else {
     MPI::COMM_WORLD.Allreduce(&tmp,&tmin,1,MPI::DOUBLE,MPI::MIN);
  }
  tmp=*max_element(t.begin(),t.end());
  if (size==1){
    tmax=tmp;
  } else {
     MPI::COMM_WORLD.Allreduce(&tmp,&tmax,1,MPI::DOUBLE,MPI::MAX);
  }

  double ttotal=tmax-tmin;

  for (int i=0; i<nsize; i++){
    t[i]-=tmin;
  }


  if (settime){
    time->setTimeWindowStart(0); 
    time->setTimeWindowLength(ttotal); 
    time->finishInit(setup);
  } 


  if (rank==0) {
    cout << "Analysis of the imported distribution" << endl;  
    cout << "   Total Bunch Length  (microns): " << ttotal*1e6 << endl;
  }

  this->analyse(ttotal,nsize);

 
  // step 5 - match/center distribution

  if (match) {
    if (rank==0){cout << "Matching external distribution..." << endl; }
    for (int i=0; i<nsize; i++){
      if (!center){ gamma=-gavg;}
      double ratio=sqrt(gavg/gamma);
      g[i]+=gamma-gavg;  // take out center so that the rematching is correct
      x[i]-=xavg;
      y[i]-=yavg;
      px[i]-=pxavg;
      py[i]-=pyavg;
      x[i]*=ratio;   // rescaling is needed to preserve emittance
      y[i]*=ratio;
      px[i]*=ratio;
      py[i]*=ratio;
      px[i]+=(ax/bx)*x[i];
      py[i]+=(ay/by)*y[i];
      x[i]*=sqrt(betax/bx);
      y[i]*=sqrt(betay/by);
      px[i]*=sqrt(bx/betax);
      py[i]*=sqrt(by/betay);
      px[i]-=(alphax/betax)*x[i];
      py[i]-=(alphay/betay)*y[i];
      g[i]-=gamma-gavg;  // apply initial offset again.
      x[i]+=xavg;
      y[i]+=yavg;
      px[i]+=pxavg;
      py[i]+=pyavg;
    }
  } 

  if (center) {
    if (rank==0){cout << "Centering external distribution..." << endl; }
    double ratio=sqrt(gavg/gamma);
    gavg=gamma-gavg;
    xavg=xcen-xavg;
    yavg=ycen-yavg;
    pxavg=pxcen-pxavg;
    pyavg=pycen-pyavg;
    for (int i=0; i<nsize; i++){
      g[i]+=gavg;
      x[i]+=xavg;
      y[i]+=yavg;
      px[i]+=pxavg;
      py[i]+=pyavg;
    }
  }



  if ((match)||(center)){
    if (rank==0) {cout << "Reanalysing matched and aligned distribution..." << endl; }
    this->analyse(ttotal,nsize);
  }


  // step 6 - sort distribution


  vector<double> s;
  int nslice=time->getPosition(&s);
  int node_off=time->getNodeOffset();
  int node_len=time->getNodeNSlice();
  beam->init(time->getNodeNSlice(),nbins,lambda,sample*lambda,s[0],one4one);  // init beam

  double smin=s[node_off];
  double smax=s[node_off+node_len-1];
  
  if (rank==0) {cout << "Sorting external distribution..." << endl; }

  double dslen=ds*ttotal;  // ds is the relative width to extract the samples (equivalent to 1/NDCUT)

  vector<vector<Particle> > dist;
  dist.resize(1);
  dist[0].clear();

  // copying all particles into the dist vector to enable sorting
  Particle part;
  for (int i=0; i<nsize; i++){
      part.theta=t[i];
      part.gamma=g[i];
      part.x=x[i];
      part.y=y[i];
      part.px=px[i]*g[i];
      part.py=py[i]*g[i];
      dist[0].push_back(part);
  }
  t.clear();
  g.clear();
  x.clear();
  y.clear();
  px.clear();
  py.clear();


  
  Sorting sort;
  sort.init(rank,size,false,true);  
  sort.configure(0,0,smin+0.5*dslen,smax-0.5*dslen,smin-0.5*dslen,smax+0.5*dslen,true); 
  sort.globalSort(&dist);


  // step 7 - populate internal distribution

  // now each node has all the particles, which is needed for the phase space reconstruction
  if (rank==0) {cout << "Generating internal particle distribution..." << endl; }


  this->initRandomSeq(setup->getSeed());
  ShotNoise sn;
  sn.init(setup->getSeed(),rank);

  int nwork=100;
  Particle *work;
  work=new Particle [nwork];
  for (int islice=0; islice<node_len;islice++){

    // step 1 - select all particles needed for the reconstruction of a given slice
    double sloc=s[islice+node_off];
    for (int i=0; i<dist[0].size();i++){
      if ((dist[0].at(i).theta>(sloc-0.5*dslen))&&(dist[0].at(i).theta<(sloc+0.5*dslen))){
	beam->beam.at(islice).push_back(dist[0].at(i));
      }
    }
   
    // step 2 - calculate the current and number of particles.
    int ncount = beam->beam.at(islice).size();
    int mpart;
    beam->current[islice]=static_cast<double>(ncount)*dQ*3e8/dslen;
    if (one4one){
      npart=static_cast<int>(round(beam->current[islice]*lambda*sample/ce));
      mpart=npart;
      nbins=1;
    } else {
      mpart=npart/nbins;
    }

    // step 3 - bring initial distribution to the right size

    if (beam->beam.at(islice).size() >= mpart){
      this->removeParticles(&beam->beam.at(islice),mpart);
    } else {
      this->addParticles(&beam->beam.at(islice),mpart);  // check for empty slices
    }

    // step 4 - refill particle phase completely new
    for (int i=0;i<beam->beam.at(islice).size();i++){
      beam->beam.at(islice).at(i).theta=theta0*ran->getElement();  // for one2one this should be the correct shot noise
    }

    if (!one4one){   // needs mirroring and shotnoise
      mpart=beam->beam.at(islice).size();
      beam->beam.at(islice).resize(mpart*nbins);
      for (int i=mpart; i>0; i--){
        int i1=i-1;
        int i2=nbins*i1;
        for (int j=0;j<nbins;j++){
          beam->beam.at(islice).at(i2+j).gamma=beam->beam.at(islice).at(i1).gamma;
  	  beam->beam.at(islice).at(i2+j).x    =beam->beam.at(islice).at(i1).x;
          beam->beam.at(islice).at(i2+j).y    =beam->beam.at(islice).at(i1).y;
          beam->beam.at(islice).at(i2+j).px   =beam->beam.at(islice).at(i1).px;
          beam->beam.at(islice).at(i2+j).py   =beam->beam.at(islice).at(i1).py;
          beam->beam.at(islice).at(i2+j).theta=beam->beam.at(islice).at(i1).theta+j*theta0;     
	}
      }
      double ne=round(beam->current[islice]*lambda*sample/ce);
      if (mpart*nbins>nwork){
      	nwork=mpart*nbins;
        delete[] work;
        work=new Particle [nwork];
      }
      for (int i=0;i<mpart*nbins;i++){
      	work[i].theta=beam->beam.at(islice).at(i).theta;  
      }
      sn.applyShotNoise(work,mpart*nbins,nbins,ne);
      for (int i=0;i<mpart*nbins;i++){
      	beam->beam.at(islice).at(i).theta=work[i].theta;  
      }

    }
   
  }


  dist[0].clear();
  delete ran;
  delete [] work;
  
  


  return true;


}


void SDDSBeam::addParticles(vector<Particle> *beam, int mpart){
   

  // check for error if there are only one or none particle to fill up distribution 
  int  ndist=beam->size();
  Particle par;

  if (ndist==0){
     par.x=0;
     par.y=0;
     par.px=0;
     par.py=0;
     par.theta=0;
     par.gamma=gamma;
     beam->push_back(par);
     ndist++;
  }
  if (ndist==1){   // mirror particle for algorithm to work 
    par.x    =beam->at(0).x;
    par.y    =beam->at(0).y;
    par.px   =beam->at(0).px;
    par.py   =beam->at(0).py;
    par.theta=beam->at(0).theta;
    par.gamma=beam->at(0).gamma;
    beam->push_back(par);
    ndist++;
  }

  // step 1 - calculate the center and the rms beam size
  double g1,x1,y1,px1,py1,g2,x2,y2,px2,py2;
  g1=0;
  x1=0;
  y1=0;
  px1=0;
  py1=0;
  g2=0;
  x2=0;
  y2=0;
  px2=0;
  py2=0;
  for (int i=0; i<ndist; i++){
    g1 +=beam->at(i).gamma;
    g2 +=beam->at(i).gamma * beam->at(i).gamma;
    x1 +=beam->at(i).x;
    x2 +=beam->at(i).x * beam->at(i).x;
    px1+=beam->at(i).px;
    px2+=beam->at(i).px * beam->at(i).px;
    y1 +=beam->at(i).y;
    y2 +=beam->at(i).y * beam->at(i).y;
    py1+=beam->at(i).py;
    py2+=beam->at(i).py * beam->at(i).py;
  }
  double scl=1/static_cast<double>(ndist);

  g1*=scl;
  g2=sqrt(fabs(g2*scl-g1*g1));
  x1*=scl;
  x2=sqrt(fabs(x2*scl-x1*x1));
  px1*=scl;
  px2=sqrt(fabs(px2*scl-px1*px1));
  y1*=scl;
  y2=sqrt(fabs(y2*scl-y1*y1));
  py1*=scl;
  py2=sqrt(fabs(py2*scl-py1*py1));

  // step 2 - invert the beam size for normalization and check for "cold" dimensions, e.g. zero energy spread
  if (g2==0) { g2=1; } else { g2=1/g2; }
  if (x2==0) { x2=1; } else { x2=1/x2; }
  if (y2==0) { y2=1; } else { y2=1/y2; }
  if (px2==0) { px2=1; } else { px2=1/px2; }
  if (py2==0) { py2=1; } else { py2=1/py2; }

  // step 3 - normalize distribution so that it is aligned to the origin and has an rms size of unity in all dimensions
  for (int i=0; i<ndist; i++){
    beam->at(i).gamma=(beam->at(i).gamma - g1)*g2;
    beam->at(i).x    =(beam->at(i).x     - x1)*x2;
    beam->at(i).y    =(beam->at(i).y     - y1)*y2;
    beam->at(i).px   =(beam->at(i).px    - px1)*px2;
    beam->at(i).py   =(beam->at(i).py    - py1)*py2;
  }
  
  // step 4 - add particles
  int ndist0=ndist;
  while (ndist<mpart){
    int n1=static_cast<int>(floor(static_cast<double>(ndist0)*ran->getElement()));
    double rmin=1e9;
    int n2=n1;
    for (int i=0; i<ndist0;i++){
       double r=this->distance(beam->at(n1),beam->at(i)); 
       if ((r<rmin) && ( i!=n1 )) {
           n2=i;
           rmin=r;
       }
    }
    par.gamma=0.5*(beam->at(n1).gamma+beam->at(n2).gamma)+(2*ran->getElement()-1)*(beam->at(n1).gamma-beam->at(n2).gamma);
    par.x =0.5*(beam->at(n1).x +beam->at(n2).x) +(2*ran->getElement()-1)*(beam->at(n1).x -beam->at(n2).x);
    par.px=0.5*(beam->at(n1).px+beam->at(n2).px)+(2*ran->getElement()-1)*(beam->at(n1).px-beam->at(n2).px);
    par.y =0.5*(beam->at(n1).y +beam->at(n2).y) +(2*ran->getElement()-1)*(beam->at(n1).y -beam->at(n2).y);
    par.py=0.5*(beam->at(n1).py+beam->at(n2).py)+(2*ran->getElement()-1)*(beam->at(n1).py-beam->at(n2).py);
    beam->push_back(par);    
    ndist++;
  }
 

  // step 5 - scale back

  for (int i=0; i<beam->size(); i++){
    beam->at(i).gamma=beam->at(i).gamma/g2 + g1;
    beam->at(i).x    =beam->at(i).x/x2 + x1;
    beam->at(i).y    =beam->at(i).y/y2 + y1;
    beam->at(i).px   =beam->at(i).px/px2 + px1; 
    beam->at(i).py   =beam->at(i).py/py2 + py1;
  }


  return;
}


double  SDDSBeam::distance(Particle p1, Particle p2){

  double tmp=p1.gamma-p2.gamma;  
  double r=tmp*tmp*ran->getElement();
  tmp=p1.x-p2.x;  
  r+=tmp*tmp*ran->getElement();
  tmp=p1.y-p2.y;  
  r+=tmp*tmp*ran->getElement();
  tmp=p1.px-p2.px;  
  r+=tmp*tmp*ran->getElement();
  tmp=p1.py-p2.py;  
  r+=tmp*tmp*ran->getElement();
  return r;
}


void SDDSBeam::removeParticles(vector<Particle> *beam,int mpart)
{
  int ndist=beam->size();
  while(ndist>mpart){
    int idx=static_cast<int>(floor(static_cast<double>(ndist)*ran->getElement()));
    beam->at(idx)=beam->at(ndist-1);
    ndist--;
  }
  beam->resize(mpart);
  return;
}

void SDDSBeam::analyse(double ttotal,int nsize)
{

  int ncount = 0;
  int nmean;
  double mt0=matchs0*ttotal;
  double mt1=matchs1*ttotal;

  double a2,a1,b2,b1,ab; 
  double c2,c1,d2,d1,cd; 
  double e1=0;
  a2=0;a1=0;b2=0;b1=0;ab=0;
  c2=0;c1=0;d2=0;d1=0;cd=0;
  
  //  cout << "Rank: " << rank << " Size: " << nsize << endl;

  for (int i=0; i <nsize;i++){
    if ((t[i]>mt0)&&(t[i]<mt1)){
      ncount++;
      a1+=x[i];
      a2+=x[i]*x[i];
      b1+=px[i];
      b2+=px[i]*px[i];
      ab+=x[i]*px[i];
      c1+=y[i];
      c2+=y[i]*y[i];
      d1+=py[i];
      d2+=py[i]*py[i];
      cd+=y[i]*py[i];
      e1+=g[i];
    }
  }

  if (size==1){
    nmean=ncount;
    gavg=e1;
    xavg=a1;
    pxavg=b1;
    yavg=c1;
    pyavg=e1;
    xvar=a2;
    pxvar=b2;
    xpx=ab;
    yvar=c2;
    pyvar=d2;
    ypy=cd;
  } else {
    MPI::COMM_WORLD.Allreduce(&ncount,&nmean,1,MPI::INT,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&e1,&gavg, 1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&a1,&xavg, 1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&b1,&pxavg,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&c1,&yavg, 1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&d1,&pyavg,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&a2,&xvar, 1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&b2,&pxvar,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&ab,&xpx,  1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&c2,&yvar, 1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&d2,&pyvar,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&cd,&ypy,  1,MPI::DOUBLE,MPI::SUM);
  }

  if (nmean>0){
    double tmp=1./static_cast<double>(nmean);
    gavg*=tmp;
    xavg*=tmp;
    yavg*=tmp;
    pxavg*=tmp;
    pyavg*=tmp;
    xvar*=tmp;
    yvar*=tmp;
    pxvar*=tmp;
    pyvar*=tmp;
    xpx*=tmp;
    ypy*=tmp;
  }

  ex=sqrt(fabs((xvar-xavg*xavg)*(pxvar-pxavg*pxavg)-(xpx-xavg*pxavg)*(xpx-xavg*pxavg)))*gavg;
  ey=sqrt(fabs((yvar-yavg*yavg)*(pyvar-pyavg*pyavg)-(ypy-yavg*pyavg)*(ypy-yavg*pyavg)))*gavg;
  bx=(xvar-xavg*xavg)/ex*gavg;
  by=(yvar-yavg*yavg)/ey*gavg;
  ax=-(xpx-xavg*pxavg)*gavg/ex;
  ay=-(ypy-yavg*pyavg)*gavg/ey;

  if (rank==0) {
       cout << "   Length for Matching (microns): " << (mt1-mt0)*1e6 << endl; 
       cout << "   Energy                  (MeV): " << gavg*eev*1e-6 << endl;
       cout << "   Norm. Emittance in x (micron): " << ex*1e6 << endl;
       cout << "   Norm. Emittance in y (micron): " << ey*1e6 << endl;
       cout << "   Beta Function in x        (m): " << bx << endl;
       cout << "   Beta Function in y        (m): " << by << endl;
       cout << "   Alpha Function in x          : " << ax << endl;
       cout << "   Alpha Function in y          : " << ay << endl;
       cout << "   Beam center in x     (micron): " << xavg*1e6 << endl;
       cout << "   Beam center in y     (micron): " << yavg*1e6 << endl;
       cout << "   Beam center in px            : " << pxavg << endl;
       cout << "   Beam center in py            : " << pyavg << endl;
  }
  return;


}

void SDDSBeam::initRandomSeq(int base)
{
     RandomU rseed(base);
     double val;
     for (int i=0; i<=rank+10000;i++){
        val=rseed.getElement();
     }
     val*=1e9;
     int locseed=static_cast<int> (round(val));
     ran  =  new RandomU (locseed);
     return;
}
::::::::::::::
src/IO/writeBeamHDF5.cpp
::::::::::::::

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





::::::::::::::
src/IO/writeFieldHDF5.cpp
::::::::::::::

#include "writeFieldHDF5.h"

extern bool MPISingle;

// constructor destructor
WriteFieldHDF5::WriteFieldHDF5()
{
}

WriteFieldHDF5::~WriteFieldHDF5()
{
}

void WriteFieldHDF5::write(string fileroot, vector<Field *> *field){

  string file;
  MPI::Status status;

  size=MPI::COMM_WORLD.Get_size(); // get size of cluster
  rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node
  if (MPISingle){
    size=1;
    rank=0;
  }

  for (int i=0; i<field->size();i++){
    int harm=field->at(i)->harm;
    char charm[10];
    sprintf(charm,".h%d",harm);
    if (harm==1){
      file=fileroot;
    } else {
      file=fileroot+string(charm);
    }
    this->writeMain(file,field->at(i));
  }

  return;
}





void WriteFieldHDF5::writeMain(string fileroot, Field *field){


 

  char filename[100];
  sprintf(filename,"%s.fld.h5",fileroot.c_str()); 
  if (rank == 0) { cout << "Writing field distribution to file: " <<filename << " ..." << endl;} 

  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  if (size>1){
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  fid=H5Fcreate(filename,H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
  H5Pclose(pid);

  s0=rank;
  int ntotal=size*field->field.size();

  // write global data
  this->writeGlobal(field->xlambda,field->slicelength,field->s0,field->dgrid,field->ngrid,ntotal);

  // loop through slices
  
  int smin=rank*field->field.size();
  int smax=smin+field->field.size();

  int ngrid=field->ngrid;
  vector<double> work;
  work.resize(ngrid*ngrid);

  double ks=4.*asin(1)/field->xlambda;
  double scl=field->dgrid*eev/ks/sqrt(vacimp);


  for (int i=0; i<ntotal;i++){
    s0=-1;
    char name[16];
    sprintf(name,"slice%6.6d",i+1);
    hid_t gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    if ((i>=smin) && (i<smax)){
      s0=0;    // select the slice which is writing
    }

    int islice= (i+field->first) % field->field.size() ;   // include the rotation due to slippage

    if (s0==0){
      for (int j=0; j<ngrid*ngrid;j++){ 
	work[j]=scl*field->field.at(islice).at(j).real();
      }  
    }
    this->writeSingleNode(gid,"field-real",&work);     

    if (s0==0){
      for (int j=0; j<ngrid*ngrid;j++){ 
	work[j]=scl*field->field.at(islice).at(j).imag();
      }  
    }
    this->writeSingleNode(gid,"field-imag",&work);     


        
    H5Gclose(gid);
  }


  H5Fclose(fid);

 
  return;
}

void WriteFieldHDF5::writeGlobal(double reflen, double slicelen, double s0, double dx, int nx, int count)
{

  
  vector<double> tmp;
  tmp.resize(1);

  tmp[0]=reflen;
  this->writeSingleNode(fid,"wavelength",&tmp);
  tmp[0]=slicelen;
  this->writeSingleNode(fid,"slicespacing",&tmp);
  tmp[0]=s0;
  this->writeSingleNode(fid,"refposition",&tmp);
  tmp[0]=dx;
  this->writeSingleNode(fid,"gridsize",&tmp);
  
  vector<int> itmp;
  itmp.resize(1);

  itmp[0]=nx;
  this->writeSingleNodeInt(fid,"gridpoints",&itmp);

  itmp[0]=count;
  this->writeSingleNodeInt(fid,"slicecount",&itmp);
  

  
  return;
}


::::::::::::::
src/Lattice/AlterLattice.cpp
::::::::::::::
#include "AlterLattice.h"

AlterLattice::AlterLattice()
{
  zmatch=0;
  err_aw=0;
  err_ax=0;
  err_ay=0;
  err_qx=0;
  err_qy=0;
  nlat=1;
}

AlterLattice::~AlterLattice(){}

void AlterLattice::usage(){

  cout << "List of keywords for Lattice" << endl;
  cout << "&lattice" << endl;
  cout << " double err_aw = 0" << endl;
  cout << " double err_ax = 0" << endl;
  cout << " double err_ay = 0" << endl;
  cout << " double err_qx = 0" << endl;
  cout << " double err_qy = 0" << endl;
  cout << " double zmatch = 0" << endl;
  cout << " int nlat = 1" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool AlterLattice::init(int inrank, int insize, map<string,string> *arg, Lattice *lat, Setup *setup)
{

  rank=inrank;
  size=insize;

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("zmatch")!=end)  {zmatch= atof(arg->at("zmatch").c_str());  arg->erase(arg->find("zmatch"));}
  if (arg->find("err_aw")!=end)  {err_aw= atof(arg->at("err_aw").c_str());  arg->erase(arg->find("err_aw"));}
  if (arg->find("err_ax")!=end)  {err_aw= atof(arg->at("err_ax").c_str());  arg->erase(arg->find("err_ax"));}
  if (arg->find("err_ay")!=end)  {err_aw= atof(arg->at("err_ay").c_str());  arg->erase(arg->find("err_ay"));}
  if (arg->find("err_qx")!=end)  {err_aw= atof(arg->at("err_qx").c_str());  arg->erase(arg->find("err_qx"));}
  if (arg->find("err_qy")!=end)  {err_aw= atof(arg->at("err_qy").c_str());  arg->erase(arg->find("err_qy"));}
  if (arg->find("nlat")!=end)    {nlat  = atof(arg->at("nlat").c_str());    arg->erase(arg->find("nlat"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &lattice" << endl; this->usage();}
    return false;
  }
  
  if (zmatch>0) {
    lat->match(rank, zmatch, setup->getReferenceEnergy());
  }

  return true;


}

::::::::::::::
src/Lattice/Lattice.cpp
::::::::::::::
#include "Lattice.h"

Lattice::Lattice()
{
  matched=false;
}


Lattice::~Lattice()
{
  for (int i=0;i<lat.size();i++){
    delete lat[i];
  }
  lat.clear();
}


bool Lattice::parse(string filename, string beamline, int rank)
{

  // release old lattice

  for (int i=0;i<lat.size();i++){
    delete lat[i];
  }
  lat.clear();

  LatticeParser parser;
  matched=false;

  if (rank == 0) { cout << "Parsing lattice file..." << endl; }
  bool err=parser.parse(filename,beamline,rank, lat);
  if (err==false) { 
    return err; 
  }
  
  layout.clear();
  
  for(int i=0; i<lat.size();i++){
    double z0=lat[i]->z;    
    if (lat[i]->type.compare("Drift")==0){
      continue;
    }
    layout[z0]=1;
    layout[z0+lat[i]->l]=1;
  }
  int last = lat.size()-1;
  if (lat[last]->type.compare("Drift")==0){
    double z0=lat[last]->z;
    layout[z0]=1;
    layout[z0+lat[last]->l]=1;
  }
  return true;
}


// generate from the input lattice the explicit representation, including some modification with 'alt'

bool Lattice::generateLattice(double delz, double lambda, double gamma, AlterLattice *alt,Undulator *und)
{

  this->unrollLattice(delz);  
  this->calcSlippage(lambda,gamma);

  und->setGammaRef(gamma);

  int ndata=lat_aw.size();

  und->aw.resize(ndata);
  und->ax.resize(ndata);
  und->ay.resize(ndata);
  und->ku.resize(ndata);
  und->kx.resize(ndata);
  und->ky.resize(ndata);
  und->gradx.resize(ndata);
  und->grady.resize(ndata);
  und->qf.resize(ndata);
  und->qx.resize(ndata);
  und->qy.resize(ndata);
  und->cx.resize(ndata);
  und->cy.resize(ndata);
  und->chic_angle.resize(ndata);
  und->chic_lb.resize(ndata);
  und->chic_ld.resize(ndata);
  und->chic_lt.resize(ndata);
  und->slip.resize(ndata);
  und->phaseshift.resize(ndata);
  und->z.resize(ndata);
  und->dz.resize(ndata);
  und->helical.resize(ndata);
  und->marker.resize(ndata+1);

  for (int i=0; i<ndata;i++){
      und->aw[i]=lat_aw[i];
      und->ax[i]=lat_ax[i];
      und->ay[i]=lat_ay[i];
      und->ku[i]=lat_ku[i];
      und->kx[i]=lat_kx[i];
      und->ky[i]=lat_ky[i];
      und->gradx[i]=lat_gradx[i];
      und->grady[i]=lat_grady[i];
      und->qf[i]=lat_qf[i];
      und->qx[i]=lat_qx[i];
      und->qy[i]=lat_qx[i];
      und->cx[i]=lat_cx[i];
      und->cy[i]=lat_cy[i];

      und->chic_angle[i]=lat_delay[i];  // here it is the delay but will converted to angle
      und->chic_lb[i]=lat_lb[i];
      und->chic_ld[i]=lat_ld[i];
      und->chic_lt[i]=lat_lt[i];       // here it is the total length but it will change to the time delay.
      
      if (und->chic_angle[i]!=0){
       double delay=fabs(und->chic_angle[i]);
       double tmin=0;
       double tmax=asin(1)-0.001;
       bool converged=false;
       double theta,d;
       while (!converged){
         theta=0.5*(tmax+tmin);
         d=4*und->chic_lb[i]*(theta/sin(theta)-1)+2*und->chic_ld[i]*(1/cos(theta)-1);
         if (d>delay) {
           tmax=theta;
         } else {
          tmin=theta;
	 }
         if (fabs(delay-d)<1e-15) { converged=true; }
       }
       und->chic_angle[i]=theta;
      }


      und->slip[i]=lat_slip[i];
      und->phaseshift[i]=lat_phase[i];
      und->z[i]=lat_z[i];
      und->dz[i]=lat_dz[i];
      und->helical[i]=lat_helical[i];
      und->marker[i]=lat_mk[i];

  }
  int lastmark=this->findMarker(und->z[ndata-1]+und->dz[ndata-1],"Marker");
  if (lastmark < 0){
    und->marker[ndata]=0;
  } else {
        Marker *mark=(Marker *)lat[lastmark];
        und->marker[ndata]=mark->action; 
  }


  return true;




}



void Lattice::calcSlippage(double lambda, double gamma)
{

  int nz=lat_aw.size();
  lat_slip.resize(nz);                                  // this needs to be improved for chicanes
  lat_phase.resize(nz);
  
  // calc the path length for a chicane

  double Lz=0;    // projected path
  double tmp;

  for (int i=0; i< nz;i++){
    if (lat_aw[i]>0){ // within undulator
      tmp=2*gamma*gamma*lambda/(1+lat_aw[i]*lat_aw[i]);
      lat_slip[i]=lat_dz[i]/tmp;
      lat_phase[i]=0;

      if ((Lz>0) && (i>0)){ // apply previous drift section
        tmp=Lz/2/gamma/gamma/lambda;  
        lat_slip[i-1]+=floor(tmp)+1;  //auto phasing would always add some slippage   
	//        tmp-=floor(tmp);
	//        lat_phase[i-1]+=tmp*4*asin(1);   // auto phasing
	Lz=0; // clear last intra beam section 
      } 
      // reset drift parameters
    } else {
      Lz   +=lat_dz[i];
      //      if (lat_delay[i]>0) { cout << "Delay: " << lat_delay[i] << " Lambda: " << lambda << endl; } 
      tmp=lat_delay[i]/lambda;  // affect of chicane is always autophasing!!!
      lat_slip[i]=floor(tmp);
      lat_phase[i]=lat_ps[i];     // phase shifter goes here
    }
  }
  // correct for the end of the lattice that autophasing is applied in case for second, suceeding run
  tmp=Lz/2/gamma/gamma/lambda;  
  lat_slip[nz-1]+=floor(tmp)+1;   
  //  tmp-=floor(tmp);
  // lat_phase[nz-1]+=tmp*4*asin(1);

}



int Lattice::findElement(double z0, double z1,string type)
{
  double zmid=0.5*(z0+z1);
  for(int i=0;i<lat.size();i++){
    double zz0 =lat[i]->z;
    double zz1 =zz0+lat[i]->l;
    if ((zmid>zz0)&&(zmid<zz1)&&(lat[i]->type.compare(type)==0)){
      return i;
    }
  }
  return -1; 
}

int Lattice::findMarker(double z0,string type)
{

  for(int i=0;i<lat.size();i++){
    double zz0 =lat[i]->z-z0;
    if ((zz0*zz0<1e-6) && (lat[i]->type.compare(type)==0)){  
      return i;
    }
  }
  return -1; 
}




void Lattice::getMatchedOptics(double *betax, double *alphax, double *betay, double *alphay)
{

  if (matched){
    *betax =mbetax;
    *alphax=malphax;
    *betay =mbetay;
    *alphay=malphay;
  }

  return;

}

void Lattice::match(int rank, double z0, double gammaref)
{

  double phix,phiy;

  this->unrollLattice(20); 
  Optics opt;
  opt.init();

  int i=0;

  while((lat_z[i]<z0)&&(i<lat_z.size())){
 
    double qf=lat_qf[i];
    double qx=lat_aw[i]*lat_aw[i]/gammaref/gammaref;
    double qy=qx*lat_ky[i];
    qx*=lat_kx[i];
    double dz=lat_dz[i];
    if ((lat_z[i]+lat_dz[i])>z0){ dz=z0-lat_z[i]; }
    opt.addElement(dz,qf,qx,qy);
    i++; 
  }

  if ((i==lat_z.size())&&((lat_z[i-1]+lat_dz[i-1])<z0)){
      if (rank==0){ cout << "*** Warning: Matching position beyond end of lattice. Ignoring matching command" << endl; }
      matched=false;
      return; 
  }


  if (matched=opt.match(&mbetax,&malphax,&mbetay,&malphay,&phix,&phiy)){
    if (rank==0){
        cout << "Matching for periodic solution between z = 0 and z = "<<z0 << " :" << endl;
        cout << "   betax (m) : " << mbetax  << endl;
        cout << "   alphax    : " << malphax << endl;
        cout << "   phix (deg): " << phix    << endl; 
        cout << "   betay (m) : " << mbetay  << endl;
        cout << "   alphay    : " << malphay << endl;
        cout << "   phiy (deg): " << phiy    << endl; 
    }
    matched=true;
  } else {
    if (rank==0){ cout << "*** Warning: Cannot find a matching solution" << endl; }
    matched=false;
  }
  return;
}
  

void Lattice::unrollLattice(double delz)
{

  lat_aw.clear();
  lat_ku.clear();
  lat_kx.clear();
  lat_ky.clear();
  lat_gradx.clear();
  lat_grady.clear();
  lat_ax.clear();
  lat_ay.clear();
  lat_helical.clear();
  lat_dz.clear();
  lat_z.clear();



  map<double,int>::iterator it=layout.begin();
  map<double,int>::iterator iend=layout.end();
  iend--;
  

  // first fill lattice with undulator

  while(it !=iend){
      //default

      double z0=it->first;
      it++;
      double z1=it->first;

      double dz=z1-z0;

      int iele=this->findElement(z0,z1,"Undulator");
      if (iele==-1){                                    // outside of an undulator
        lat_dz.push_back(z1-z0);
        lat_aw.push_back(0);
        lat_ku.push_back(0);
        lat_kx.push_back(0);
        lat_ky.push_back(0);
        lat_gradx.push_back(0);
        lat_grady.push_back(0);
        lat_ax.push_back(0);
        lat_ay.push_back(0);
        lat_helical.push_back(0);

      } else {
        ID *und=(ID *)lat[iele];
        int nz=round(dz/delz);
        if (nz==0) {nz=1;}
        dz=dz/static_cast<double>(nz);
        for (int iz=0;iz<nz;iz++){
	  lat_dz.push_back(dz);
          lat_aw.push_back(und->aw);
	  double ku=4.*asin(1)/und->lambdau;
          lat_ku.push_back(ku);
          lat_kx.push_back(ku*ku*und->kx);
          lat_ky.push_back(ku*ku*und->ky);
          lat_gradx.push_back(ku*und->gradx);
          lat_grady.push_back(ku*und->grady);
          lat_ax.push_back(und->ax);
          lat_ay.push_back(und->ay);
          lat_helical.push_back(static_cast<int>(und->helical));
        }
      }
  }

  int nz=lat_aw.size();
  lat_z.resize(nz);
  lat_z[0]=0; 
  for (int i=1; i<nz;i++){
    lat_z[i]=(lat_z[i-1]+lat_dz[i-1]);
  }   

  // resize vector and assign default zero values 
  lat_qf.resize(nz);
  lat_qx.resize(nz);
  lat_qy.resize(nz);
  lat_delay.resize(nz);
  lat_lb.resize(nz);
  lat_ld.resize(nz);
  lat_lt.resize(nz);
  lat_cx.resize(nz);
  lat_cy.resize(nz);
  lat_mk.resize(nz);
  lat_ps.resize(nz);
  for (int i=0; i<nz;i++){
    lat_qf[i]=0;
    lat_qx[i]=0;
    lat_qy[i]=0;
    lat_delay[i]=0;
    lat_lb[i]=0;
    lat_ld[i]=0;
    lat_lt[i]=0;
    lat_cx[i]=0;
    lat_cy[i]=0;
    lat_mk[i]=0;
    lat_ps[i]=0;
  }

  int lastChicane=-1;
  for (int i=0;i<nz;i++){
    double z0=lat_z[i];
    double z1=z0+lat_dz[i];

    bool inUnd=(lat_aw[i]>0);
    bool inQuad=0;

    int iele=this->findElement(z0,z1,"Quadrupole");

    if (iele!=-1){                                     // outside of an undulator
        Quadrupole *quad=(Quadrupole *)lat[iele];
        lat_qf[i]=quad->k1; 
        lat_qx[i]=quad->dx; 
        lat_qy[i]=quad->dy; 
	inQuad=true;
     }

    iele=this->findElement(z0,z1,"Chicane");
    if (iele!=-1){ 
      if (iele!=lastChicane){
	if ((!inUnd)&&(!inQuad)){
          Chicane *chicane=(Chicane *)lat[iele];
          lat_delay[i]=chicane->delay; 
          lat_lb[i]=chicane->lb; 
          lat_ld[i]=chicane->ld; 
          lat_lt[i]=chicane->l;    // save length of chicane in case that the chicane is split in parts
	  lastChicane=iele;
	} else {
	  cout << "*** Warning: Chicane inside Undulator or Quadrupole. Ignoring Chicane" << endl; 
	}
      }
    } 

    iele=this->findElement(z0,z1,"Corrector");
    if (iele!=-1){                                     // outside of a bend
        Corrector *cor=(Corrector *)lat[iele];
        lat_cx[i]=cor->cx; 
        lat_cy[i]=cor->cy; 
    } 

    iele=this->findElement(z0,z1,"Phaseshifter");
    if (iele!=-1){                                     // outside of a bend
        Phaseshifter *cor=(Phaseshifter *)lat[iele];
        lat_ps[i]=cor->phi; 
    } 

    // note that markers are only find at the beginning of the integration step. Thus a marker as the last element is ignored.
    iele=this->findMarker(z0,"Marker");
    if (iele!=-1){   
        Marker *mark=(Marker *)lat[iele];
        lat_mk[i]=mark->action; 
    } 
  }

  return;  

}



























::::::::::::::
src/Lattice/LatticeElements.cpp
::::::::::::::

::::::::::::::
src/Lattice/LatticeParser.cpp
::::::::::::::
#include "LatticeParser.h"

LatticeParser::LatticeParser()
{
}


LatticeParser::~LatticeParser()
{
}


bool LatticeParser::parse(string file, string line, int rank, vector<Element *> &lat)
{

  ifstream fin(file.c_str(),ios_base::in);

  if (!fin){
    if (rank==0) {cout << "*** Error: Cannot open magnetic lattice file: " << file << endl;}
    return false;
  }

  //------------------------------------------------------
  // step one - coarse parsing of the input deck

  string instring;
  string comstring="";
  vector<string> content;

  while(getline(fin,instring,'\n')){    // read line
    this->trim(instring);
    if ((!instring.compare(0,1,"#") || instring.length() < 1)){ continue; } // skip comment and empty rows
    comstring.append(" ");
    comstring.append(instring);  // add all content into one string
  }
  fin.close();

  for (int i=0; i<comstring.size();i++){ // convert to lower case
    comstring[i]=tolower(comstring[i]);
  }

  size_t pos;
  while ((pos=comstring.find_first_of(";")) !=string::npos){  // split into individual lines
     instring=comstring.substr(0,pos);
     this->trim(instring);
     content.push_back(instring);
     comstring.erase(0,pos+1);
  }
  
  this->trim(comstring);
  if ((comstring.length()>0) && (rank==0)){
    cout << "*** Warning: Ignoring incomplete last line in magnetic lattice" << endl;
  } 

  // -----------------------------------------------------------------------
  // step two - parse each individual line of the lattice according to the format
  // label: type =(content);
  // e.g.   "QF1: Quadrupole = {L=0.2, k1=0.8}; 

  label.clear();
  type.clear();
  argument.clear();

  string inlabel, intype, inargument;
  bool error=false;
 
  for(int i=0;i<content.size();i++){
    
    if ((pos=content[i].find_first_of(':'))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }
    inlabel=content[i].substr(0,pos);
    this->trim(inlabel);
    content[i].erase(0,pos+1);
    if ((pos=content[i].find_first_of('='))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }
    intype=content[i].substr(0,pos);
    this->trim(intype);
    content[i].erase(0,pos+1);
    
    if ((pos=content[i].find_first_of('{'))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }
    content[i].erase(0,pos+1);
    if ((pos=content[i].find_first_of('}'))==string::npos){
      if (rank==0){ cout<< "*** Error: Invalid Format in lattice file: " << content[i] << endl;}
      error=true;
      continue;
    }

    content[i].erase(pos,content[i].size());
    inargument=content[i];
    this->trim(inargument);
    label.push_back(inlabel);
    type.push_back(intype.substr(0,4));
    argument.push_back(inargument);
    
  }

  if (error){ return false; }

  // -------------------------------------------------------------
  // step 3 - resolving all references

  
  int recursion = 20;
  for (int i=0;i<label.size();i++){
    if (type[i].compare("line")!=0){
      error=this->resolve(i,recursion-1,rank);
      if (error==false){ return false; }
    }   
  } 


  // --------------------------------------------------------------
  // step 4 - unrolling the line
  
  for (int i=0;i<line.size();i++){
    line[i]=tolower(line[i]);      
  }

  int idx=this->findIndex(&label,line);

  sequence.clear();
  zref.clear();
  refele=-1;
  
  if ((idx>-1) && (type[idx].compare("line")==0)) {
    error=this->unroll(idx, recursion-1, rank);
  } else {
      if (rank==0) {cout << "*** Error: Lattice file does not contain the beamline: " << line << endl;}
      return false;  
  }


  if (error==false){ return false; }

  //---------------------------------------------------------------------------------
  // step 5 - initiate

  double z;
  
  for (int i=0;i<sequence.size();i++){
    if (zref[i]<0){
      z=0;
    } else {
      z=lat[zref[i]]->z+lat[zref[i]]->l;
    } 
    idx=this->findIndex(&label,sequence[i]);
    if (type[idx].compare("quad")==0){ error=false; lat.push_back(this->parseQuad(idx,rank,z));}
    if (type[idx].compare("undu")==0){ error=false; lat.push_back(this->parseID(idx,rank,z)); }
    if (type[idx].compare("drif")==0){ error=false; lat.push_back(this->parseDrift(idx,rank,z));}
    if (type[idx].compare("corr")==0){ error=false; lat.push_back(this->parseCorrector(idx,rank,z));}
    if (type[idx].compare("chic")==0){ error=false; lat.push_back(this->parseChicane(idx,rank,z)); }
    if (type[idx].compare("mark")==0){ error=false; lat.push_back(this->parseMarker(idx,rank,z)); }
    if (type[idx].compare("phas")==0){ error=false; lat.push_back(this->parsePhaseshifter(idx,rank,z)); }
    if (error) { 
      if (rank==0) {cout << "*** Error: Unknown element type " << type[idx] << " in lattice" << endl;}
      return false;
    }
  }
  return true;
}


ID *LatticeParser::parseID(int idx,int rank, double zin)
{
  ID *ele=new ID;
  
  ele->type="Undulator";
  ele->z=zin;
  ele->aw=0;
  ele->lambdau=0;
  ele->nwig=0;
  ele->kx=0;
  ele->ky=1;
  ele->ax=0;
  ele->ay=0;
  ele->gradx=0;
  ele->grady=0;
  ele->helical=false;

  ele->paw=0;
  ele->pkx=0;
  ele->pky=1;
  ele->pdadx=0;
  ele->pdady=0;
  ele->phase=0;


  ele->helical=false;

  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0)    { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("lambdau")==0){ ele->lambdau=atof(val.c_str()); found=true; };
    if (fld.compare("aw")==0)   { ele->aw=atof(val.c_str()); found=true; };
    if (fld.compare("aw_perp")==0)  { ele->paw=atof(val.c_str()); found=true; };
    if (fld.compare("nwig")==0) { ele->nwig=atof(val.c_str()); found=true; };
    if (fld.compare("kx")==0)   { ele->kx=atof(val.c_str()); found=true; };
    if (fld.compare("ky")==0)   { ele->ky=atof(val.c_str()); found=true; };
    if (fld.compare("kx_perp")==0)  { ele->pkx=atof(val.c_str()); found=true; };
    if (fld.compare("ky_perp")==0)  { ele->pky=atof(val.c_str()); found=true; };
    if (fld.compare("ax")==0)   { ele->ax=atof(val.c_str()); found=true; };
    if (fld.compare("ay")==0)   { ele->ay=atof(val.c_str()); found=true; };
    if (fld.compare("gradx")==0) { ele->gradx=atof(val.c_str()); found=true; };
    if (fld.compare("grady")==0) { ele->grady=atof(val.c_str()); found=true; };
    if (fld.compare("gradx_perp")==0){ ele->pdadx=atof(val.c_str()); found=true; };
    if (fld.compare("grady_perp")==0){ ele->pdady=atof(val.c_str()); found=true; };
    if (fld.compare("phase_perp")==0){ ele->phase=atof(val.c_str()); found=true; };
    if (fld.compare("helical")==0)   { ele->helical=atob(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknow parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  ele->l=ele->lambdau*ele->nwig;
  return ele;
}


Corrector *LatticeParser::parseCorrector(int idx,int rank, double zin)
{
  Corrector *ele=new Corrector;

  ele->type="Corrector";
  ele->z=zin;
  ele->l=0;
  ele->cx=0;
  ele->cy=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0)    { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("cx")==0)   { ele->cx=atof(val.c_str()); found=true; };
    if (fld.compare("cy")==0)   { ele->cy=atof(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknow parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}


Chicane *LatticeParser::parseChicane(int idx,int rank, double zin)
{
  Chicane *ele=new Chicane;

  ele->type="Chicane";
  ele->z=zin;
  ele->l=0;
  ele->delay=0;
  ele->lb=0;
  ele->ld=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0)    { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("delay")==0){ ele->delay=atof(val.c_str()); found=true; };
    if (fld.compare("lb")==0)   { ele->lb=atof(val.c_str()); found=true; };
    if (fld.compare("ld")==0)   { ele->ld=atof(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknow parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}

Marker *LatticeParser::parseMarker(int idx,int rank, double zin)
{
  Marker *ele=new Marker;

  ele->type="Marker";
  ele->z=zin;
  ele->l=0;
  ele->action=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    int tag=atoi(val.c_str());
    if (tag !=0) { tag=1;}
    if (fld.compare("dumpfield")==0){ele->action|=1*tag;found=true;}
    if (fld.compare("dumpbeam")==0){ele->action|=2*tag;found=true;}
    if (fld.compare("sort")==0){ele->action|=4*tag;found=true;}
    if (fld.compare("stop")==0){ele->action|=8*tag;found=true;}
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknow parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}



Drift *LatticeParser::parseDrift(int idx,int rank, double zin)
{
  Drift *ele=new Drift;

  ele->type="Drift";
  ele->z=zin;
  ele->l=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0) { ele->l=atof(val.c_str());  found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknow parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}


Quadrupole *LatticeParser::parseQuad(int idx,int rank, double zin)
{
  Quadrupole *ele=new Quadrupole;

  ele->type="Quadrupole";
  ele->z=zin;
  ele->l=0;
  ele->k1=0;
  ele->dx=0;
  ele->dy=0;
  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0) { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("k1")==0){ ele->k1=atof(val.c_str()); found=true; };
    if (fld.compare("dx")==0){ ele->dx=atof(val.c_str()); found=true; };
    if (fld.compare("dy")==0){ ele->dy=atof(val.c_str()); found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknow parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}

Phaseshifter *LatticeParser::parsePhaseshifter(int idx,int rank, double zin)
{
  Phaseshifter *ele=new Phaseshifter;

  ele->type="Phaseshifter";
  ele->z=zin;
  ele->l=0;
  ele->phi=0;

  vector <string> par;
  string fld,val;

  this->chop(argument[idx],&par);
  for (int i=0;i<par.size();i++){
    size_t pos=par[i].find_first_of("=");
    if (pos==string::npos){
      if (rank==0){cout << "*** Warning: Ignoring invalid format: " << par[i] << " for element " << label[idx]<< endl;}
      continue;  
    }
    fld=par[i].substr(0,pos);
    val=par[i].erase(0,pos+1);
    this->trim(fld);
    bool found=false;
    if (fld.compare("l")==0) { ele->l=atof(val.c_str());  found=true; };
    if (fld.compare("phi")==0){ele->phi=atof(val.c_str());found=true; };
    if (found==false){
      if (rank==0){cout << "*** Warning: Ignoring unknow parameter: " << fld << " for element " << label[idx]<< endl;}
    }
  }
  return ele;
}

bool LatticeParser::resolve(int idx, int recursion, int rank)
{
  if (recursion < 0){
    if (rank==0){  cout << "*** Error: Too many nested references" << endl;}
   return false;
  }
  
  size_t pos=argument[idx].find("ref");
  if (pos==string::npos){ return true;}
  vector<string> arg1;
  this->chop(argument[idx], &arg1);

  string ref;
  for (int i=0; i< arg1.size(); i++){
    pos=arg1[i].find("ref");
    if (pos !=string::npos){
      pos=arg1[i].find_first_of("=");
      ref=arg1[1].erase(0,pos+1);
      arg1.erase(arg1.begin()+i);
      break;
    }
  } 
  this->trim(ref); 
  int iref=findIndex(&label,ref);
  if (iref<0){
      if (rank==0) {cout << "*** Error: Unresolved reference: " << ref <<" for element: " << label[idx] << endl;}
      return false;
  }  

  if (this->resolve(iref,recursion-1,rank)==false) {return false;}

  argument[idx]=argument[iref];  
  for (int i=0; i< arg1.size(); i++){
    argument[idx].append(",");
    argument[idx].append(arg1[i]);
  }
 
  return true;



}

bool LatticeParser::unroll(int idx, int recursion,int rank){

 
  if (recursion < 0){
    if (rank==0){  cout << "*** Error: Too many nested elements in selected beamline" << endl;}
   return false;
  }
  vector<string> line;
  this->chop(argument[idx],&line);

  int refelesave=refele;

  for (int i=0; i<line.size();i++){

    int count = checkMultiplier(&line[i]);
    bool resetpos = checkResetPosition(&line[i]);
    int ix=this->findIndex(&label,line[i]);
    if (ix <0){
      if (rank==0) {cout << "*** Error: Undefined element in beamline: " << line[i] << endl;}
      return false;
    }
    for (int j=0; j < count ; j++){
      if (type[ix].compare("line")==0){
        bool error=this->unroll(ix,recursion-1,rank);
        if (error==false) { return error; }
        if (resetpos) { refele=refelesave; }
      } else {        
        sequence.push_back(line[i]); 
        zref.push_back(refele);
        refele++;
        if (resetpos) { refele--; }

      }
    }
  }

  line.clear();
  return true;

}


int LatticeParser::checkMultiplier(string *element){

  size_t pos=element->find_first_of("*");
  if (pos==string::npos){
    return 1;
  } else {
    string num=element->substr(0,pos);
    element->erase(0,pos+1);
    this->trim(num);
    this->trim(*element);
    return atoi(num.c_str());
  }
}

bool LatticeParser::checkResetPosition(string *element){
  size_t pos=element->find_first_of("@");
  if (pos==string::npos){
    return false;
  }
  element->erase(pos);
  this->trim(*element);
  return true; 
} 

int LatticeParser::findIndex(vector<string> *list, string element){

  for (int i=list->size()-1;i>-1;i--){
    if (list->at(i).compare(element)==0) { return i; }
  }
    return -1;
}
::::::::::::::
src/Lattice/Optics.cpp
::::::::::::::
#include "Optics.h"

Optics::Optics(){}
Optics::~Optics(){}

//----------- migrate into optics class

void Optics::init(){
  Dx11=1;
  Dx12=0;
  Dx21=0;
  Dx22=1;
  Dy11=1;
  Dy12=0;
  Dy21=0;
  Dy22=1;
  return;
}

bool Optics::match(double *bx, double *ax, double *by, double *ay,double *phix, double *phiy)
{
   
  double arg=2-Dx11*Dx11-2*Dx12*Dx21-Dx22*Dx22;
  if (arg<=0){ return false;}
  *bx=2*Dx12/sqrt(arg);
  *ax=(Dx11-Dx22)/sqrt(arg);
  *phix=acos(0.5*(Dx11+Dx22))*90/asin(1);
  arg=2-Dy11*Dy11-2*Dy12*Dy21-Dy22*Dy22;
  if (arg<=0){ return false;}
  *by=2*Dy12/sqrt(arg);
  *ay=(Dy11-Dy22)/sqrt(arg);
  *phiy=acos(0.5*(Dy11+Dy22))*90./asin(1.);

  return true;


}

void Optics::addElement(double dz, double qf, double qx, double qy)
{
  
  if ((qf+qx)==0){
    getDrift(dz);
  } else {
    getQuad(qf+qx,dz);
  }  
  this->MatMult(true);

  if ((-qf+qy)==0){
    getDrift(dz);
  } else {
    getQuad(-qf+qy,dz);
  }  
  this->MatMult(false);

  return;
}


void Optics::getDrift(double L)
{
  
  M11=1;
  M12=L;
  M21=0;
  M22=1;
  return;
}


void Optics::getQuad(double k1,double L)
{
  if (k1==0){
    this->getDrift(L);
    return;
  }
  double omg=sqrt(fabs(k1))*L;
  if (k1>0){
    M11=cos(omg);
    M12=sin(omg)/sqrt(fabs(k1));
    M21=-sin(omg)*sqrt(fabs(k1));
    M22=cos(omg);
  } else {
    M11=cosh(omg);
    M12=sinh(omg)/sqrt(fabs(k1));
    M21=sinh(omg)*sqrt(fabs(k1));
    M22=cosh(omg);
  }
  return;
}


void Optics::MatMult(bool isX)
{
  
  double D11,D12,D21,D22;
  if (isX){
   D11=Dx11;
   D12=Dx12;
   D21=Dx21;
   D22=Dx22;
   Dx11=M11*D11+M12*D21;
   Dx12=M11*D12+M12*D22;
   Dx21=M21*D11+M22*D21;
   Dx22=M21*D12+M22*D22;
  } else {
   D11=Dy11;
   D12=Dy12;
   D21=Dy21;
   D22=Dy22;
   Dy11=M11*D11+M12*D21;
   Dy12=M11*D12+M12*D22;
   Dy21=M21*D11+M22*D21;
   Dy22=M21*D12+M22*D22;
  }
 
  return;
}
::::::::::::::
src/Loading/ImportBeam.cpp
::::::::::::::
#include "ImportBeam.h"
#include "readBeamHDF5.h"


ImportBeam::ImportBeam()
{
  offset=0;
  dotime=true;
}

ImportBeam::~ImportBeam(){}

void ImportBeam::usage(){

  cout << "List of keywords for IMPORTBEAM" << endl;
  cout << "&importbeam" << endl;
  cout << " string file = <empty>" << endl;
  cout << " double offset = 0" << endl;
  cout << " bool time = true" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool ImportBeam::init(int rank, int size, map<string,string> *arg, Beam *beam, Setup *setup, Time *time)
{


  if (beam->beam.size()>0){
    if (rank==0) {cout << "*** Error: Cannot import beam, because beam is already defined" << endl; }
    return false;
  }
    
  double lambda=setup->getReferenceLength();   // reference length for theta
  bool one4one=setup->getOne4One();            // check for one4one simulations
  double gamma=setup->getReferenceEnergy();           // get default energy from setup input deck
  int nbins=setup->getNbins();


  map<string,string>::iterator end=arg->end();

  if (arg->find("file")!=end    ){file=arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("offset")!=end  ){offset=atof(arg->at("offset").c_str());  arg->erase(arg->find("offset"));}
  if (arg->find("time")!=end)    {dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &importbeam" << endl; this->usage();}
    return false;
  }

  

  if (rank==0){
    cout << "Importing particle distribution from file: " << file << " ..." << endl; 
  }

  ReadBeamHDF5 import;

  bool check=import.readGlobal(rank, size, file, setup, time, offset,dotime);
  if (!check) { 
    import.close();
    return check; 
  }



  // sample rate and time dependent run could have changed when taken by externaldistribution in readGlobal
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  dotime=time->isTime();                                            // check for time simulation

  vector<double> s;
  int nslice=time->getPosition(&s);
  beam->init(time->getNodeNSlice(),nbins,lambda,sample*lambda,s[0],one4one);



  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    double sloc=s[i];
    import.readSlice(s[i],&beam->beam[j],&beam->current[j]);
  }
  import.close();

  return true;

}


::::::::::::::
src/Loading/ImportField.cpp
::::::::::::::
#include "ImportField.h"
#include "readFieldHDF5.h"


ImportField::ImportField()
{
  offset=0;
  dotime=true;
}

ImportField::~ImportField(){}

void ImportField::usage(){

  cout << "List of keywords for IMPORTFIELD" << endl;
  cout << "&importfield" << endl;
  cout << " string file = <empty>" << endl;
  cout << " double offset = 0" << endl;
  cout << " bool time = true" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool ImportField::init(int rank, int size, map<string,string> *arg, vector<Field *> *beam, Setup *setup, Time *time)
{

  
  if (beam->size()>0){
    if (rank==0) {cout << "*** Error: Cannot import field, because field is already defined" << endl; }
    return false;
  }
    
  double lambda=setup->getReferenceLength();   // reference length for theta
  double gamma=setup->getReferenceEnergy();           // get default energy from setup input deck


  map<string,string>::iterator end=arg->end();

  if (arg->find("file")!=end    ){file=arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("offset")!=end  ){offset=atof(arg->at("offset").c_str());  arg->erase(arg->find("offset"));}
  if (arg->find("time")!=end)    {dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &importfield" << endl; this->usage();}
    return false;
  }

  if (rank==0){
    cout << "Importing field distribution from file: " << file << " ..." << endl; 
  }

  ReadFieldHDF5 import;

  bool check=import.readGlobal(rank, size, file, setup, time, offset,dotime);
  if (!check) { 
    import.close();
    return check; 
  }

  /*

  // sample rate and time dependent run could have changed when taken by externaldistribution in readGlobal
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  dotime=time->isTime();                                            // check for time simulation

  vector<double> s;
  int nslice=time->getPosition(&s);
  beam->init(time->getNodeNSlice(),nbins,lambda,sample*lambda,s[0],one4one);



  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    double sloc=s[i];
    import.readSlice(s[i],&beam->beam[j],&beam->current[j]);
  }
  import.close();
  
  */ 

  return true;

}



::::::::::::::
src/Loading/LoadBeam.cpp
::::::::::::::
#include "LoadBeam.h"

LoadBeam::LoadBeam()
{
  gamma=5800/0.511; gammaref="";
  delgam=0; delgamref="";
  current=0; currentref="";
  ex=0.3e-6; exref="";
  ey=0.3e-6; eyref="";
  betax=15; betaxref="";
  betay=15; betayref="";
  alphax=0; alphaxref="";
  alphay=0; alphayref="";
  xcen=0; xcenref="";
  ycen=0; ycenref="";
  pxcen=0; pxcenref="";
  pycen=0; pycenref="";
  bunch=0; bunchref="";
  bunchphase=0; bunchphaseref="";
  emod=0; emodref="";
  emodphase=0; emodphaseref="";
}

LoadBeam::~LoadBeam(){}

void LoadBeam::usage(){

  cout << "List of keywords for BEAM" << endl;
  cout << "&beam" << endl;
  cout << " double gamma = gammaref / reference" << endl;
  cout << " double delgam = 0 / reference" << endl;
  cout << " double current = 1000 / reference" << endl;
  cout << " double ex = 0.3e-6 / reference" << endl;
  cout << " double ey = 0.3e-6 / reference" << endl;
  cout << " double betax = 15 / reference / matched" << endl;
  cout << " double betay = 15 / reference / matched" << endl;
  cout << " double alphax = 0 / reference / matched" << endl;
  cout << " double alphay = 0 / reference / matched" << endl;
  cout << " double xcenter = 0 / reference" << endl;
  cout << " double ycenter = 0 / reference" << endl;
  cout << " double pxcenter = 0 / reference" << endl;
  cout << " double pycenter = 0 / reference" << endl;
  cout << " double bunch = 0 / reference" << endl;
  cout << " double bunchphase = 0 / reference" << endl;
  cout << " double emod = 0 / reference" << endl;
  cout << " double emodphase = 0 / reference" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool LoadBeam::init(int rank, int size, map<string,string> *arg, Beam *beam, Setup *setup, Time *time, Profile *prof, Lattice *lat)
{


  if (beam->beam.size()>0){
    if (rank==0) {cout << "*** Error: Cannot generat beam, because beam is already defined" << endl; }
    return false;
  }


  double lambda=setup->getReferenceLength();   // reference length for theta
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  bool one4one=setup->getOne4One();            // check for one4one simulations
  bool shotnoise=setup->getShotNoise();
  int npart=setup->getNpart();
  int nbins=setup->getNbins();
  bool dotime=time->isTime();                  // check for time simulation

  gamma=setup->getReferenceEnergy();           // get default energy from setup input deck
  lat->getMatchedOptics(&betax,&alphax,&betay,&alphay);  // use matched value if calculated

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("gamma")!=end    ){this->reference(arg->at("gamma"),&gamma,&gammaref); arg->erase(arg->find("gamma"));}
  if (arg->find("delgam")!=end   ){this->reference(arg->at("delgam"),&delgam,&delgamref); arg->erase(arg->find("delgam"));}
  if (arg->find("current")!=end  ){this->reference(arg->at("current"),&current,&currentref); arg->erase(arg->find("current"));}
  if (arg->find("ex")!=end       ){this->reference(arg->at("ex"),&ex,&exref); arg->erase(arg->find("ex"));}
  if (arg->find("ey")!=end       ){this->reference(arg->at("ey"),&ey,&eyref); arg->erase(arg->find("ey"));}
  if (arg->find("betax")!=end    ){this->reference(arg->at("betax"),&betax,&betaxref); arg->erase(arg->find("betax"));}
  if (arg->find("betay")!=end    ){this->reference(arg->at("betay"),&betay,&betayref); arg->erase(arg->find("betay"));}
  if (arg->find("alphax")!=end   ){this->reference(arg->at("alphax"),&alphax,&alphaxref); arg->erase(arg->find("alphax"));}
  if (arg->find("alphay")!=end   ){this->reference(arg->at("alphay"),&alphay,&alphayref); arg->erase(arg->find("alphay"));}
  if (arg->find("xcenter")!=end  ){this->reference(arg->at("xcenter"),&xcen,&xcenref); arg->erase(arg->find("xcenter"));}
  if (arg->find("ycenter")!=end  ){this->reference(arg->at("ycenter"),&ycen,&ycenref); arg->erase(arg->find("ycenter"));}
  if (arg->find("pxcenter")!=end ){this->reference(arg->at("pxcenter"),&pxcen,&pxcenref); arg->erase(arg->find("pxcenter"));}
  if (arg->find("pycenter")!=end ){this->reference(arg->at("pycenter"),&pycen,&pycenref); arg->erase(arg->find("pycenter"));}
  if (arg->find("bunch")!=end    ){this->reference(arg->at("bunch"),&bunch,&bunchref); arg->erase(arg->find("bunch"));}
  if (arg->find("bunchphase")!=end ){this->reference(arg->at("bunchphase"),&bunchphase,&bunchphaseref); arg->erase(arg->find("bunchphase"));}
  if (arg->find("emod")!=end     ){this->reference(arg->at("emod"),&emod,&emodref); arg->erase(arg->find("emod"));}
  if (arg->find("emodphase")!=end ){this->reference(arg->at("emodphase"),&emodphase,&emodphaseref); arg->erase(arg->find("emodphase"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &beam" << endl; this->usage();}
    return false;
  }



  if (rank==0){cout << "Generating input particle distribution..." << endl; }


  double theta0=4.*asin(1.);

  if (one4one){
    nbins=1;
    theta0*=sample;
  }  
  
  if ( (npart % nbins) != 0){
    if (rank==0) { cout << "*** Error: NPART is not a multiple of NBINS" << endl; } 
    return false;
  }
 

  vector<double> s;
  int nslice=time->getPosition(&s);

  beam->init(time->getNodeNSlice(),nbins,lambda,sample*lambda,s[0],one4one);
  beam->initSorting(rank,size,false,one4one);  // sorting routine is initialized, with default values to suppress field slippage but do sorting if one4one is enabled

  int nbeam=1024;
  Particle *beamslice = new Particle [nbeam]; 
  BeamSlice slice;


  QuietLoading ql;
  ShotNoise sn;
  

  // choice of distribution and seed needs to be done here
  
  if (one4one){
    int base[2]={setup->getSeed(),rank};
    ql.init(one4one,&base[0]);
  } else {
    int base[6] = { 0,1,2,3,4,5};
    ql.init(one4one,&base[0]);
  }

  
  sn.init(setup->getSeed(),rank);


  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    slice.current=prof->value(s[i],current,currentref);
    double ne=slice.current*lambda*sample/ce;
    slice.gamma =prof->value(s[i],gamma,gammaref);
    slice.delgam=prof->value(s[i],delgam,delgamref);
    slice.ex    =prof->value(s[i],ex,exref);
    slice.ey    =prof->value(s[i],ey,eyref);
    slice.betax =prof->value(s[i],betax,betaxref);
    slice.betay =prof->value(s[i],betay,betayref);
    slice.alphax=prof->value(s[i],alphax,alphaxref);
    slice.alphay=prof->value(s[i],alphay,alphayref);
    slice.xcen  =prof->value(s[i],xcen,xcenref);
    slice.ycen  =prof->value(s[i],ycen,ycenref);
    slice.pxcen =prof->value(s[i],pxcen,pxcenref);
    slice.pycen =prof->value(s[i],pycen,pycenref);
    slice.bunch =prof->value(s[i],bunch,bunchref);
    slice.emod  =prof->value(s[i],emod,emodref);
    slice.bunchphase =prof->value(s[i],bunchphase,bunchphaseref);
    slice.emodphase  =prof->value(s[i],emodphase,emodphaseref);
    
    int npartloc=npart;
    if (one4one) { npartloc=static_cast<int>(round(ne)); }
 
    if (npartloc>nbeam){
        delete [] beamslice;
        nbeam=npartloc;
        beamslice = new Particle [nbeam];
    }

    ql.loadQuiet(beamslice, &slice, npartloc, nbins,theta0,i);
    if ((shotnoise)&&(!one4one)&&(dotime)){ sn.applyShotNoise(beamslice,npartloc,nbins,ne); }

    beam->beam[j].resize(npartloc);
    beam->current[j]=slice.current;
    for (int k=0;k<npartloc;k++){
      beam->beam[j].at(k).gamma=beamslice[k].gamma;
      beam->beam[j].at(k).theta=beamslice[k].theta;
      beam->beam[j].at(k).x    =beamslice[k].x;
      beam->beam[j].at(k).y    =beamslice[k].y;
      beam->beam[j].at(k).px   =beamslice[k].px;
      beam->beam[j].at(k).py   =beamslice[k].py;
      }
  }



  delete [] beamslice;  



 
  return true;

}

::::::::::::::
src/Loading/LoadField.cpp
::::::::::::::
#include "LoadField.h"

LoadField::LoadField()
{
  lambda=0;
  lambdaref="";
  power=0;
  powerref="";
  phase=0;
  phaseref="";
  z0=0;
  z0ref="";
  w0=100e-6;
  w0ref="";
  dgrid=1e-3;
  ngrid=151;
  harm=1;
  nx=0;
  ny=0;
  xcen=0;
  ycen=0;
  xangle=0;
  yangle=0;
  add=false;
}

LoadField::~LoadField(){}



void LoadField::usage(){

  cout << "List of keywords for FIELD" << endl;
  cout << "&field" << endl;
  cout << " double lambda = lambdaref" << endl;
  cout << " double power = 0 / reference" << endl;
  cout << " double phase = 0 / reference" << endl;
  cout << " double waist_pos = 0 / reference" << endl;
  cout << " double waist_size = 100e-6 / reference" << endl;
  cout << " double xcenter = 0" << endl;
  cout << " double ycenter = 0" << endl;
  cout << " double xangle = 0" << endl;
  cout << " double yangle = 0" << endl;
  cout << " double dgrid = 1e-3 / from existing field" << endl;
  cout << " int ngrid = 151 / from existing field" << endl;
  cout << " int harm = 1" << endl;
  cout << " int nx = 0" << endl;
  cout << " int ny = 0" << endl;
  cout << " bool accumulate = false " << endl;
  cout << "&end" << endl << endl;
  return;
}

bool LoadField::init(int rank, int size, map<string,string> *arg, vector<Field *> *fieldin,  Setup *setup, Time *time, Profile *prof)
{

  bool dotime=time->isTime();                  // check for time simulation
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  lambda=setup->getReferenceLength();

  map<string,string>::iterator end=arg->end();

  if (arg->find("lambda")!=end){this->reference(arg->at("lambda"),&lambda,&lambdaref); arg->erase(arg->find("lambda"));}
  if (arg->find("power")!=end) {this->reference(arg->at("power"),&power,&powerref); arg->erase(arg->find("power"));}
  if (arg->find("phase")!=end) {this->reference(arg->at("phase"),&phase,&phaseref); arg->erase(arg->find("phase"));}
  if (arg->find("waist_pos")!=end) {this->reference(arg->at("waist_pos"),&z0,&z0ref); arg->erase(arg->find("waist_pos"));}
  if (arg->find("waist_size")!=end){this->reference(arg->at("waist_size"),&w0,&w0ref); arg->erase(arg->find("waist_size"));}  
  if (arg->find("xcenter")!=end) {xcen = atof(arg->at("xcenter").c_str()); arg->erase(arg->find("xcenter"));}
  if (arg->find("ycenter")!=end) {ycen = atof(arg->at("ycenter").c_str()); arg->erase(arg->find("ycenter"));}
  if (arg->find("xangle")!=end)  {xangle = atof(arg->at("xangle").c_str()); arg->erase(arg->find("xangle"));}
  if (arg->find("yangle")!=end)  {yangle = atof(arg->at("yangle").c_str()); arg->erase(arg->find("yangle"));}
  if (arg->find("dgrid")!=end) {dgrid = atof(arg->at("dgrid").c_str()); arg->erase(arg->find("dgrid"));}
  if (arg->find("ngrid")!=end) {ngrid = atoi(arg->at("ngrid").c_str());  arg->erase(arg->find("ngrid"));}
  if (arg->find("harm")!=end)  {harm  = atoi(arg->at("harm").c_str());  arg->erase(arg->find("harm"));}
  if (arg->find("nx")!=end)    {nx  = atoi(arg->at("nx").c_str());  arg->erase(arg->find("nx"));}
  if (arg->find("ny")!=end)    {ny  = atoi(arg->at("ny").c_str());  arg->erase(arg->find("ny"));}
  if (arg->find("accumulate")!=end) {add  = atob(arg->at("accumulate"));  arg->erase(arg->find("accumulate"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &field" << endl; this->usage();}
    return false;
  }

  // check for existing field record;
  int idx=-1;
  Field *field;
  for (int i=0; i<fieldin->size();i++){
    if (fieldin->at(i)->harm==harm){
      field=fieldin->at(i);
      idx=i;
      if (add){
         dgrid=field->gridmax;    // take grid size from existing field
         ngrid=field->ngrid; 
      }
    }
  }

  if (idx<0){
    field=new Field;
    add=false;
  }

  // field points now to the record which is filled


  if (rank==0){ 
    if (add) { cout << "Adding "; } else  { cout << "Generating "; }
    cout << "input radiation field for HARM = " << harm <<  " ..." << endl; 
  }

  vector<double> s;
  int nslice=time->getPosition(&s);

  field->init(time->getNodeNSlice(),ngrid,dgrid,lambda,sample*lambda,s[0],harm);
  
  if (idx<0){
    fieldin->push_back(field);
    idx=fieldin->size()-1;
  }

  complex< double >  *fieldslice = new complex<double> [ngrid*ngrid];
  FieldSlice slice;
  GaussHermite gh;

  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    slice.lambda=prof->value(s[i],lambda,lambdaref);
    slice.power=prof->value(s[i],power,powerref);
    slice.phase=prof->value(s[i],phase,phaseref);
    slice.z0=prof->value(s[i],z0,z0ref);
    slice.w0=prof->value(s[i],w0,w0ref);
    slice.xcen=xcen;
    slice.ycen=ycen;
    slice.xangle=xangle;
    slice.yangle=yangle;
    slice.nx=nx;
    slice.ny=ny;
    slice.harm=harm;
    gh.loadGauss(fieldslice,&slice,dgrid,ngrid);
    if (add){
      for (int k=0; k<ngrid*ngrid;k++){
        fieldin->at(idx)->field[j].at(k)+=fieldslice[k];
      } 
    } else {
      for (int k=0; k<ngrid*ngrid;k++){
        fieldin->at(idx)->field[j].at(k)=fieldslice[k]; 
      }      
    }

  }


  delete [] fieldslice;  

  // if (idx<0){fieldin->push_back(field);}

  

  return true;

}
::::::::::::::
src/Loading/Profile.cpp
::::::::::::::
#include "Profile.h"


Profile::Profile()
{
}

Profile::~Profile()
{
  prof.clear();
}

bool Profile::init(int rank, map<string,string> *arg,string element)
{

  ProfileBase *p;
  string label;

  if (element.compare("&profile_const")==0){
    p=(ProfileBase *)new ProfileConst();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_gauss")==0){
    p=(ProfileBase *)new ProfileGauss();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_polynom")==0){
    p=(ProfileBase *)new ProfilePolynom();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_step")==0){
    p=(ProfileBase *)new ProfileStep();
    label=p->init(rank,arg);
  } 
  if (element.compare("&profile_file")==0){
    p=(ProfileBase *)new ProfileFile();
    label=p->init(rank,arg);
  } 

  if (label.size()<1){
    return false;
  } else {
    prof[label]=p;
  }
  if (rank==0) {cout << "Adding profile with label: " << label << endl;}
  return true;
}

double Profile::value(double s, double val, string label)
{
  if ((label.size()<1)||(prof.find(label)==prof.end())){  
    return val;
  } else {
    return  prof[label]->value(s);
  }
}






//------------------------------------
// individual profiles


string ProfileConst::init(int rank, map<string,string>*arg)
{

  string label="";
  c0=0;
  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_const" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_const" << endl; this->usage();
  }
  return label;
}

double ProfileConst::value(double z)
{
  return c0;
}

void ProfileConst::usage(){
  cout << "List of keywords for PROFILE_CONST" << endl;
  cout << "&profile_const" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

//-----------------------

string ProfilePolynom::init(int rank, map<string,string>*arg)
{
  string label="";
  map<string,string>::iterator end=arg->end();
  c.resize(5);
  for (int i=0; i< c.size();i++){ c[i]=0;}
  

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c[0]    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("c1")!=end)   {c[1]    = atof(arg->at("c1").c_str());  arg->erase(arg->find("c1"));}
  if (arg->find("c2")!=end)   {c[2]    = atof(arg->at("c2").c_str());  arg->erase(arg->find("c2"));}
  if (arg->find("c3")!=end)   {c[3]    = atof(arg->at("c3").c_str());  arg->erase(arg->find("c3"));}
  if (arg->find("c4")!=end)   {c[4]    = atof(arg->at("c4").c_str());  arg->erase(arg->find("c4"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown element in &profile_polynom" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_polynom" << endl; this->usage();
  }
  return label;
}


double ProfilePolynom::value(double z)
{
  double val=0;
  double zsave=1;
  for (int i=0;i<c.size();i++){
    val+=c[i]*zsave;
    zsave*=z;
  }
  return val;
}

void ProfilePolynom::usage(){
  cout << "List of keywords for PROFILE_POLYNOM" << endl;
  cout << "&profile_polynom" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double c1 = 0" << endl;
  cout << " double c2 = 0" << endl;
  cout << " double c3 = 0" << endl;
  cout << " double c4 = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

//-----------------------

 string ProfileStep::init(int rank, map<string,string>*arg)
{

  string label="";
  c0=0;
  sstart=0;
  send=0;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end)  {label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)     {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("s_start")!=end){sstart= atof(arg->at("s_start").c_str());  arg->erase(arg->find("s_start"));}
  if (arg->find("s_end")!=end)  {send  = atof(arg->at("s_end").c_str());  arg->erase(arg->find("s_end"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_step" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_step" << endl; this->usage();
  }
  return label;
}

double ProfileStep::value(double z)
{
  if ((z>=sstart) && (z <=send)){
    return c0;
  }
  return 0;
}

void ProfileStep::usage(){

  cout << "List of keywords for PROFILE_STEP" << endl;
  cout << "&profile_step" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double s_start = 0" << endl;
  cout << " double s_end = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

//-----------------------------------

string ProfileGauss::init(int rank, map<string,string>*arg)
{
  string label="";
  c0=0;
  s0=0;
  sig=1;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("s0")!=end)   {s0    = atof(arg->at("s0").c_str());  arg->erase(arg->find("s0"));}
  if (arg->find("sig")!=end)  {sig   = atof(arg->at("sig").c_str());  arg->erase(arg->find("sig"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_gauss" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_gauss" << endl; this->usage();
  }
  return label;
}

double ProfileGauss::value(double z)
{
  return c0*exp(-0.5*(z-s0)*(z-s0)/sig/sig);
}

void ProfileGauss::usage(){ 
  cout << "List of keywords for PROFILE_GAUSS" << endl;
  cout << "&profile_gauss" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double s0 = 0" << endl;
  cout << " double sig= 1" << endl;
  cout << "&end" << endl << endl;
  return;
  return;
}


//-----------------------------------

string ProfileFile::init(int rank, map<string,string>*arg)
{
  string label="";
  xdataset="";
  ydataset="";
  isTime=false;
  revert=false;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("xdata")!=end)   {xdataset = arg->at("xdata");  arg->erase(arg->find("xdata"));}
  if (arg->find("ydata")!=end)   {ydataset = arg->at("ydata");  arg->erase(arg->find("ydata"));}
  if (arg->find("isTime")!=end)  {isTime = atob(arg->at("isTime").c_str()); arg->erase(arg->find("isTime"));}
  if (arg->find("reverse")!=end) {revert = atob(arg->at("reverse").c_str()); arg->erase(arg->find("reverse"));}
  

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_file" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_file" << endl; this->usage();
  }

  int ndata=-1;
  bool success;

  success=this->simpleReadDouble1D(xdataset,&xdat);
  if (!success){
    if (rank==0){
      cout << "*** Error: Cannot read the HDF5 dataset: " << xdataset << endl;
    }
    return "";
  }

  success=this->simpleReadDouble1D(ydataset,&ydat);
  if (!success){
    if (rank==0){
      cout << "*** Error: Cannot read the HDF5 dataset: " << ydataset << endl;
    }
    return "";
  }


  if (isTime){ 
    for (int i=0; i<xdat.size();i++){
      xdat[i]*=3e8;         // scale time variable to space varial by multiplying the speed of light
    }  
  }
  
  if (revert){
    double xmin=xdat[0];
    double xmax=xdat[xdat.size()-1];
    reverse(xdat.begin(),xdat.end());
    reverse(ydat.begin(),ydat.end());
    for (int i=0;i<xdat.size();i++){
      xdat[i]=-xdat[i]+xmin+xmax;    // get the correct time window
    }
  }
  return label;
}

double ProfileFile::value(double z)
{
  if (z<xdat[0]){ return ydat[0]; }
  if (z>xdat[xdat.size()-1]){ return ydat[xdat.size()-1]; }
  int idx=0;
  while(z>=xdat[idx]){
    idx++;
  }
  idx--;
  double wei=(z-xdat[idx])/(xdat[idx+1]-xdat[idx]);
  double val=ydat[idx]*(1-wei)+wei*ydat[idx+1];
  return val;
}

void ProfileFile::usage(){ 
  cout << "List of keywords for PROFILE_FILE" << endl;
  cout << "&profile_file" << endl;
  cout << " string label = <empty>" << endl;
  cout << " string xdata = <empty>" << endl;
  cout << " string ydata = <empty>" << endl;
  cout << " bool isTime  = false" << endl;
  cout << " bool reverse = false" << endl;
  cout << "&end" << endl << endl;
  return;
  return;
}
::::::::::::::
src/Loading/QuietLoading.cpp
::::::::::::::
#include "QuietLoading.h"

QuietLoading::QuietLoading(){
  sx=NULL;
  sy=NULL;
  st=NULL;
  spx=NULL;
  spy=NULL;
  sg=NULL;
}
QuietLoading::~QuietLoading(){}



void QuietLoading::init(bool one4one, int *base)
{
  if (sx !=NULL) { delete sx; }
  if (sy !=NULL) { delete sy; }
  if (spx!=NULL) { delete spx; }
  if (spy!=NULL) { delete spy; }
  if (st !=NULL) { delete st; }
  if (sg !=NULL) { delete sg; }

  if (one4one){
     RandomU rseed(base[0]);
     double val;
     for (int i=0; i<=base[1];i++){
        val=rseed.getElement();
     }
     val*=1e9;
     int locseed=static_cast<int> (round(val));
     st  = (Sequence *) new RandomU (locseed);
     sg  = st;
     sx  = st;
     sy  = st;
     spx = st;
     spy = st;

  } else {
    st  = (Sequence *) new Hammerslay (base[0]);
    sg  = (Sequence *) new Hammerslay (base[1]);
    sx  = (Sequence *) new Hammerslay (base[2]);
    sy  = (Sequence *) new Hammerslay (base[3]);
    spx = (Sequence *) new Hammerslay (base[4]);
    spy = (Sequence *) new Hammerslay (base[5]);
  }

}

void QuietLoading::loadQuiet(Particle *beam, BeamSlice *slice, int npart, int nbins, double theta0, int islice)
{


  // resets Hammersley sequence but does nothing for random sequence; 
  Sequence *seed = (Sequence *) new RandomU(islice);
  int iseed=static_cast<int>(round(seed->getElement()*1e9));

  // initialize the sequence to new values to avoid that all core shave the same distribution
  st->set(iseed);
  sg->set(iseed);
  sx->set(iseed);
  sy->set(iseed);
  spx->set(iseed);
  spy->set(iseed);

  int mpart=npart/nbins; 
  
  Inverfc erf;

  double dtheta=1./static_cast<double>(nbins);  
  // raw distribution  
   for (int i=0;i <mpart; i++){
    beam[i].theta=st->getElement()*dtheta;
    beam[i].gamma=erf.value(2*sg->getElement());
    beam[i].x    =erf.value(2*sx->getElement());
    beam[i].y    =erf.value(2*sy->getElement());
    beam[i].px   =erf.value(2*spx->getElement());
    beam[i].py   =erf.value(2*spy->getElement());
  }


  // normalization

  double z=0;
  double zz=0;

  double norm=0;
  if (mpart>0){
    norm=1./static_cast<double>(mpart);
  } 

  // energy

  for (int i=0; i<mpart; i++){
    z +=beam[i].gamma;
    zz+=beam[i].gamma*beam[i].gamma;
  }

  z*=norm;
  zz=sqrt(fabs(zz*norm-z*z));

  if (zz>0){ zz=1/zz;}
  zz*=slice->delgam;

  for (int i=0; i<mpart; i++){
    beam[i].gamma= (beam[i].gamma-z)*zz+slice->gamma; 
  }
  
  // x and px correlation;
  z=0;
  zz=0;
  double p=0;
  double zp=0;

  for (int i=0; i<mpart;i++){
    z +=beam[i].x;
    zz+=beam[i].x*beam[i].x;
    p +=beam[i].px;
    zp+=beam[i].x*beam[i].px;    
  }
  z*=norm;
  p*=norm;
  zz=sqrt(fabs(zz*norm-z*z));
  if (zz>0) {zz=1./zz;}
  zp=(zp*norm-z*p)*zz*zz;
  
  for (int i=0; i<mpart;i++){
    beam[i].px-=zp*beam[i].x;
    beam[i].x=(beam[i].x-z)*zz;
  }

  // y and py correlation;
  z=0;
  zz=0;
  p=0;
  zp=0;

  for (int i=0; i<mpart;i++){
    z +=beam[i].y;
    zz+=beam[i].y*beam[i].y;
    p +=beam[i].py;
    zp+=beam[i].y*beam[i].py;    
  }
  z*=norm;
  p*=norm;
  zz=sqrt(fabs(zz*norm-z*z));
  if (zz>0) {zz=1./zz;}
  zp=(zp*norm-z*p)*zz*zz;
  
  for (int i=0; i<mpart;i++){
    beam[i].py-=zp*beam[i].y;
    beam[i].y=(beam[i].y-z)*zz;
  }


  // px and py size

  z=0;
  zz=0;
  p=0;
  double pp=0;

  for (int i=0; i<mpart;i++){
    z +=beam[i].px;
    zz+=beam[i].px*beam[i].px;
    p +=beam[i].py;
    pp+=beam[i].py*beam[i].py;    
  }
  z*=norm;
  p*=norm;
  zz=sqrt(fabs(zz*norm-z*z));
  if (zz>0) {zz=1./zz;}
  pp=sqrt(fabs(pp*norm-p*p));
  if (pp>0) {pp=1./pp;}
  
  for (int i=0; i<mpart;i++){
    beam[i].px=(beam[i].px-z)*zz;
    beam[i].py=(beam[i].py-p)*pp;
  }
   
  // scale to physical size
 
  double sigx=sqrt(slice->ex*slice->betax/slice->gamma);
  double sigy=sqrt(slice->ey*slice->betay/slice->gamma);
  double sigpx=sqrt(slice->ex/slice->betax/slice->gamma);
  double sigpy=sqrt(slice->ey/slice->betay/slice->gamma);
  double corx=-slice->alphax/slice->betax;
  double cory=-slice->alphay/slice->betay;


  for (int i=0; i<mpart;i++){
    beam[i].x *=sigx;
    beam[i].y *=sigy;
    beam[i].px=sigpx*beam[i].px+corx*beam[i].x;
    beam[i].py=sigpy*beam[i].py+cory*beam[i].y;
    beam[i].px*=beam[i].gamma;
    beam[i].py*=beam[i].gamma;
    beam[i].x +=slice->xcen;
    beam[i].y +=slice->ycen;
    beam[i].px+=slice->pxcen;
    beam[i].py+=slice->pycen;
  }

  // mirror particles for quiet loading

  for (int i=mpart; i>0; i--){
    int i1=i-1;
    int i2=nbins*i1;
    for (int j=0;j<nbins;j++){
      beam[i2+j].gamma=beam[i1].gamma;
      beam[i2+j].x    =beam[i1].x;
      beam[i2+j].y    =beam[i1].y;
      beam[i2+j].px   =beam[i1].px;
      beam[i2+j].py   =beam[i1].py;
      beam[i2+j].theta=beam[i1].theta+j*dtheta;     
    } 
  }


  // scale phase

  for (int i=0;i<npart;i++){
    beam[i].theta*=theta0;
  }

  // bunching and energy modulation

  if ((slice->bunch!=0)||(slice->emod!=0)){
    for (int i=0; i<npart;i++){
      beam[i].gamma-=slice->emod*sin(beam[i].theta-slice->emodphase);
      beam[i].theta-=2*slice->bunch*sin(beam[i].theta-slice->bunchphase);
    }
  }

  return;
}
::::::::::::::
src/Loading/ShotNoise.cpp
::::::::::::::
#include "ShotNoise.h"

ShotNoise::ShotNoise(){
  sran=NULL;
  nwork=1000;
  work=new double [nwork];
}

ShotNoise::~ShotNoise()
{
  delete [] work;
}



void ShotNoise::init(int base,int rank)
{
  RandomU rseed(base);
  double val;
  for (int i=0; i<=rank;i++){
    val=rseed.getElement();
  }
  val*=1e9;
  int locseed=static_cast<int> (round(val));
  if (sran !=NULL) { delete sran; }
  sran  = new RandomU (locseed);
}

void ShotNoise::applyShotNoise(Particle *beam, int npart, int nbins, double ne)
{


  if (npart>nwork){
    delete [] work;
    nwork=npart;
    work = new double [nwork];
  }

  int mpart=npart/nbins;

  double nbl=ne/static_cast<double>(mpart);  // number of simulated electrons per beamlet 

  for (int i=0; i< npart; i++){
    work[i]=0;
  }
  for (int ih=0; ih< (nbins-1)/2; ih++){
    for (int i1=0;i1<mpart;i1++){
      double phi=4*asin(1)*sran->getElement();
      double an=sqrt(-log(sran->getElement())/nbl)*2/static_cast<double>(ih+1);
      for (int i2=0; i2<nbins; i2++){
        int idx=i1*nbins+i2;
        work[idx]-=an*sin(beam[idx].theta*(ih+1)+phi);
      }
    }
  }
  
  for (int i=0; i< npart; i++){
    beam[i].theta+=work[i];
  }

  return;
}
::::::::::::::
src/Main/AlterSetup.cpp
::::::::::::::
#include "AlterSetup.h"


AlterSetup::AlterSetup()
{
  rootname="";
  lattice="";
  beamline="";
  harmonic=1;
  subharmonic=1;
  resample=false;
  disable=false;

}

AlterSetup::~AlterSetup(){}

void AlterSetup::usage(){

  cout << "List of keywords for ALTER_SETUP" << endl;
  cout << "&alter_setup" << endl;
  cout << " string rootname = <taken from SETUP + Increment>" << endl;
  //  cout << " string lattice = <taken from SETUP>" << endl;
  cout << " string beamline = <empty>" << endl;
  cout << " double delz = <taken from SETUP>" << endl;
  cout << " int harmonic = 1" << endl;
  cout << " int subharmonic = 1" << endl;
  cout << " bool resample = false" << endl;
  cout << " bool disable = false" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool AlterSetup::init(int inrank, map<string,string> *arg, Setup *setup, Lattice *lat, Time *time, Beam *beam, vector<Field *> *field)
{

  beamline="";  // clear beamline name for multiple calls of altersetup
  lattice=setup->getLattice();
  delz=setup->getStepLength();
  rank=inrank;
  
  map<string,string>::iterator end=arg->end();

  if (arg->find("rootname")!=end){rootname = arg->at("rootname"); arg->erase(arg->find("rootname"));}
  //  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice"); arg->erase(arg->find("lattice"));}
  if (arg->find("beamline")!=end){beamline = arg->at("beamline"); arg->erase(arg->find("beamline"));}
  if (arg->find("delz")!=end)    {delz     = atof(arg->at("delz").c_str());  arg->erase(arg->find("delz"));}
  if (arg->find("subharmonic")!=end){subharmonic  = atoi(arg->at("subharmonic").c_str());  arg->erase(arg->find("subharmonic"));}
  if (arg->find("harmonic")!=end){harmonic  = atoi(arg->at("harmonic").c_str());  arg->erase(arg->find("harmonic"));}
  if (arg->find("resample")!=end){resample  = atob(arg->at("resample"));  arg->erase(arg->find("resample"));}
  if (arg->find("disable")!=end){disable  = atob(arg->at("disable"));  arg->erase(arg->find("disable"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &alter_setup" << endl; this->usage();}
    return false;
  }

  setup->setStepLength(delz);
  if (!time->isTime()){ resample=false; }  // no resampling in non-timedependent runs allowed

  // step one: Select new lattice if selected
  if (beamline!="") {
     bool status = lat->parse(lattice,beamline,rank);
     if (status==false) { return status ; }
  }

  // step two: Add increment for rootname or assign new root name
  if (rootname!=""){
    setup->setRootName(&rootname); 
  } else {
    setup->incrementCount();
  }

  //step three: do subharmonic conversion
  // step 3.1 - time window adjustment
  if (subharmonic>1){
    if (rank==0){ cout << "Converting to subharmonic: " << subharmonic << "..." << endl;}
    double lam0=setup->getReferenceLength()*static_cast<double>(subharmonic);
    setup->setReferenceLength(lam0);
    if (!resample) {
      double samp=time->getSampleRate()/static_cast<double>(subharmonic);
      time->setSampleRate(samp);
    }
    time->finishInit(setup);
    // step 3.2 - beam
    if (!beam->subharmonicConversion(subharmonic,resample)){ 
      if (rank==0) {cout << "*** Error: Cannot convert beam distribution to lower harmonic" << endl;}
      return false;
    }
    // step 3.3 - field
    for (int i=0; i<field->size();i++){
      if (field->at(i)->getHarm()!=1){
	if (rank==0) {cout << "Deleting higher radiation harmonic: " << field->at(i)->getHarm() << " in subharmonic conversion" << endl;}
	field->erase(field->begin()+i);
	i--; // reseting the count
      }else{
	if (rank==0) {cout << "Converting radiation fundamental to harmonic: " << subharmonic << " and keep it in memory" << endl;}
	 if (!(field->at(i)->subharmonicConversion(subharmonic,resample))){
	    if (rank==0) {cout << "*** Error: Cannot convert field distribution to higher harmonic" << endl;}
	    return false;	    
	 }
      }
    }
  }

  // step four: do harmonic conversion
  // step 4.1 - time window
  if (harmonic>1){
    if (rank==0){ cout << "Converting to harmonic: " << harmonic << "..." << endl;}
    double lam0=setup->getReferenceLength()/static_cast<double>(harmonic);
    setup->setReferenceLength(lam0);
    if (!resample) {
      double sam=time->getSampleRate()*static_cast<double>(harmonic);
      time->setSampleRate(sam);
    }
    time->finishInit(setup);

    // step 4.2 - beam
    if (!beam->harmonicConversion(harmonic,resample)){ 
      if (rank==0) {cout << "*** Error: Cannot convert beam distribution to higher harmonic" << endl;}
      return false;
    }
    
    // step 4.3 - field
    for (int i=0; i<field->size();i++){
      if (field->at(i)->getHarm()!=harmonic){
        if (disable) {
	    if (rank==0) {cout << "Disabling non-matching radiation harmonic: " << field->at(i)->getHarm() << " in harmonic conversion" << endl;}
	    field->at(i)->disable(1./static_cast<double>(harmonic));

	} else {
	    if (rank==0) {cout << "Deleting non-matching radiation harmonic: " << field->at(i)->getHarm() << " in harmonic conversion" << endl;}
	    field->erase(field->begin()+i);
	    i--; // reseting the count
        }
      }else{
	if (rank==0) {cout << "Converting radiation harmonic: " << harmonic << " to fundamental and keep it in memory" << endl;}
	 if (!(field->at(i)->harmonicConversion(harmonic,resample))){
	    if (rank==0) {cout << "*** Error: Cannot convert field distribution to higher harmonic" << endl;}
	    return false;	    
	 }
      }
  
    }
  }

  return true;


}


::::::::::::::
src/Main/Dump.cpp
::::::::::::::
#include "dump.h"
#include "writeFieldHDF5.h"
#include "writeBeamHDF5.h"

Dump::Dump()
{
  dumpfield="";
  dumpbeam ="";
}

Dump::~Dump(){}

void Dump::usage(){

  cout << "List of keywords for WRITE" << endl;
  cout << "&write" << endl;
  cout << " string field = <empty>" << endl;
  cout << " string beam  = <empty>" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool Dump::init(int inrank, int insize, map<string,string> *arg, Setup *setup, Beam *beam, vector<Field *> *field)
{

  rank=inrank;
  size=insize;

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("field")!=end){dumpfield=arg->at("field"); arg->erase(arg->find("field"));}
  if (arg->find("beam")!=end) {dumpbeam =arg->at("beam");  arg->erase(arg->find("beam"));}
  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &write" << endl; this->usage();}
    return false;
  }

  
  if (dumpfield.size()>0){
   WriteFieldHDF5 dump;
   dump.write(dumpfield,field);
  }
  if (dumpbeam.size()>0){
   WriteBeamHDF5 dump;
   dump.write(dumpbeam,beam);
  }

  return true;


}


::::::::::::::
src/Main/EField.cpp
::::::::::::::
#include "EField.h"
#include "Beam.h"

EField::EField(){}
EField::~EField(){}

void EField::usage(){

  cout << "List of keywords for EFIELD" << endl;
  cout << "&efield" << endl;
  cout << " double rmax = 0 " << endl;
  cout << " int nz   = 0" << endl;
  cout << " int nphi = 0" << endl;
  cout << " int ngrid = 100" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool EField::init(int rank, int size, map<string,string> *arg,  Beam *beam, Setup *setup, Time *time)
{

  bool dotime=time->isTime();                  // check for time simulation
  lambda=setup->getReferenceLength();

  map<string,string>::iterator end=arg->end();

  if (arg->find("rmax")!=end)  {rmax  = atof(arg->at("rmax").c_str()); arg->erase(arg->find("rmax"));}
  if (arg->find("ngrid")!=end) {ngrid = atoi(arg->at("ngrid").c_str());  arg->erase(arg->find("ngrid"));}
  if (arg->find("nz")!=end)    {nz    = atoi(arg->at("nz").c_str());  arg->erase(arg->find("nz"));}
  if (arg->find("nphi")!=end)  {nphi  = atoi(arg->at("nphi").c_str());  arg->erase(arg->find("nphi"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &efield" << endl; this->usage();}
    return false;
  }


  beam->initEField(rmax,ngrid,nz,nphi,lambda);


  return true;
}
 

 
::::::::::::::
src/Main/GenMain.cpp
::::::::::::::
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <ctime>


#include <fenv.h>
#include <signal.h>

#include "mpi.h"



// genesis headerfiles & classes


#include "Beam.h"
#include "Field.h"
#include "EField.h"

#include "Parser.h"
#include "Profile.h"
#include "Setup.h"
#include "AlterSetup.h"
#include "Lattice.h"
#include "Time.h"
#include "Gencore.h"
#include "LoadField.h"
#include "LoadBeam.h"
#include "AlterLattice.h"
#include "Track.h"
#include "SDDSBeam.h"
#include "SponRad.h"
#include "dump.h"
#include "ImportBeam.h"
#include "ImportField.h"
#include "writeBeamHDF5.h"
#include "writeFieldHDF5.h"


#include "Collective.h"

using namespace std;

const double vacimp = 376.73;
const double eev    = 510999.06; 
const double ce     = 4.8032045e-11;


const int versionmajor = 4;
const int versionminor = 0;
const int versionrevision = 5;
const bool versionbeta=true;

string *meta_inputfile;
string *meta_latfile;

bool MPISingle;  // global variable to do mpi or not

int genmain (char *inputfile, bool split) {

 
  MPISingle=split;       

        int rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node
        int size=MPI::COMM_WORLD.Get_size(); // get size of cluster
	if (MPISingle){
	  rank=0;
	  size=1;
        }


        time_t timer;
	if (rank==0) {
          time(&timer);
          cout << "---------------------------------------------" << endl;
          cout << "GENESIS - Version " <<  versionmajor <<"."<< versionminor << "." << versionrevision ;
	  if (versionbeta) {cout << " (beta)";}
	  cout << " has started..." << endl;			
	  cout << "Starting Time: " << ctime(&timer)<< endl;
          cout << "MPI-Comm Size: " << size << " nodes" << endl << endl;
        }

        //---------------------------------------------------------
        // Instance of beam and field class to carry the distribution

        vector<Field *> field;   // an vector of various field components (harmonics, vertical/horizonthal components)
        Beam  *beam =new Beam;


        //----------------------------------------------------------
        // main loop extracting one element with arguments at a time
      
        Parser parser; 
        string element;
        map<string,string> argument;

        Setup *setup=new Setup;
	AlterLattice *alt=new AlterLattice;
        Lattice *lattice=new Lattice;
        Profile *profile=new Profile;
        Time *timewindow=new Time;

	meta_inputfile=new string (inputfile);

        parser.open(inputfile,rank);

        while(parser.parse(&element,&argument)){
           
          //----------------------------------------------
	  // setup & parsing the lattice file

          if (element.compare("&setup")==0){
            if (!setup->init(rank,&argument,lattice)){ break;}
	    meta_latfile=new string (setup->getLattice());
            continue;  
          }  

          //----------------------------------------------
	  // modifying run

          if (element.compare("&alter_setup")==0){
	    AlterSetup *altersetup= new AlterSetup;
            if (!altersetup->init(rank,&argument,setup,lattice,timewindow,beam,&field)){ break;}
	    delete altersetup;
            continue;  
          }  

          //----------------------------------------------
	  // modifying the lattice file

          if (element.compare("&lattice")==0){
            if (!alt->init(rank,size,&argument,lattice,setup)){ break;}
            continue;  
          }  

          //---------------------------------------------------
          // adding profile elements

          if ((element.compare("&profile_const")==0)||
              (element.compare("&profile_gauss")==0)||
              (element.compare("&profile_file")==0)||
              (element.compare("&profile_polynom")==0)||
              (element.compare("&profile_step")==0)){            
            if (!profile->init(rank,&argument,element)){ break; }
            continue;
	  }

          //----------------------------------------------------
          // defining the time window of simulation

	  if (element.compare("&time")==0){
            if (!timewindow->init(rank,size,&argument,setup)){ break;}
            continue;  
          }  

          //----------------------------------------------------
          // internal generation of the field

	  if (element.compare("&field")==0){
	    LoadField *loadfield=new LoadField;
            if (!loadfield->init(rank,size,&argument,&field,setup,timewindow,profile)){ break;}
	    delete loadfield;
            continue;  
          }  
	 
          //----------------------------------------------------
          // setup of space charge field

	  if (element.compare("&efield")==0){
   	    EField *efield=new EField;
            if (!efield->init(rank,size,&argument,beam,setup,timewindow)){ break;}
	    delete efield;
            continue;  
          }  

          //----------------------------------------------------
          // setup of spontaneous radiation

	  if (element.compare("&sponrad")==0){
            SponRad *sponrad=new SponRad;
            if (!sponrad->init(rank,size,&argument,beam)){ break;}
	    delete sponrad;
            continue;  
          }  

          //----------------------------------------------------
          // internal generation of beam

	  if (element.compare("&beam")==0){
            LoadBeam *loadbeam=new LoadBeam;
            if (!loadbeam->init(rank,size,&argument,beam,setup,timewindow,profile,lattice)){ break;}
	    delete loadbeam;
            continue;  
          }  

          //----------------------------------------------------
          // external generation of beam with an sdds file

	  if (element.compare("&importdistribution")==0){
            SDDSBeam *sddsbeam=new SDDSBeam;
            if (!sddsbeam->init(rank,size,&argument,beam,setup,timewindow,lattice)){ break;}
	    delete sddsbeam;
            continue;  
          }  

          //----------------------------------------------------
          // tracking - the very core part of Genesis

	  if (element.compare("&track")==0){
            Track *track=new Track;
	    if (!track->init(rank,size,&argument,beam,&field,setup,lattice,alt,timewindow)){ break;}
            delete track;
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&sort")==0){
	    beam->sort();
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&write")==0){
            Dump *dump=new Dump;
	    if (!dump->init(rank,size,&argument,setup,beam,&field)){ break;}
            delete dump;
            continue;  
          }  


          //----------------------------------------------------
          // import beam from a particle dump

	  if (element.compare("&importbeam")==0){
            ImportBeam *import=new ImportBeam;
	    if (!import->init(rank,size,&argument,beam,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  


          //----------------------------------------------------
          // import field from a field dump

	  if (element.compare("&importfield")==0){
            ImportField *import=new ImportField;
	    if (!import->init(rank,size,&argument,&field,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  



          //-----------------------------------------------------
          // error because the element typ is not defined

          if (rank==0){
            cout << "*** Error: Unknow element in input file: " << element << endl; 
	  }
          break;
        } 



 	if (rank==0) {
          time(&timer);
          cout << endl<< "Program is terminating..." << endl;
	  cout << "Ending Time: " << ctime(&timer);
          cout << "-------------------------------------" << endl;

        }



        return 0;

}
::::::::::::::
src/Main/maincore.cpp
::::::::::::::
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <ctime>


#include <fenv.h>
#include <signal.h>

#include "mpi.h"



// genesis headerfiles & classes


#include "Beam.h"
#include "Field.h"
#include "EField.h"

#include "Parser.h"
#include "Profile.h"
#include "Setup.h"
#include "AlterSetup.h"
#include "Lattice.h"
#include "Time.h"
#include "Gencore.h"
#include "LoadField.h"
#include "LoadBeam.h"
#include "AlterLattice.h"
#include "Track.h"
#include "SDDSBeam.h"
#include "SponRad.h"
#include "dump.h"
#include "ImportBeam.h"
#include "ImportField.h"
#include "writeBeamHDF5.h"
#include "writeFieldHDF5.h"


#include "Collective.h"

using namespace std;

const double vacimp = 376.73;
const double eev    = 510999.06; 
const double ce     = 4.8032045e-11;


const int versionmajor = 4;
const int versionminor = 0;
const int versionrevision = 5;
const bool versionbeta=true;

string *meta_inputfile;
string *meta_latfile;


int main (int argc, char *argv[]) {


        //-------------------------------------------------------
        // init MPI and get size etc.
        //

        MPI::Status status; //MPI
        MPI::Init(argc, argv); //MPI

       
        int size=MPI::COMM_WORLD.Get_size(); // get size of cluster
        int rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node


	Collective col;
	//        col.WakeRes();
	//	col.WakeGeo();

        time_t timer;
	if (rank==0) {
          time(&timer);
          cout << "---------------------------------------------" << endl;
          cout << "GENESIS - Version " <<  versionmajor <<"."<< versionminor << "." << versionrevision ;
	  if (versionbeta) {cout << " (beta)";}
	  cout << " has started..." << endl;			
	  cout << "Starting Time: " << ctime(&timer)<< endl;
          cout << "MPI-Comm Size: " << size << " nodes" << endl << endl;
        }


        //---------------------------------------------------------
        // Instance of beam and field class to carry the distribution

        vector<Field *> field;   // an vector of various field components (harmonics, vertical/horizonthal components)
        Beam  *beam =new Beam;


        //----------------------------------------------------------
        // main loop extracting one element with arguments at a time
      
        Parser parser; 
        string element;
        map<string,string> argument;

        Setup *setup=new Setup;
	AlterLattice *alt=new AlterLattice;
        Lattice *lattice=new Lattice;
        Profile *profile=new Profile;
        Time *timewindow=new Time;

	meta_inputfile=new string (argv[argc-1]);

        parser.open(argv[argc-1],rank);

        while(parser.parse(&element,&argument)){
           
          //----------------------------------------------
	  // setup & parsing the lattice file

          if (element.compare("&setup")==0){
            if (!setup->init(rank,&argument,lattice)){ break;}
	    meta_latfile=new string (setup->getLattice());
            continue;  
          }  

          //----------------------------------------------
	  // modifying run

          if (element.compare("&alter_setup")==0){
	    AlterSetup *altersetup= new AlterSetup;
            if (!altersetup->init(rank,&argument,setup,lattice,timewindow,beam,&field)){ break;}
	    delete altersetup;
            continue;  
          }  

          //----------------------------------------------
	  // modifying the lattice file

          if (element.compare("&lattice")==0){
            if (!alt->init(rank,size,&argument,lattice,setup)){ break;}
            continue;  
          }  

          //---------------------------------------------------
          // adding profile elements

          if ((element.compare("&profile_const")==0)||
              (element.compare("&profile_gauss")==0)||
              (element.compare("&profile_file")==0)||
              (element.compare("&profile_polynom")==0)||
              (element.compare("&profile_step")==0)){            
            if (!profile->init(rank,&argument,element)){ break; }
            continue;
	  }

          //----------------------------------------------------
          // defining the time window of simulation

	  if (element.compare("&time")==0){
            if (!timewindow->init(rank,size,&argument,setup)){ break;}
            continue;  
          }  

          //----------------------------------------------------
          // internal generation of the field

	  if (element.compare("&field")==0){
	    LoadField *loadfield=new LoadField;
            if (!loadfield->init(rank,size,&argument,&field,setup,timewindow,profile)){ break;}
	    delete loadfield;
            continue;  
          }  
	 
          //----------------------------------------------------
          // setup of space charge field

	  if (element.compare("&efield")==0){
   	    EField *efield=new EField;
            if (!efield->init(rank,size,&argument,beam,setup,timewindow)){ break;}
	    delete efield;
            continue;  
          }  

          //----------------------------------------------------
          // setup of spontaneous radiation

	  if (element.compare("&sponrad")==0){
            SponRad *sponrad=new SponRad;
            if (!sponrad->init(rank,size,&argument,beam)){ break;}
	    delete sponrad;
            continue;  
          }  

          //----------------------------------------------------
          // internal generation of beam

	  if (element.compare("&beam")==0){
            LoadBeam *loadbeam=new LoadBeam;
            if (!loadbeam->init(rank,size,&argument,beam,setup,timewindow,profile,lattice)){ break;}
	    delete loadbeam;
            continue;  
          }  

          //----------------------------------------------------
          // external generation of beam with an sdds file

	  if (element.compare("&importdistribution")==0){
            SDDSBeam *sddsbeam=new SDDSBeam;
            if (!sddsbeam->init(rank,size,&argument,beam,setup,timewindow,lattice)){ break;}
	    delete sddsbeam;
            continue;  
          }  

          //----------------------------------------------------
          // tracking - the very core part of Genesis

	  if (element.compare("&track")==0){
            Track *track=new Track;
	    if (!track->init(rank,size,&argument,beam,&field,setup,lattice,alt,timewindow)){ break;}
            delete track;
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&sort")==0){
	    beam->sort();
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&write")==0){
            Dump *dump=new Dump;
	    if (!dump->init(rank,size,&argument,setup,beam,&field)){ break;}
            delete dump;
            continue;  
          }  


          //----------------------------------------------------
          // import beam from a particle dump

	  if (element.compare("&importbeam")==0){
            ImportBeam *import=new ImportBeam;
	    if (!import->init(rank,size,&argument,beam,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  


          //----------------------------------------------------
          // import field from a field dump

	  if (element.compare("&importfield")==0){
            ImportField *import=new ImportField;
	    if (!import->init(rank,size,&argument,&field,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  



          //-----------------------------------------------------
          // error because the element typ is not defined

          if (rank==0){
            cout << "*** Error: Unknow element in input file: " << element << endl; 
	  }
          break;
        } 



 	if (rank==0) {
          time(&timer);
          cout << endl<< "Program is terminating..." << endl;
	  cout << "Ending Time: " << ctime(&timer);
          cout << "-------------------------------------" << endl;

        }


        MPI::Finalize(); // node turned off

        return 0;

}
::::::::::::::
src/Main/main.cpp
::::::::::::::
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <ctime>


#include <fenv.h>
#include <signal.h>

#include "mpi.h"



// genesis headerfiles & classes


#include "Beam.h"
#include "Field.h"
#include "EField.h"

#include "Parser.h"
#include "Profile.h"
#include "Setup.h"
#include "AlterSetup.h"
#include "Lattice.h"
#include "Time.h"
#include "Gencore.h"
#include "LoadField.h"
#include "LoadBeam.h"
#include "AlterLattice.h"
#include "Track.h"
#include "SDDSBeam.h"
#include "SponRad.h"
#include "dump.h"
#include "ImportBeam.h"
#include "ImportField.h"
#include "writeBeamHDF5.h"
#include "writeFieldHDF5.h"


#include "Collective.h"

using namespace std;

const double vacimp = 376.73;
const double eev    = 510999.06; 
const double ce     = 4.8032045e-11;


const int versionmajor = 4;
const int versionminor = 0;
const int versionrevision = 5;
const bool versionbeta=true;

string *meta_inputfile;
string *meta_latfile;


int main (int argc, char *argv[]) {


        //-------------------------------------------------------
        // init MPI and get size etc.
        //

        MPI::Status status; //MPI
        MPI::Init(argc, argv); //MPI

       
        int size=MPI::COMM_WORLD.Get_size(); // get size of cluster
        int rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node


	Collective col;
	//        col.WakeRes();
	//	col.WakeGeo();

        time_t timer;
	if (rank==0) {
          time(&timer);
          cout << "---------------------------------------------" << endl;
          cout << "GENESIS - Version " <<  versionmajor <<"."<< versionminor << "." << versionrevision ;
	  if (versionbeta) {cout << " (beta)";}
	  cout << " has started..." << endl;			
	  cout << "Starting Time: " << ctime(&timer)<< endl;
          cout << "MPI-Comm Size: " << size << " nodes" << endl << endl;
        }


        //---------------------------------------------------------
        // Instance of beam and field class to carry the distribution

        vector<Field *> field;   // an vector of various field components (harmonics, vertical/horizonthal components)
        Beam  *beam =new Beam;


        //----------------------------------------------------------
        // main loop extracting one element with arguments at a time
      
        Parser parser; 
        string element;
        map<string,string> argument;

        Setup *setup=new Setup;
	AlterLattice *alt=new AlterLattice;
        Lattice *lattice=new Lattice;
        Profile *profile=new Profile;
        Time *timewindow=new Time;

	meta_inputfile=new string (argv[argc-1]);

        parser.open(argv[argc-1],rank);

        while(parser.parse(&element,&argument)){
           
          //----------------------------------------------
	  // setup & parsing the lattice file

          if (element.compare("&setup")==0){
            if (!setup->init(rank,&argument,lattice)){ break;}
	    meta_latfile=new string (setup->getLattice());
            continue;  
          }  

          //----------------------------------------------
	  // modifying run

          if (element.compare("&alter_setup")==0){
	    AlterSetup *altersetup= new AlterSetup;
            if (!altersetup->init(rank,&argument,setup,lattice,timewindow,beam,&field)){ break;}
	    delete altersetup;
            continue;  
          }  

          //----------------------------------------------
	  // modifying the lattice file

          if (element.compare("&lattice")==0){
            if (!alt->init(rank,size,&argument,lattice,setup)){ break;}
            continue;  
          }  

          //---------------------------------------------------
          // adding profile elements

          if ((element.compare("&profile_const")==0)||
              (element.compare("&profile_gauss")==0)||
              (element.compare("&profile_file")==0)||
              (element.compare("&profile_polynom")==0)||
              (element.compare("&profile_step")==0)){            
            if (!profile->init(rank,&argument,element)){ break; }
            continue;
	  }

          //----------------------------------------------------
          // defining the time window of simulation

	  if (element.compare("&time")==0){
            if (!timewindow->init(rank,size,&argument,setup)){ break;}
            continue;  
          }  

          //----------------------------------------------------
          // internal generation of the field

	  if (element.compare("&field")==0){
	    LoadField *loadfield=new LoadField;
            if (!loadfield->init(rank,size,&argument,&field,setup,timewindow,profile)){ break;}
	    delete loadfield;
            continue;  
          }  
	 
          //----------------------------------------------------
          // setup of space charge field

	  if (element.compare("&efield")==0){
   	    EField *efield=new EField;
            if (!efield->init(rank,size,&argument,beam,setup,timewindow)){ break;}
	    delete efield;
            continue;  
          }  

          //----------------------------------------------------
          // setup of spontaneous radiation

	  if (element.compare("&sponrad")==0){
            SponRad *sponrad=new SponRad;
            if (!sponrad->init(rank,size,&argument,beam)){ break;}
	    delete sponrad;
            continue;  
          }  

          //----------------------------------------------------
          // internal generation of beam

	  if (element.compare("&beam")==0){
            LoadBeam *loadbeam=new LoadBeam;
            if (!loadbeam->init(rank,size,&argument,beam,setup,timewindow,profile,lattice)){ break;}
	    delete loadbeam;
            continue;  
          }  

          //----------------------------------------------------
          // external generation of beam with an sdds file

	  if (element.compare("&importdistribution")==0){
            SDDSBeam *sddsbeam=new SDDSBeam;
            if (!sddsbeam->init(rank,size,&argument,beam,setup,timewindow,lattice)){ break;}
	    delete sddsbeam;
            continue;  
          }  

          //----------------------------------------------------
          // tracking - the very core part of Genesis

	  if (element.compare("&track")==0){
            Track *track=new Track;
	    if (!track->init(rank,size,&argument,beam,&field,setup,lattice,alt,timewindow)){ break;}
            delete track;
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&sort")==0){
	    beam->sort();
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&write")==0){
            Dump *dump=new Dump;
	    if (!dump->init(rank,size,&argument,setup,beam,&field)){ break;}
            delete dump;
            continue;  
          }  


          //----------------------------------------------------
          // import beam from a particle dump

	  if (element.compare("&importbeam")==0){
            ImportBeam *import=new ImportBeam;
	    if (!import->init(rank,size,&argument,beam,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  


          //----------------------------------------------------
          // import field from a field dump

	  if (element.compare("&importfield")==0){
            ImportField *import=new ImportField;
	    if (!import->init(rank,size,&argument,&field,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  



          //-----------------------------------------------------
          // error because the element typ is not defined

          if (rank==0){
            cout << "*** Error: Unknow element in input file: " << element << endl; 
	  }
          break;
        } 



 	if (rank==0) {
          time(&timer);
          cout << endl<< "Program is terminating..." << endl;
	  cout << "Ending Time: " << ctime(&timer);
          cout << "-------------------------------------" << endl;

        }


        MPI::Finalize(); // node turned off

        return 0;

}
::::::::::::::
src/Main/mainwrap.cpp
::::::::::::::


#include "mpi.h"
#include "genesis.h"

// very basic wrapper for genesis. Most of the genesis stuff is moved into genmain.


int main (int argc, char *argv[]) {


        //-------------------------------------------------------
        // init MPI and get size etc.
        //

        MPI::Status status; //MPI
        MPI::Init(argc, argv); //MPI
	
	genmain(argv[argc-1],false);

        MPI::Finalize(); // node turned off

        return 0;

}
::::::::::::::
src/Main/Parser.cpp
::::::::::::::
#include "Parser.h"

Parser::Parser()
{
}

Parser::~Parser()
{
}

bool Parser::open(string file, int inrank)
{

  rank=inrank;
  fin.open(file.c_str(),ios_base::in);

  if (!fin){
    if (rank==0) {cout << "*** Error: Cannot open main input file: " << file << endl;}
    return false;
  }
  return true;

}



bool Parser::parse(string *element, map<string,string> *argument)
{

  string instring, list;
  bool accumulate=false;
  
  
  while(getline(fin,instring,'\n')){    // read line
    this->trim(instring);

    // check for empty lines or comments
    if ((!instring.compare(0,1,"#") || instring.length() < 1)){ continue; } // skip comment and empty rows

    // check for terminating element and then return
    if ((instring.compare("&end")==0)|| (instring.compare("&END")==0) || (instring.compare("&End")==0)){
      if (!accumulate){
        if (rank==0){cout << "*** Error: Termination string outside element definition in input file" << endl;}    
	return false;
      }
      return this->fillMap(&list,argument);
    }

    // check whether element starts
    if (!instring.compare(0,1,"&")){
      if (accumulate){
        if (rank==0) {cout << "*** Error: Nested elements in main input file" << endl;}
	return false;
      }
      accumulate=true;
      for (int i=0; i<instring.size();i++){
	instring[i]=tolower(instring[i]);
      }
      *element=instring;   
      continue;

    }

    // add content if in element
    if (accumulate){
       list.append(instring);  // add all content into one string
       list.append(";");
    }
  }
  
  fin.close();
  return false;

}

bool Parser::fillMap(string *list,map<string,string> *map){

  // splits the string into a map with pairs : key and value.
  // keys are converted to lower case
  
  map->clear();
  string key,val;
 
  size_t pos,tpos;
  while ((pos=list->find_first_of(";")) !=string::npos){  // split into individual lines
     val=list->substr(0,pos);
     list->erase(0,pos+1);
     this->trim(key);
     if (val.length()>0){
       tpos=val.find_first_of("=");
       if (tpos ==string::npos){
         if (rank==0){ cout << "*** Error: Invalid format " << val << " in input file" << endl;}
         return false;
       } else{
         key=val.substr(0,tpos);
         val.erase(0,tpos+1);
         this->trim(key);
         this->trim(val);
         map->insert(pair<string,string>(key,val));
       }
     }
  }
  
  return true;
}
::::::::::::::
src/Main/Setup.cpp
::::::::::::::
#include <sstream>
#include "Setup.h"

Setup::Setup()
{
  rootname="";
  lattice="";
  beamline="";
  partfile="";
  fieldfile="";
  one4one=false;
  shotnoise=true;
  nbins=4;
  npart=8192;
  gamma0=5800/0.511;
  lambda0=1e-10;
  delz=0.015; 
  seed=123456789;
  runcount = 0 ;  // count of runs in conjunction of calls of altersetup 
}

Setup::~Setup(){}

void Setup::usage(){

  cout << "List of keywords for SETUP" << endl;
  cout << "&setup" << endl;
  cout << " string rootname = <empty>" << endl;
  cout << " string lattice = <empty>" << endl;
  cout << " string beamline = <empty>" << endl;
  cout << " string partfile = <empty>" << endl;
  cout << " string fieldfile = <empty>" << endl;
  cout << " double gamma0 = 5800/0.511" << endl;
  cout << " double lambda0 = 1e-10" << endl;
  cout << " double delz = 0.015" << endl;
  cout << " int seed = 123456789" << endl;
  cout << " int npart = 8192" << endl;
  cout << " int nbins = 4" << endl;
  cout << " bool one4one = false" << endl;
  cout << " bool shotnoise = true" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool Setup::init(int inrank, map<string,string> *arg, Lattice *lat)
{

  rank=inrank;
  map<string,string>::iterator end=arg->end();

  if (arg->find("rootname")!=end){rootname = arg->at("rootname"); arg->erase(arg->find("rootname"));}
  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice");  arg->erase(arg->find("lattice"));}
  if (arg->find("beamline")!=end){beamline = arg->at("beamline"); arg->erase(arg->find("beamline"));}
  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice");  arg->erase(arg->find("lattice"));}
  if (arg->find("partfile")!=end){partfile   = arg->at("partfile");  arg->erase(arg->find("partfile"));}
  if (arg->find("fieldfile")!=end){fieldfile = arg->at("fieldfile");  arg->erase(arg->find("fieldfile"));}
  if (arg->find("gamma0")!=end)  {gamma0   = atof(arg->at("gamma0").c_str());  arg->erase(arg->find("gamma0"));}
  if (arg->find("lambda0")!=end) {lambda0  = atof(arg->at("lambda0").c_str()); arg->erase(arg->find("lambda0"));}
  if (arg->find("delz")!=end)    {delz     = atof(arg->at("delz").c_str());  arg->erase(arg->find("delz"));}
  if (arg->find("seed")!=end)    {seed     = atoi(arg->at("seed").c_str());  arg->erase(arg->find("seed"));}
  if (arg->find("one4one")!=end) {one4one  = atob(arg->at("one4one"));  arg->erase(arg->find("one4one"));}
  if (arg->find("npart")!=end)    {npart  = atoi(arg->at("npart").c_str());  arg->erase(arg->find("npart"));}
  if (arg->find("nbins")!=end)    {nbins  = atoi(arg->at("nbins").c_str());  arg->erase(arg->find("nbins"));}
  if (arg->find("shotnoise")!=end){shotnoise  = atob(arg->at("shotnoise"));  arg->erase(arg->find("shotnoise"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &setup" << endl; this->usage();}
    return false;
  }

  return lat->parse(lattice,beamline,rank);

}



bool Setup::getRootName(string *filename)
{
  if (rootname.size()<1){
    return false;
  }
  *filename=rootname;
  if (runcount>0) {
    stringstream ss;
    ss << ".Run" << (runcount+1) ;
    *filename+=ss.str();
  }
  return true; 


}
::::::::::::::
src/Main/SponRad.cpp
::::::::::::::
#include "SponRad.h"
#include "Beam.h"

SponRad::SponRad(){}
SponRad::~SponRad(){}

void SponRad::usage(){

  cout << "List of keywords for SponRad" << endl;
  cout << "&sponrad" << endl;
  cout << " int seed = 1234 " << endl;
  cout << " bool doLoss   = false" << endl;
  cout << " bool doSpread = false" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool SponRad::init(int rank, int size, map<string,string> *arg,  Beam *beam)
{


  map<string,string>::iterator end=arg->end();

  if (arg->find("seed")!=end) {seed = atoi(arg->at("seed").c_str());  arg->erase(arg->find("seed"));}
  if (arg->find("doLoss")!=end)    {doLoss    = atob(arg->at("doLoss").c_str());  arg->erase(arg->find("doLoss"));}
  if (arg->find("doSpread")!=end)  {doSpread  = atob(arg->at("doSpread").c_str());  arg->erase(arg->find("doSpread"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &sponrad" << endl; this->usage();}
    return false;
  }


  //  beam->initEField(rmax,ngrid,nz,nphi,lambda);


  return true;
}
 


::::::::::::::
src/Main/Time.cpp
::::::::::::::
#include "Time.h"

Time::Time()
{
  dotime=false;
  doscan=false;
  s0=0;
  slen=0;
  ds=1;
  sample=1;
  nslice=1;
  ns_node=1;
  noff_node=0;
  initialized=false;
}

Time::~Time(){}

void Time::usage(){

  cout << "List of keywords for TIME" << endl;
  cout << "&time" << endl;
  cout << " double s0 = 0" << endl;
  cout << " double slen   = 0" << endl;
  cout << " int sample = 1" << endl;
  cout << " bool time = true" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool Time::init(int inrank, int insize, map<string,string> *arg, Setup *setup)
{

  rank=inrank;
  size=insize;
  dotime=true;

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("s0")!=end)  {s0   = atof(arg->at("s0").c_str());  arg->erase(arg->find("s0"));}
  if (arg->find("slen")!=end){slen = atof(arg->at("slen").c_str()); arg->erase(arg->find("slen"));}
  if (arg->find("sample")!=end){sample = atoi(arg->at("sample").c_str()); arg->erase(arg->find("sample"));}
  if (arg->find("time")!=end){dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}
  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &time" << endl; this->usage();}
    return false;
  }

  initialized=true;  
  this->finishInit(setup);
  return true;


}


void Time::finishInit(Setup *setup)
{
  if (!initialized){return;}

  doscan=!dotime;

  ds=setup->getReferenceLength()*sample;
  nslice=static_cast<int> (round(slen/ds));
  if (nslice < size) { nslice = size ; }

  ns_node=nslice/size;
  if ((nslice % size)!=0) {ns_node++;}
  nslice=ns_node*size;
  noff_node=ns_node*rank;
  // needs to adjust the time window for harmonic conversions
  slen=ds*nslice;
  
  if (rank==0){
    cout << "Setting up time window of " << slen*1e6 << " microns with " << nslice <<" sample points..." << endl;
  }
}




int Time::getPosition(vector<double> *s)
{
  if (nslice<1){nslice=1;} 
  s->resize(nslice);
  for (int i=0;i <nslice; i++){
    s->at(i)=s0+i*ds;
  }
  return nslice;
}
::::::::::::::
src/Main/Track.cpp
::::::::::::::
#include "Track.h"
#include "Gencore.h"

Track::Track()
{
  zstop=1e9;
  output_step=1;
  sort_step=0;
  dumpFieldStep=0;
  dumpBeamStep=0;
  bunchharm=1;
}

Track::~Track(){}

void Track::usage(){

  cout << "List of keywords for TRACK" << endl;
  cout << "&track" << endl;
  cout << " double zstop = 1e9" << endl;
  cout << " double s0    = <taken from TIME module>" << endl;
  cout << " double slen  = <taken from TIME module>" << endl;
  cout << " int output_step  = 1" << endl;
  cout << " int field_dump_step  = 0" << endl;
  cout << " int beam_dump_step  = 0" << endl;
  cout << " int sort_step = 0" << endl;
  cout << " int bunchharm = 1" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool Track::init(int inrank, int insize, map<string,string> *arg, Beam *beam, vector<Field *> *field,Setup *setup, Lattice *lat, AlterLattice *alt,Time *time)
{
 
  rank=inrank;
  size=insize;
  bunchharm=1; //reset to default for each tracking
  
  bool isTime=time->isTime();
  bool isScan=time->isScan();
  double sample=time->getSampleRate();
  s0=time->getTimeWindowStart();
  slen=time->getTimeWindowLength();

  
  map<string,string>::iterator end=arg->end();

  if (arg->find("zstop")!=end)  {zstop= atof(arg->at("zstop").c_str());  arg->erase(arg->find("zstop"));}
  if (arg->find("s0")!=end)     {s0= atof(arg->at("s0").c_str());  arg->erase(arg->find("s0"));}
  if (arg->find("slen")!=end)   {slen= atof(arg->at("slen").c_str());  arg->erase(arg->find("slen"));}
  if (arg->find("output_step")!=end)   {output_step= atoi(arg->at("output_step").c_str());  arg->erase(arg->find("output_step"));}
  if (arg->find("field_dump_step")!=end)  {dumpFieldStep= atoi(arg->at("field_dump_step").c_str()); arg->erase(arg->find("field_dump_step"));}
  if (arg->find("beam_dump_step")!=end)   {dumpBeamStep = atoi(arg->at("beam_dump_step").c_str());  arg->erase(arg->find("beam_dump_step"));}
  if (arg->find("sort_step")!=end)   {sort_step= atoi(arg->at("sort_step").c_str());  arg->erase(arg->find("sort_step"));}
  if (arg->find("bunchharm")!=end)   {bunchharm= atoi(arg->at("bunchharm").c_str());  arg->erase(arg->find("bunchharm"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &track" << endl; this->usage();}
    return false;
  }


  if (output_step < 1) { output_step=1; }

  string file;
  setup->getRootName(&file);
  file.append(".out.h5");


 
  Undulator *und = new Undulator;

  lat->generateLattice(setup->getStepLength(),setup->getReferenceLength(),setup->getReferenceEnergy(),alt, und);  
  und->updateOutput(zstop,output_step);
  und->updateMarker(dumpFieldStep,dumpBeamStep,sort_step,zstop);

  beam->setBunchingHarmonicOutput(bunchharm);
  // call to gencore to do the actual tracking.  

  Gencore core;
  core.run(file.c_str(),beam,field,und,isTime,isScan);


  delete und;
   
  if  (rank==0) { cout << "End of Track" << endl;}
 
  return true;

}
::::::::::::::
src/Main/Wake.cpp
::::::::::::::
#include "Wake.h"
#include "Beam.h"

Wake::Wake()
{
  radius=2.5e-3;
  conductivity=0;
  relaxation=0;
  roundpipe=true;
}
Wake::~Wake(){}


void Wake::usage(){

  cout << "List of keywords for Wake" << endl;
  cout << "&wake" << endl;
  cout << " double radius = 2.5e-3" << endl;
  cout << " bool   roundpipe   = true" << endl;
  cout << " string material  = <empty>" << endl;
  cout << " double conductivity = 0e-6" << endl;
  cout << " double relaxation  = 0e-6" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool SponRad::init(int rank, int size, map<string,string> *arg,  Beam *beam)
{

  string material='';
  map<string,string>::iterator end=arg->end();

  if (arg->find("radius")!=end) {radius = atof(arg->at("radius").c_str());  arg->erase(arg->find("radius"));}
  if (arg->find("conductivity")!=end) {conductivity= atof(arg->at("conductivity").c_str());  arg->erase(arg->find("conductivity"));}
  if (arg->find("relaxation")!=end) {relaxation = atof(arg->at("relaxation").c_str());  arg->erase(arg->find("relaxation"));}
  if (arg->find("roundpipe")!=end)    {roundpipe    = atob(arg->at("roundpipe").c_str());  arg->erase(arg->find("roundpipe"));}
  if (arg->find("material")!=end){material = arg->at("material"); arg->erase(arg->find("material"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &wake" << endl; this->usage();}
    return false;
  }
  
  if ((material=="CU") || (material=="Cu") || (material=="cu")){
    conductivity=5.813e7;
    relaxation=8.1e-6;
  }

  if ((material=="AL") || (material=="Al") || (material=="al")){
    conductivity=3.571e7;
    relaxation=2.4e-6;
  }

  return true;
}
 
::::::::::::::
src/Util/BesselJ.cpp
::::::::::::::
/*
 *  BesselJ.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 12/5/11.
 *  Copyright 2011 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "BesselJ.h"

BesselJ::BesselJ(){}
BesselJ::~BesselJ(){}


double BesselJ::BesselJ0(double x)
{
    const double     p1=1.e0;
    const double     p2=-.1098628627e-2;
    const double     p3=.2734510407e-4;
    const double     p4= -.2073370639e-5;
    const double     p5=.2093887211e-6;
    const double     q1=-.1562499995e-1;
    const double     q2=.1430488765e-3;
    const double     q3=-.6911147651e-5;
    const double     q4=.7621095161e-6;
    const double     q5=-.934945152e-7;
    const double     r1=57568490574.e0;
    const double     r2=-13362590354.e0;
    const double     r3=651619640.7e0;
    const double     r4=-11214424.18e0;
    const double     r5=77392.33017e0;
    const double     r6=-184.9052456e0;
    const double     s1=57568490411.e0;
    const double     s2=1029532985.e0;
    const double     s3=9494680.718e0;
    const double     s4=59272.64853e0;
    const double     s5=267.8532712e0;
    const double     s6=1.e0;   
    
	if (fabs(x) < 8){
		double y=x*x;
		return (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));        
	}
	else
	{
		double ax=fabs(x);
		double z=8/ax;
		double y=z*z;
		double xx=ax-0.785398164;
		return sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-
									 z*sin(xx)* (q1+y*(q2+y*(q3+y*(q4+y*q5)))));       
	}
}


double BesselJ::BesselJ1(double x)
{
    const double    r1=72362614232.e0;
    const double    r2=-7895059235.e0;
    const double    r3=242396853.1e0;
    const double    r4=-2972611.439e0;
    const double    r5=15704.48260e0;
    const double    r6=-30.16036606e0;
    const double    s1=144725228442.e0;
    const double    s2=2300535178.e0;
    const double    s3=18583304.74e0;
    const double    s4=99447.43394e0;
    const double    s5=376.9991397e0;
    const double    s6=1.e0;
    const double    p1=1.e0;
    const double    p2=.183105e-2;
    const double    p3=-.3516396496e-4;
    const double    p4=.2457520174e-5;
    const double    p5=-.240337019e-6;
    const double    q1=.04687499995e0;
    const double    q2=-.2002690873e-3;
    const double    q3=.8449199096e-5;
    const double    q4=-.88228987e-6;
    const double    q5=.105787412e-6;
	
	if (fabs(x) < 8){
		double y=x*x;
		return x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));        
	}
	else
	{
		double ax=fabs(x);
		double z=8/ax;
		double y=z*z;
		double xx=ax-2.356194491;
		return sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-
									 z*sin(xx)* (q1+y*(q2+y*(q3+y*(q4+y*q5)))));       
	}
}





double BesselJ::value(int n, double x)
{

	
	double bessj=0;
	if (n==0){
  	  return this->BesselJ0(x);
	}
	if (n==1) {
	  return this->BesselJ1(x);
	}
	
	if (x==0) {
		return 0;
	}
	
    double ax = fabs(x);
	
	if (ax > static_cast<double> (n)){
		double tox=2/ax;
		double bjm=this->BesselJ0(ax);
		double bj =this->BesselJ1(ax);
		for (int j=1; j<n; j++) {
			double bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
	    bessj=bj;
	} else {
		double tox=2./ax;
		int m=2*(n+floor(sqrt(static_cast<double>(40*n)))/2);
	    bessj=0;
		int jsum=0;
		double sum=0;
		double bjp=0;
		double bj=1;
		for (int j=m; j>0; j--) {
			double bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj)>1e10){
				bj*=1e-10;
				bjp*=1e-10;
				bessj*=1e-10;
				sum*=1e-10;
			}
			if (jsum!=0){sum+=bj;}
			jsum=1-jsum;
			if (j==n) {bessj=bjp;}
		}
		sum=2*sum-bj;
		bessj=bessj/sum;
				 
	}
	
    if ((x<0)&& ((n%2)==1 )) {bessj=-bessj;}
		
	return bessj;
}::::::::::::::
src/Util/GaussHermite.cpp
::::::::::::::
/*
 *  GaussHermite.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 5/14/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "GaussHermite.h"

// Constructor/Destructor

GaussHermite::GaussHermite(){}

GaussHermite::~GaussHermite(){}



void GaussHermite::loadGauss(complex<double> *field, FieldSlice *slice, double dgrid, int ngrid)
{

	
  	
  double k=4.*asin(1)/slice->lambda*static_cast<double>(slice->harm);
  double z0=-slice->z0;
  double w0=slice->w0;
  double zr=w0*w0*k*0.5;  // Zr=w0^2*k/2
  double f0=sqrt(k*zr/(zr*zr+z0*z0));          // is same as  sqrt(2)/w(z)


  complex<double> qz=complex<double>(-z0,zr);  // q(z)= z-z0 + i zr , see Siegman p.664
  complex<double> q0=complex<double>(0,  zr);
 
  complex<double> coef=complex<double>(0,1)*0.5*k/qz;		

  // some normalization crap
  double unit=sqrt(slice->power*vacimp)*k/eev;  // P=|E|^2/Z0 -> u = k (e/mc^2) E	
  double norm=f0/sqrt(2.*asin(1.)*this->fac(slice->nx)*this->fac(slice->ny)*(1<<(slice->nx+slice->ny))); // f0/sqrt(pi 2^(nx+ny) nx! ny!)	
  complex<double> zscale=unit*norm*complex<double>(cos(slice->phase),sin(slice->phase));
  


  double xmid=dgrid+slice->xcen;
  double ymid=dgrid+slice->ycen;
  double dxy=2.*dgrid/(ngrid-1.);
	
  double kx=k*slice->xangle;
  double ky=k*slice->yangle;
	
  for (int iy=0;iy<ngrid;iy++){
     double y=iy*dxy-ymid;
     double Hy=this->Hn(f0*y,slice->ny);
     double y2=y*y;
     double phiy=kx*y;
     for(int ix=0;ix<ngrid;ix++){ // x is inner loop
	  double x=ix*dxy-xmid;
	  double Hx=this->Hn(f0*x,slice->nx);
	  double r2=y2+x*x;
       	  complex<double> phi=complex<double>(0,phiy+kx*x);
	  field[iy*ngrid+ix]=zscale*exp(-coef*r2+phi)*Hx*Hy;
     } 
  }
  return;	
		
}

int GaussHermite::fac(int n)
{
	if (n<=1){
		return 1;
	}
    return n*this->fac(n-1);
}


double GaussHermite::Hn(double x, int n)
{
    if (n<=0) {
	  return 1;
    } else {
      if ( n == 1) {
	return 2*x;
      } else {
	return 2.*x*this->Hn(x,n-1)-2*(n-1)*this->Hn(x, n-2);
      }
    }
}
::::::::::::::
src/Util/Hammerslay.cpp
::::::::::::::
/*
 *  Hammerslay.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 5/26/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "Hammerslay.h"

// constructor/destructor


Hammerslay::Hammerslay(unsigned int base_in)
{
	double bases[26] = {2,3,5,7,11,13,17,19,23,29,
		31,37,41,43,47,53,59,61,67,
		71,73,79,83,89,97,101};
	
	base=bases[base_in];
	idx=0;
}

Hammerslay::~Hammerslay(){}


void Hammerslay::set(unsigned int i)
{
	idx=i;
   if (idx<0) {
	   idx=0;
   }
}

double Hammerslay::getElement()
{
	double xs=0;
	double xsi=1.0;
	unsigned int    i1;
	
	unsigned int i2=++idx;
	
	do{
		xsi=xsi/base;
		i1=i2/static_cast< unsigned int > (base); 
		xs+=(i2-base*i1)*xsi;
		i2=i1;
	} while(i2>0);
	
	return xs;
}					               

::::::::::::::
src/Util/Inverfc.cpp
::::::::::::::
/*
 *  Inverfc.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 09.12.11.
 *  Copyright 2011 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "Inverfc.h"
#include "math.h"

Inverfc::Inverfc(){}
Inverfc::~Inverfc(){}

double Inverfc::value(double y)
{

	
	/*
	 inverted error function
	 original author:  Takuya OOURA
	 Takuya OOURA, Research Institute for Mathematical Sciences // 
	 Kyoto University, Kyoto 606-01 Japan // 
	 Email : ooura@kurims.kyoto-u.ac.jp (orooura@mmm.t.u-tokyo.ac.jp ). 
	 reference: http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
	 
	 function is used for generating gaussian distribution avoiding the
	 joint-propability approach with two uniform distributed random number distributions
	 */
	
	double z = y;
	if (y > 1) { z = 2 - y; }
	
	double w = 0.916461398268964 - log(z);
	double u = sqrt(w);
	double s = (log(u) + 0.488826640273108) / w ;
	double t = 1 / (u + 0.231729200323405);
	double x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
	+   ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
		   +   0.150689047360223) * t + 0.116065025341614) * t + 
		 +   0.499999303439796) * t;
	
	t = 3.97886080735226 / (x + 3.97886080735226);
	u = t - 0.5;
	s = (((((((((0.00112648096188977922 * u + 
				 +   1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
			   +   7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
			 +   0.00339721910367775861) * u - 0.011274916933250487) * u - 
		   +   0.0118598117047771104) * u + 0.0142961988697898018) * u + 
		 +   0.0346494207789099922) * u + 0.00220995927012179067;
	s = ((((((((((((s * u - 0.0743424357241784861) * u - 
				   +   0.105872177941595488) * u + 0.0147297938331485121) * u + 
				 +   0.316847638520135944) * u + 0.713657635868730364) * u + 
			   +   1.05375024970847138) * u + 1.21448730779995237) * u + 
			 +   1.16374581931560831) * u + 0.956464974744799006) * u + 
		   +   0.686265948274097816) * u + 0.434397492331430115) * u + 
		 +   0.244044510593190935) * t - 
	+   z * exp(x * x - 0.120782237635245222);
	x = x+s * (x * s + 1);
	if (y > 1) {x = -x;}
	return x;
	

}::::::::::::::
src/Util/RandomU.cpp
::::::::::::::
/*
 *  RandomU.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 5/26/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */


// Routine Ran2 from Numerical Receipe

#include "RandomU.h"

// constructor + destructor

RandomU::RandomU(unsigned int istart)
{

const int ia1=40014;
const int im1=2147483563;
const int iq1=53668;
const int ir1=12211;
const int ntab=32;
int k;
iseed = (1 < abs(istart)) ? abs(istart) : 1 ;
iseed2=iseed;
for (int i=ntab+7;i>=0;i--){
	k=iseed/iq1;
	iseed=ia1*(iseed-k*iq1)-k*ir1;
	if (iseed < 0 ) iseed+=im1;
		if (i < ntab) iv[i]=iseed;
			}
iy=iv[0]; 
}

RandomU::~RandomU(){}



void RandomU::set(unsigned int istart)
{
  return;
}

double RandomU::getElement()
{
	
	const int ia1=40014;
	const int ia2=40692;
	const int im1=2147483563;
	const int im2=2147483399;
	const int imm1=2147483562;
	const int iq1=53668;
	const int iq2=52774;
	const int ir1=12211;
	const int ir2=3791;
	const int ndiv=67108862;
	const double am=1./2147483563;
	const double rnmx=1.-1.2e-40;	
	
	int k=iseed/iq1;
	iseed=ia1*(iseed-k*iq1)-k*ir1;
	if (iseed <0) iseed=iseed+im1;
	k=iseed2/iq2;
	iseed2=ia2*(iseed2-k*iq2)-k*ir2;
	if (iseed2 < 0) iseed2=iseed2+im2;
	int j=iy/ndiv;
	iy=iv[j]-iseed2;
	iv[j]=iseed; 
	if ( iy < 1 ) iy=iy+imm1;
	return (am*iy < rnmx) ? am*iy : rnmx;
}
::::::::::::::
src/Util/Sorting.cpp
::::::::::::::
#include "Sorting.h"

extern bool MPISingle;
 
#include <algorithm>

Sorting::Sorting()
{
  dosort=false;
  doshift=false;
}

Sorting::~Sorting(){
}



void Sorting::init(int rank_in,int size_in, bool doshift_in, bool dosort_in)
{
  rank=rank_in;  // rank of the node (0-size-1)
  size=size_in;  // size of the mpi run
  dosort=dosort_in;        // flag to do the sorting at all
  doshift=doshift_in;      // flag to apply the global shifting, meaning that the radiationfield is rather moved than the beam
  return;

}


void Sorting::configure(double s0_in, double slicelen_in, double sendmin_in, double sendmax_in, double keepmin_in, double keepmax_in, bool frame_in)
{
  s0=s0_in;      // offset in particle distribution of the given node
  slen=slicelen_in; // distance between two adjacent slice in particle distribution vector
  sendmin=sendmin_in;  // s < sendmin -> particle is send to previous node
  sendmax=sendmax_in;  // same in the other direction
  keepmin=keepmin_in;  // s < keepmin -> particle is deleted from internal slice
  keepmax=keepmax_in;  // same in the other direction
  globalframe=frame_in; // flag whether particle position is global or not
  return;

}





int Sorting::sort(vector <vector <Particle> > * recdat){

  //   cout << "Rank: " << rank << " s0: " << s0 << " slen: " << slen << " sendmin: " << sendmin << " sendmax: " << sendmax << " keepmin: " << keepmin << " keepmax: " << keepmax << " globaL : " << globalframe << endl;

  if (!dosort) {
    if (rank==0) {cout << "*** Warning: Sorting only enabled for one-2-one simulations" << endl;}
    return 0;
  }
  if (rank==0) {cout << "Sorting..." << endl; }
  int shift =0;
  // step one - calculate global shift and see whether it can be adjusted by shifting radiation instead
  //////  int shift = this->centerShift(recdat);


  // step two - push outside particles to other nodes
  this->globalSort(recdat);
  // step three - sort within the given node
  this->localSort(recdat);


  return shift;

}





void Sorting::localSort(vector <vector <Particle> > * recdat)  // most arguments are now part of the class
{
 
  Particle p;  

  double invslen=1./slen;
  cout << "Testfloor: " << floor(-0.2) << " and " << static_cast<int>(floor(-0.2)) << endl;
  // note that global sorting comes first. Therefore all particles are staying in the same domain 
  
  vector<int> count,count2;
  count.resize(recdat->size());
  count2.resize(recdat->size());
  for  (int i=0; i< count.size(); i++){ 
     count[i]=0;
     count2[i]=0;
  }


  for (int a=0;a<recdat->size();a++){  //Run over the slices 
    int b=0;
    while ( b < recdat->at(a).size()){
      double theta=recdat->at(a).at(b).theta;
      int atar=static_cast<int>(floor(theta*invslen));   // relative target slice. atar = 0 -> stays in same slice

      if (atar!=0) {     // particle not in the same slice

	  p.theta=recdat->at(a).at(b).theta-slen*(atar);
	  p.gamma=recdat->at(a).at(b).gamma;
	  p.x    =recdat->at(a).at(b).x;
	  p.y    =recdat->at(a).at(b).y;
	  p.px   =recdat->at(a).at(b).px;
	  p.py   =recdat->at(a).at(b).py;
	  recdat->at(a+atar).push_back(p);  // pushing particle in correct slice

          count[a+atar]+=1;
          count2[a]-=1;

	  // eliminating the current particle
 	  int ilast=recdat->at(a).size()-1;
      	  recdat->at(a).at(b).theta=recdat->at(a).at(ilast).theta;
	  recdat->at(a).at(b).gamma=recdat->at(a).at(ilast).gamma;
	  recdat->at(a).at(b).x    =recdat->at(a).at(ilast).x;
	  recdat->at(a).at(b).y    =recdat->at(a).at(ilast).y;
	  recdat->at(a).at(b).px   =recdat->at(a).at(ilast).px;
	  recdat->at(a).at(b).py   =recdat->at(a).at(ilast).py;
	  recdat->at(a).pop_back();
      } else {
	b++;
      }
    }
  }	  
	  //	  count++;
      
    



  for (int i=0; i<count.size(); i++){
    cout<<"Rank: " << rank << " Slice: " << i << " received: " << count[i] << " send: "<<count2[i]<<endl;
  } 
   
  return;
}



// routine which moves all particles, which are misplaced in the given domain of the node to other nodes.
// the methods is an iterative bubble sort, pushing excess particles to next node. There the fitting particles are collected the rest moved further.

void Sorting::globalSort(vector <vector <Particle> > *rec)
{

  this->fillPushVectors(rec);   // here is the actual sorting to fill the vectore pushforward and pushbackward
  if (rank==(size-1)) { pushforward.clear(); }        
  if (rank==0) { pushbackward.clear(); }
  if (size==1) { return; } // no need to transfer if only one node is used.
  
  cout << "Rank: " << rank << " - Forward: " << pushforward.size()/6 << " - Backward: " << pushbackward.size()/6 << endl;

  int maxiter=size-1;  
  int nforward=pushforward.size();
  int nbackward=pushbackward.size();
  int ntotal=nforward+nbackward;
  int nreduce=0;
  MPI::COMM_WORLD.Allreduce(&ntotal,&nreduce,1,MPI::INT,MPI::SUM);  // get the info on total number of particles shifted among nodes

  if (nreduce == 0){ return; }	

  while(maxiter>0){
    if (rank==0) {cout << "Sorting: Transferring " << nreduce/6 << " particles to other nodes at iteration " << size-maxiter << endl;}



  // step one - pairing ranks (0,1) (2,3) etc, the last rank, if uneven last element should clear its pushforward.
     bool transfer = true;
     if (((rank % 2) == 0) && (rank == (size -1))) { transfer = false; }  // last rank for an uneven count of cores -> no forward transmission

     if ((rank % 2)==0) {   // even ranks sending
       if (transfer) {this->send(rank+1,&pushforward);}        // sends its forward particles to higher node
       pushforward.clear();                                    // no need for the record, data has been sent
       if (transfer) {this->recv(rank+1,rec,&pushbackward);}   // catch particles which are sent back. In self->recv either get the particles or put them in backword array
     }  else {   // odd ranks receiving - there will be always a smaller even rank therefore always receiving
	   this->recv(rank-1,rec,&pushforward);
           this->send(rank-1,&pushbackward);
	   pushbackward.clear();
     }
      
  // step two - pairing ranks (1,2) (3,4) etc

     transfer == true;
     if (((rank % 2) == 1) && (rank == (size -1))) { transfer = false; }  // last core for an even number
     if (rank==0) { transfer = false; }                                  // as well as first core.
     
     if ((rank % 2)==1) {  // odd ranks sending
       if (transfer) {this->send(rank+1,&pushforward);}
	 pushforward.clear();
	 if (transfer) {this->recv(rank+1,rec,&pushbackward);}
     } else {  // even one receiving - here check for the very first node
       if (transfer){
	   this->recv(rank-1,rec,&pushforward); 
           this->send(rank-1,&pushbackward);
           pushbackward.clear();
       }
     }
 
     maxiter--;
     nforward=pushforward.size();
     nbackward=pushbackward.size();
     ntotal=nforward+nbackward;
     nreduce=0;
     MPI::COMM_WORLD.Allreduce(&ntotal,&nreduce,1,MPI::INT,MPI::SUM);
     if (nreduce == 0){  return; }
  }
  pushforward.clear();
  pushbackward.clear();
  return;



  

}


void Sorting::send(int target, vector<double> *data)
{
  int ndata=data->size();
  MPI::COMM_WORLD.Send( &ndata,1,MPI::INT,target,0);
  if (ndata == 0) {
    return;
  }
  MPI::COMM_WORLD.Send(&data->front(),ndata,MPI::DOUBLE,target,0);
}


 void Sorting::recv(int source, vector <vector <Particle> > *rec ,vector<double> *olddata)
{

  double shift=slen;
  if(globalframe){shift=0;}

   MPI::Status status;
   int ndata=0;


   MPI::COMM_WORLD.Recv(&ndata,1,MPI::INT,source,0,status);
   if (ndata==0) {  // no data received.
     return;
   }

   vector<double> data (ndata);
   MPI::COMM_WORLD.Recv(&data.front(),ndata,MPI::DOUBLE,source,0,status); // geting the particles from adjacent node.
 

   //Determines whether the data needs to be pushed forward or backwards or stored in the correct slices
   int np=rec->size();                      // number of slices
   if (source>rank) {                       // data are coming from higher node -> particles are pushed backward
     for (int a=0;a<ndata;a+=6) {
       double s = s0+slen*(np-1)+data[a];  // get the actual positionassume that backward the particles are placed in the last slice !
       if (s<sendmin){
         olddata->push_back(data[a]+shift*np);  // if the particle isstillpushed through than it phase is adjusted by slen*slicenumber
	 for (int b=1;b<6;b++){
           olddata->push_back(data[a+b]);
         }
       } 
       if (s>=keepmin){
         Particle par;
         par.theta=data[a];
         par.gamma=data[a+1];
         par.x    =data[a+2];
         par.y    =data[a+3];
         par.px   =data[a+4];
         par.py   =data[a+5];
         rec->at(np-1).push_back(par); // add at last slice;
       }
     }
   } else {    // particles are coming from a lower node -> particles are pushed forward
     for (int a=0;a<ndata;a+=6) {
       double s = s0+data[a];  // get the actual position assume that forward the particles are placed in the first slice !
       if (s>sendmax){
	 olddata->push_back(data[a]-shift*np);
	 for (int b=1;b<6;b++){
           olddata->push_back(data[a+b]);
         }
       } 
       if (s<=keepmax){
         Particle par;
         par.theta=data[a];
         par.gamma=data[a+1];
         par.x    =data[a+2];
         par.y    =data[a+3];
         par.px   =data[a+4];
         par.py   =data[a+5];
         rec->at(0).push_back(par); // add to first slice
       }
     }
   }

   return;
}




void Sorting::fillPushVectors(vector< vector <Particle> >*rec)
{
  cout << "Rank: " << rank << " sendmin: " << sendmin << " sendmax: " << sendmax << endl;
  cout << "Rank: " << rank << " keepmin: " << keepmin << " keepmax: " << keepmax << endl;
  //step one - fill the push vectors
  pushforward.clear();
  pushbackward.clear();
  
  int nsize=rec->size();
  double shift=slen;  // flag to indicate correction in position because each slice has its own position ( 3pi in slice 5 is pi in slice 6}
  if (globalframe) {shift = 0;} // don't change position if it is a global frame (e.g. when importing elegant distibution)
  
  int count = 0;

  for (int i = 0; i < nsize; i++){  // loop over slices
    for (int j = 0; j < rec->at(i).size(); j++){ // loop over particles in slice
      double s = s0+slen*i+rec->at(i).at(j).theta;  // get the actual position
      if (s<sendmin){
	pushbackward.push_back(rec->at(i).at(j).theta+(i+1)*shift); 
	// example: in first slice phase is -1. Particle needs to be sent to previous node to the last slice there. Therefore
        // the offset is just the slice  length (normally 2 pi). In the second slice it needs to be a value of -7 < - 2 pi.
	pushbackward.push_back(rec->at(i).at(j).gamma);
	pushbackward.push_back(rec->at(i).at(j).x);
	pushbackward.push_back(rec->at(i).at(j).y);
	pushbackward.push_back(rec->at(i).at(j).px);
	pushbackward.push_back(rec->at(i).at(j).py);
      }
      if (s>sendmax){
	pushforward.push_back(rec->at(i).at(j).theta-(nsize-i)*shift); 
	// example: last slice has value of 7 > 2 pi to be pushed to next node in the first slice there. The offset correction would be - 2 pi, which is given by
	// -(nlen-i)*slicelength  with i = nlen-1 for the last slice
	pushforward.push_back(rec->at(i).at(j).gamma);
	pushforward.push_back(rec->at(i).at(j).x);
	pushforward.push_back(rec->at(i).at(j).y);
	pushforward.push_back(rec->at(i).at(j).px);
	pushforward.push_back(rec->at(i).at(j).py);

      }
    }

    int j=0;

    while (j<rec->at(i).size()){
        double s = s0+slen*i+rec->at(i).at(j).theta;  // get the actual position
        if ((s<keepmin)||(s>keepmax)){
          count++;
 	  int ilast=rec->at(i).size()-1;
	  rec->at(i).at(j).theta=rec->at(i).at(ilast).theta;
	  rec->at(i).at(j).gamma=rec->at(i).at(ilast).gamma;
	  rec->at(i).at(j).x    =rec->at(i).at(ilast).x;
	  rec->at(i).at(j).y    =rec->at(i).at(ilast).y;
	  rec->at(i).at(j).px   =rec->at(i).at(ilast).px;
	  rec->at(i).at(j).py   =rec->at(i).at(ilast).py;
	  rec->at(i).pop_back();
	} else {
	  j++;
	}
    }

  }
  cout << "Rank: " <<rank << " Deleted: " << count << " Forward: " << pushforward.size()/6 << " Backward: " <<pushbackward.size()/6 << endl;
}


int Sorting::centerShift(vector <vector <Particle> > * recdat)
{
  if (!doshift){ return 0;}

  double invslen=1./slen;
  double shift = 0;
  double part  = 0;

  // calculate on the node the amount of total transferred particles

  for (int a=0;a<recdat->size();a++){  //Run over the slices 
    for (int b=0;b<recdat->at(a).size();b++) {  //Loop over the partiles in the slice
      part+=1.0;
      shift+=floor(recdat->at(a).at(b).theta*invslen);   // relative target slice. = 0 -> stays in same slice
    }
  }

  double tshift,tpart;
  if (MPISingle){
    tshift=shift;
    tpart=part;

  } else {
     MPI::COMM_WORLD.Allreduce(&shift,&tshift,1,MPI::DOUBLE,MPI::SUM);  // get the info on total number of particles shifted among nodes
     MPI::COMM_WORLD.Allreduce(&part,&tpart,1,MPI::DOUBLE,MPI::SUM);  // get the info on total number of particles shifted among nodes
  }

  int nshift=-static_cast<int>(round(tshift/tpart));  // instead of moving more than half the particles forward, it is easier to move the field backwards. Also theta needs to be corrected

  for (int a=0;a<recdat->size();a++){  //Run over the slices 
    for (int b=0;b<recdat->at(a).size();b++) {  //Loop over the partiles in the slice
      recdat->at(a).at(b).theta+=static_cast<double>(nshift)*slen;   // adjust theta position because the field is moved instead of particles
    }
  }

  return nshift;
}
::::::::::::::
src/Util/StringProcessing.cpp
::::::::::::::
/*
 *  StringProcessing.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 9/10/10.
 *  Copyright 2010 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "StringProcessing.h"

StringProcessing::StringProcessing(){}
StringProcessing::~StringProcessing(){}

//---------------------------------------------------------------------------
// some general string utility functions

void StringProcessing::chop(string str, vector<string> *list){
	
	size_t pos;
	
	list->clear();
	
	
	while((pos=str.find_first_of(","))!=string::npos){
	    list->push_back(str.substr(0,pos));
	    str.erase(0,pos+1);
	}
        list->push_back(str);
	
	for (int i=0; i<list->size();i++){
	  this->trim(list->at(i));
        }

	return;
}

void StringProcessing::trim(string &str){
    size_t startpos = str.find_first_not_of(" \t");
	size_t endpos = str.find_last_not_of(" \t");
	if(( string::npos == startpos ) || ( string::npos == endpos)){
		str="";
	} else {
		str = str.substr( startpos, endpos-startpos+1 );
	}
	return;
}


bool StringProcessing::atob(string in){
    
	bool ret=false;
	if ((in.compare("1")==0)||(in.compare("true")==0)||(in.compare("t")==0)) { ret=true; }
	return ret;
}


void StringProcessing::reference(string in , double *val, string *ref)
{
  size_t pos=in.find_first_of("@");
  if (pos!=string::npos){
    *ref=in.erase(0,pos+1);
    return;
  } else {
    *ref="";
    *val=atof(in.c_str());
    return;
  }
}
