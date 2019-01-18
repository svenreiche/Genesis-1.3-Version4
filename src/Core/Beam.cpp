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
  eloss.resize(nsize);
  for (int i=0; i<nsize; i++){
    eloss[i]=0;
  }

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
  efld.resize(nz*ns);

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
  col.apply(this,und,delz);         // apply effect of collective effects
  
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
    efld[ioff+is]=eloss[is];  

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

