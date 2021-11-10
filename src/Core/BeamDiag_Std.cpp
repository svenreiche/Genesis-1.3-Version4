#include <mpi.h>
#include <hdf5.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include "BeamDiag_Std.h"
#include "BeamDiag_HDF5helper.h"

BeamDiag_Std::BeamDiag_Std() {
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);

	// force call to class-specific function 'init' before 'do_diag'
	is_initialized_=false;
	// force call to class-specific function 'configure' before 'do_diag'
	is_configured_=true; // <-- no required configuration variables at the moment

	bharm=1;

	/* defaults as previously in Beam::Beam */
	do_global_stat=false;
	doCurrent=false;
	doSpatial=true;
	doEnergy=true;
	doAux=true;

	verbose_ = false;

	ns_=0;
	idx_=0;
}

BeamDiag_Std::~BeamDiag_Std() {

}

std::string BeamDiag_Std::to_str(void) const {
	stringstream ss;
	ss << "BeamDiag_Std";
	return(ss.str());
}

void BeamDiag_Std::init(int nz, int ns) {
	ns_ = ns;
	idx_ = 0;
	is_initialized_=true;

  /* code from Beam::initDiagnostics */
  zpos.resize(nz);
  if (doSpatial){
    xavg.resize(nz*ns);
    xsig.resize(nz*ns);
    yavg.resize(nz*ns);
    ysig.resize(nz*ns);
    pxavg.resize(nz*ns);
    pyavg.resize(nz*ns);
  } else {
    xavg.resize(0);
    xsig.resize(0);
    yavg.resize(0);
    ysig.resize(0);
    pxavg.resize(0);
    pyavg.resize(0);
  }
  if (doEnergy) {
    gavg.resize(nz*ns);
    gsig.resize(nz*ns);
  } else {
    gavg.resize(0);
    gsig.resize(0);
  }
  bunch.resize(nz*ns); 
  bphi.resize(nz*ns);
  if (doAux){
    efld.resize(nz*ns);
  } else {
    efld.resize(0);
  }

  bx.resize(ns);
  by.resize(ns);
  ax.resize(ns);
  ay.resize(ns);
  ex.resize(ns);
  ey.resize(ns);
  if (doCurrent) {
    cu.resize(ns*nz);
  } else {
    cu.resize(ns);
  }
  
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
  if (doEnergy){ 
    tgavg.resize(nz);
    tgsig.resize(nz);
  } else {
    tgavg.resize(0);
    tgsig.resize(0);
  }
  if (doSpatial){
    txavg.resize(nz);
    txsig.resize(nz);
    tyavg.resize(nz);
    tysig.resize(nz);
  } else {
    txavg.resize(0);
    txsig.resize(0);
    tyavg.resize(0);
    tysig.resize(0);
  }
  tbun.resize(nz);
}

/* configuration functions */
void BeamDiag_Std::setBunchingHarmonicOutput(int harm_in){bharm=harm_in;}
int BeamDiag_Std::getBunchingHarmonics(){return bharm;}
void BeamDiag_Std::set_global_stat(bool in){do_global_stat=in;}
bool BeamDiag_Std::get_global_stat(void){return(do_global_stat);}
void BeamDiag_Std::setOutput(bool noCurrent_in, bool noEnergy_in, bool noSpatial_in, bool noAux_in) {
  doCurrent = !noCurrent_in;
  doSpatial = !noSpatial_in;
  doEnergy = !noEnergy_in;
  doAux = !noAux_in;
}

/* Called when the electron beam is diagnosed, typically after every integration step */
void BeamDiag_Std::do_diag(const Beam *beam, double z) {
	int ioff=idx_*ns_;

	if((is_configured_==false) || (is_initialized_==false)){
		cout << "Error (BeamDiag_Std): do_diag called, but class instance is not configured" << endl;
	}

	if(verbose_ && (0==my_rank_)) {
		cout << "--> BeamDiag_Std::do_diag called" << endl;
	}

	/* FIXME */
	if(beam->beam.size() != ns_) {
		cout << "*** error: assertion ***" << endl;
	}

  /* Code from Beam::diagnostics */
  double acc_cur,acc_g,acc_g2,acc_x,acc_x2,acc_y,acc_y2;
  complex<double> acc_b=(0,0);
  acc_cur=0;
  acc_g=0;
  acc_g2=0;
  acc_x=0;
  acc_x2=0;
  acc_y=0;
  acc_y2=0;

  zpos[idx_]=z;

  for (int is=0; is < ns_; is++){
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

    unsigned int nsize=beam->beam.at(is).size();
    for (int i=0;i < nsize;i++){
      double xtmp=beam->beam.at(is).at(i).x;
      double ytmp=beam->beam.at(is).at(i).y;
      double pxtmp=beam->beam.at(is).at(i).px;
      double pytmp=beam->beam.at(is).at(i).py;
      double gtmp=beam->beam.at(is).at(i).gamma;
      double btmp=beam->beam.at(is).at(i).theta;
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
 
    if (do_global_stat){
      acc_cur+=beam->current[is];
      acc_g+= bgavg*beam->current[is];
      acc_g2+=bgsig*beam->current[is];
      acc_x+= bxavg*beam->current[is];
      acc_x2+=bxsig*beam->current[is];
      acc_y+= byavg*beam->current[is];
      acc_y2+=bysig*beam->current[is];
    }

    bbavg=sqrt(bi*bi+br*br)*scl;
    bbphi=atan2(bi,br);
    bgsig=sqrt(fabs(bgsig-bgavg*bgavg));
    bxsig=sqrt(fabs(bxsig-bxavg*bxavg));
    bysig=sqrt(fabs(bysig-byavg*byavg));

    if (doCurrent){
      cu[ioff+is]=beam->current[is];
    }
    if (doEnergy){
      gavg[ioff+is]=bgavg;
      gsig[ioff+is]=bgsig;
    }
    if (doSpatial){
      xavg[ioff+is]=bxavg;
      xsig[ioff+is]=bxsig;
      yavg[ioff+is]=byavg;
      ysig[ioff+is]=bysig;
      pxavg[ioff+is]=bpxavg;
      pyavg[ioff+is]=bpyavg;
    }
    bunch[ioff+is]=bbavg;
    bphi[ioff+is]=bbphi;
    if (doAux){
       efld[ioff+is]=beam->eloss[is];  
    }
    //    partcount[ioff+is]=nsize;

    for (int ih=1; ih<bharm;ih++){   // calculate the harmonics of the bunching
      br=0;
      bi=0;
      for (int i=0;i < nsize;i++){
        double btmp=static_cast<double>(ih+1)*beam->beam.at(is).at(i).theta;
        br+=cos(btmp);
        bi+=sin(btmp);
      }
      bh[ih-1][ioff+is]=sqrt(bi*bi+br*br)*scl;
      ph[ih-1][ioff+is]=atan2(bi,br);
    }
  }


  // accumulate all data from the cores
  double temp=0;
  int size;      
  if(do_global_stat) {
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      if (size>1){
	MPI_Allreduce(&acc_cur, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_cur=temp;
	MPI_Allreduce(&acc_g, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_g=temp;
	MPI_Allreduce(&acc_g2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_g2=temp;
	MPI_Allreduce(&acc_x, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_x=temp;
	MPI_Allreduce(&acc_x2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_x2=temp;
	MPI_Allreduce(&acc_y, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_y=temp;
	MPI_Allreduce(&acc_y2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	acc_y2=temp;
      }
      if (acc_cur==0){
         acc_cur=1.; 
      }
      acc_g/= acc_cur;
      acc_g2/=acc_cur;
      acc_x/= acc_cur;
      acc_x2/=acc_cur;
      acc_y/= acc_cur;
      acc_y2/=acc_cur;
      acc_g2=sqrt(abs(acc_g2-acc_g*acc_g));
      acc_x2=sqrt(abs(acc_x2-acc_x*acc_x));
      acc_y2=sqrt(abs(acc_y2-acc_y*acc_y));
      if (doEnergy){
        tgavg[idx_]=acc_g;
	tgsig[idx_]=acc_g2;
      }
      if (doSpatial){
        txavg[idx_]=acc_x;
        txsig[idx_]=acc_x2;
        tyavg[idx_]=acc_y;
        tysig[idx_]=acc_y2;
     }
  }

	idx_++;
}
void BeamDiag_Std::do_initial_diag(const Beam *beam)
{
  // double gx,gy,gammax,gammay;
  double x1,y1,x2,y2,px1,py1,px2,py2,g1,xpx,ypy;

  int ds=beam->beam.size();

	if((is_configured_==false) || (is_initialized_==false)){
		cout << "Error (BeamDiag_Std): do_initial_diag called, but class instance is not configured" << endl;
	}

	if(verbose_ && (0==my_rank_)) {
		cout << "--> BeamDiag_Std::do_initial_diag called" << endl;
	}

  for (int is=0; is<ds;is++){
    if (!doCurrent){
      cu[is]=beam->current[is];
    }
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
    unsigned int nsize=beam->beam.at(is).size();
    for (int i=0;i < nsize;i++){
      double xtmp=beam->beam.at(is).at(i).x;
      double ytmp=beam->beam.at(is).at(i).y;
      double pxtmp=beam->beam.at(is).at(i).px;
      double pytmp=beam->beam.at(is).at(i).py;
      double gtmp=beam->beam.at(is).at(i).gamma;
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
    ax[is]=-(xpx-x1*px1)/ex[is];
    ay[is]=-(ypy-y1*py1)/ey[is];
     
    // gx=(1+ax[is]*ax[is])/bx[is];
    // gy=(1+ay[is]*ay[is])/by[is];

  }

  return;
}


/* This function is called when the .out.h5 file is written */
void BeamDiag_Std::output(hid_t parentobj) {
	if(verbose_ && (0==my_rank_)) {
		cout << "--> BeamDiag_Std::output called" << endl;
	}

  hid_t gid=parentobj;
  hid_t gidsub;
  BeamDiag_HDF5Helper hh;

  hh.set_ds(ns_);
  hh.set_s0(my_rank_*ns_);

  /* Code from Output::writeBeamBuffer */
  // step 2 - write individual datasets
  if (doEnergy){
    hh.writeBuffer(gid, "energy"," ",&gavg);
    hh.writeBuffer(gid, "energyspread"," ", &gsig);
  }
  if (doSpatial){
    hh.writeBuffer(gid, "xposition","m",&xavg);
    hh.writeBuffer(gid, "yposition","m",&yavg);
    hh.writeBuffer(gid, "pxposition","rad", &pxavg);
    hh.writeBuffer(gid, "pyposition","rad", &pyavg);
    hh.writeBuffer(gid, "xsize","m", &xsig);
    hh.writeBuffer(gid, "ysize","m", &ysig);
  }
  hh.writeBuffer(gid, "bunching"," ",&bunch);
  hh.writeBuffer(gid, "bunchingphase","rad", &bphi);
  if (doAux){
    hh.writeBuffer(gid, "efield","eV/m", &efld);
  }
  
  hh.writeBuffer(gid, "betax","m",&bx);
  hh.writeBuffer(gid, "betay","m",&by);
  hh.writeBuffer(gid, "alphax","rad",&ax);
  hh.writeBuffer(gid, "alphay","rad",&ay);
  hh.writeBuffer(gid, "emitx","m",&ex);
  hh.writeBuffer(gid, "emity","m",&ey);
  hh.writeBuffer(gid, "current","A",&cu);

  // int bh=beam->getBunchingHarmonics();
  char bgroup[20];
  for (int i=1; i<bharm;i++){
    sprintf(bgroup,"bunching%d",(i+1));
    hh.writeBuffer(gid, bgroup, " ",  &bh[i-1]);
    sprintf(bgroup,"bunchingphase%d",(i+1));
    hh.writeBuffer(gid, bgroup,"rad",  &ph[i-1]);
  }

  if(do_global_stat) {
    gidsub=H5Gcreate(gid,"Global",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if (doEnergy){
      hh.writeSingleNode(gidsub,"energy"," ", &tgavg);
      hh.writeSingleNode(gidsub,"energyspread"," ", &tgsig);
    }
    if (doSpatial){
      hh.writeSingleNode(gidsub,"xposition","m", &txavg);
      hh.writeSingleNode(gidsub,"xsize","m", &txsig);
      hh.writeSingleNode(gidsub,"yposition","m", &tyavg);
      hh.writeSingleNode(gidsub,"ysize","m", &tysig);
    }
    H5Gclose(gidsub);  
  }
}
