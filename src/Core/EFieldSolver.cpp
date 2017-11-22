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
