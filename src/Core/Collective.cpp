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


