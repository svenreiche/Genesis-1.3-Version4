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
