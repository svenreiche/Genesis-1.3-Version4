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
