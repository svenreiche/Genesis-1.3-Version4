#ifndef __GENESIS_OPTICS__
#define __GENESIS_OPTICS__

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <map>




using namespace std;

class Optics{
 public:
   Optics();
   virtual ~Optics();
   void init();
   void addElement(double, double, double, double);
   bool match(double *, double *, double *, double *, double *, double *);
 private:
   double Dx11,Dx12,Dx21,Dx22,Dy11,Dy12,Dy21,Dy22;
   double M11,M12,M21,M22;
   void getDrift(double);
   void getQuad(double, double);
   void MatMult(bool);
};


#endif
