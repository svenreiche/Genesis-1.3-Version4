#ifndef __GENESIS_COLLECTIVE__
#define __GENESIS_COLLECTIVE__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <math.h>

extern const double vacimp;
extern const double ce;

using namespace std;


class Collective{
 public:
   Collective();
   virtual ~Collective();

   void WakeRes();
   void WakeGeo();
   void WakeRou();

 private:
   bool doWakes,doSpaceCharge,doCSR;

   vector<double> current;  // holds current profile
   double ds;               // sample rate
   double ns;               // sample points

   // for wakefields:
   vector<double> wakeres,wakegeo,wakerou;
   
   double r,sigma,tau;  // beam pipe radius,, conductivity and relaxation time
   double gap,lgap;    // gap and the normalization length for the gap.
   double hrough,lrough; // amplitude and periodlength of roughness
   bool roundpipe;



};

#endif
