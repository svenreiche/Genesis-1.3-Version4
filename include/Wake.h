#ifndef __GENESIS_WAKE__
#define __GENESIS_WAKE__

#include <vector>
#include <iostream>
#include <string>
#include <complex>


#include "StringProcessing.h"
#include "Setup.h"
#include "GenTime.h"
#include "Particle.h"
#include "GenProfile.h"

class Beam;

using namespace std;

extern const double vacimp; 
extern const double ce; 

class Wake: public StringProcessing{
 public:
   Wake();
   virtual ~Wake();
   bool init(int,int,map<string,string> *, Time *, Setup *, Beam *, Profile *);


 private:   
   void usage();
   void singleWakeResistive(int);
   void singleWakeGeometric(int);
   void singleWakeRoughness(int);
   void KernelRoughness(vector<complex<double> > *, complex<double>, complex<double>);
   double TrapIntegrateRoughness(vector< complex<double> > *, complex<double> , complex<double> , double);

   double radius, relaxation,conductivity,ztrans,gap,lgap,hrough,lrough,rrough;
   double loss;
   string lossref;
   bool roundpipe,transient, hasWake;

   unsigned int ns;
   double slen,ds;
   double *wakeres, *wakegeo, *wakerou, *wakeext, *wakeint;


};

#endif
