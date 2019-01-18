#ifndef __GENESIS_WAKE__
#define __GENESIS_WAKE__

#include <vector>
#include <iostream>
#include <string>
#include <complex>


#include "StringProcessing.h"
#include "Setup.h"
#include "Time.h"
#include "Particle.h"

class Beam;

using namespace std;

extern const double vacimp; 
extern const double ce; 

class Wake: public StringProcessing{
 public:
   Wake();
   virtual ~Wake();
   bool init(int,int,map<string,string> *, Time *, Setup *, Beam *);


 private:   
   void usage();
   void singleWakeResistive(int);


   double radius, relaxation,conductivity,ztrans;
   bool roundpipe,transient, hasWake;

   unsigned int ns;
   double slen,ds;
   double *wakeres;

};

#endif
