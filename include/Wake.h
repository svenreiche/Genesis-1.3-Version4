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


class Wake: public StringProcessing{
 public:
   Wake();
   virtual ~Wake();
   bool init(int,int,map<string,string> *, Beam *);


 private:
   void usage();
   double radius, relaxation,conductivity;
   bool roundpipe;


};

#endif
