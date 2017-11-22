#ifndef __GENESIS_SPONRAD__
#define __GENESIS_SPONRAD__

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


class SponRad: public StringProcessing{
 public:
   SponRad();
   virtual ~SponRad();
   bool init(int,int,map<string,string> *, Beam *);


 private:
   void usage();
   bool doLoss,doSpread;
   int seed;

};

#endif
