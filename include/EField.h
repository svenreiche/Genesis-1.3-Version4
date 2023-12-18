#ifndef __GENESIS_EFIELD__
#define __GENESIS_EFIELD__

#include <vector>
#include <iostream>
#include <string>
#include <complex>


#include "StringProcessing.h"
#include "Setup.h"
#include "GenTime.h"
#include "Particle.h"

class Beam;

using namespace std;


class EField: public StringProcessing{
 public:
   EField();
   virtual ~EField();
   bool init(int,map<string,string> *, Beam *, Setup *);


 private:
   static void usage();
   double rmax{0};
   long nz{0},nphi{0},ngrid{100};
};

#endif
