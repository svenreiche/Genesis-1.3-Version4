#ifndef __GENESIS_INCOHERENT__
#define __GENESIS_INCOHERENT__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <math.h>

#include "Undulator.h"
#include "Particle.h"
#include "Sequence.h"
#include "RandomU.h"

class Beam;

using namespace std;

extern const double vacimp;
extern const double eev;



class Incoherent{
 public:
   Incoherent();
   virtual ~Incoherent();
   void init(int,bool,bool);
   void apply(Beam *,Undulator *und, double );

 private:
   // Builds one generator per slice held on this core, each keyed on the index
   // of its slice in the full time window. Called on first use because
   // &sponrad may be parsed before the beam has been loaded.
   void ensureStreams(const Beam *);

   bool doLoss,doSpread;
   unsigned int base;
   int noff;
   vector<RandomU> sran;
};

#endif
