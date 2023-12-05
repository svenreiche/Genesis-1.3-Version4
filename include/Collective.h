#ifndef __GENESIS_COLLECTIVE__
#define __GENESIS_COLLECTIVE__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <math.h>

#include <mpi.h>

class Beam;
class Undulator;

extern bool MPISingle; 
extern const double ce;   

using namespace std;


class Collective{
public:
   Collective();
   virtual ~Collective();
   //   void initWake(unsigned int, double, double *, double *, double *, double, double, bool);
   void initWake(unsigned int, unsigned int, double, double *, double *, double *, double *, double, double, bool);
   void apply(Beam *,Undulator *, double );
   void update(Beam *, double);
   void forceUpdate();

private:
   bool transient,hasWake,needsUpdate;
   double ztrans,radius;
   double ds,dscur;
   unsigned int ns;
   int size,rank,ncur;
   double *wakeext, *wakeint, *wakeres, *wakegeo, *wakerou, *wake, *current, *dcurrent;
   // double *cur;
   std::vector<double> cur;
   int *count;
};


inline void Collective::forceUpdate()
{
  needsUpdate=true;
}
#endif
