#ifndef __GENESIS_COLLECTIVE__
#define __GENESIS_COLLECTIVE__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>

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
   void clearWake();
   void apply(Beam *,Undulator *, double );
   void update(Beam *, double);
   void forceUpdate();
   [[nodiscard]] bool hasWakeDefined() const;

private:
   bool transient,hasWake,needsUpdate;
   double ztrans,radius;
   double ds,dscur;
   unsigned int ns;
   int size,rank,ncur;
   double *wakeext, *wakeres, *wakegeo, *wakerou, *wake, *current, *dcurrent;
   std::vector<double> wakeint;
   std::vector<double> cur;
   std::vector<int> count;

   void resize_and_zero(std::vector<double>&, size_t);
   void resize_and_zero_i(vector<int>& v, size_t n);
};

inline bool Collective::hasWakeDefined() const{
    return hasWake;
}
inline void Collective::forceUpdate()
{
  needsUpdate=true;
}
#endif
