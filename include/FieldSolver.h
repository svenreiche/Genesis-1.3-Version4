#ifndef __GENESIS_FIELDSOLVER__
#define __GENESIS_FIELDSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>


class Field;
class Beam;

#include "Particle.h"
#include "Undulator.h"


using namespace std;


class FieldSolver{
 public:
    virtual ~FieldSolver() {};
    virtual void init(double,double,double,int) = 0;
    virtual void advance(double, Field *, Beam *, Undulator *) = 0;
    virtual void initSourceFilter(double,double,double,bool) = 0;
};


#endif
