#ifndef __GENESIS_ALTERLATTICE__
#define __GENESIS_ALTERLATTICE__

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <map>
#include <stdlib.h>

#include "Lattice.h"
#include "Setup.h"
#include "Series.h"
#include "SeriesManager.h"

class Setup;
class Lattice;

using namespace std;

class AlterLattice:  public StringProcessing{
 public:
    AlterLattice();
    virtual ~AlterLattice();
    bool init(int, int, map<string,string> *, Lattice *, Setup *, SeriesManager *);
    double getFieldError() const;
    bool getOrbitError() const;
    unsigned long getSeed() const;
 private:
    void usage();
    string element,field;
    double value;
    int instance;
    unsigned long iseed;
    bool resolve,add, orbiterror;
    double zmatch, fielderror;
    int rank{}, size{};
    string valueref;
};
inline double AlterLattice::getFieldError() const{return fielderror;}
inline bool AlterLattice::getOrbitError() const{return orbiterror;}
inline unsigned long AlterLattice::getSeed() const{return iseed;}
#endif
