#ifndef __GENESIS_LOADFIELD__
#define __GENESIS_LOADFIELD__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <complex>

#include "StringProcessing.h"
#include "Setup.h"
#include "Time.h"
#include "Profile.h"
#include "GaussHermite.h"
#include "Field.h"

using namespace std;



class LoadField : public StringProcessing{
 public:
   LoadField();
   virtual ~LoadField();
   bool init(int, int, map<string,string> *,vector<Field *> *, Setup *, Time *, Profile *);

 private:
   void usage();
   double lambda,power,phase,z0,w0;
   string lambdaref,powerref,phaseref,z0ref,w0ref;
   double dgrid,xcen,ycen,xangle,yangle;
   int ngrid;
   int harm,nx,ny;
   bool add;
};


#endif
