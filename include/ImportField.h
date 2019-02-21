#ifndef __GENESIS_IMPORTFIELD__
#define __GENESIS_IMPORTFIELD__

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
#include "Lattice.h"
#include "Field.h"



using namespace std;

extern const double ce;

class ImportField : public StringProcessing{
 public:
   ImportField();
   virtual ~ImportField();
   bool init(int, int, map<string,string> *,vector<Field *> *,Setup *, Time *);

 private:
   void usage();
   string file;
   int harm;
   double offset;
   bool dotime;
};


#endif
