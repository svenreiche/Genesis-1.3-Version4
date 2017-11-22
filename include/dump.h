#ifndef __GENESIS_DUMP__
#define __GENESIS_DUMP__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <cctype>

#include "StringProcessing.h"
#include "Setup.h"



using namespace std;

class Beam;
class Field;

class Dump: public StringProcessing{
 public:
   Dump();
   virtual ~Dump();
   bool init(int, int, map<string,string> *, Setup *, Beam *beam, vector<Field *> *field);
 private:
   void usage();
   string dumpfield,dumpbeam;
   int rank,size;
};



#endif
