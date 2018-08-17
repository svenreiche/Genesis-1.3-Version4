#ifndef __GENESIS_LATTICEPARSER__
#define __GENESIS_LATTICEPARSER__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>

#include "StringProcessing.h"
#include "LatticeElements.h"


using namespace std;



class LatticeParser : public StringProcessing {
 public:
   LatticeParser();
   virtual ~LatticeParser();
   bool parse(string,string, int, vector<Element *> &,bool);
 private:

   Quadrupole *parseQuad(int,int,double);
   Drift *parseDrift(int,int,double);
   Marker *parseMarker(int,int,double);
   Chicane *parseChicane(int,int,double);
   Corrector *parseCorrector(int,int,double);
   ID *parseID(int,int,double);
   Phaseshifter *parsePhaseshifter(int,int,double);

   int findIndex(vector<string> *,string);
   bool unroll(int, int,int);
   bool resolve(int, int,int);
   int checkMultiplier(string *);
   bool checkResetPosition(string *);
   vector<string> label,type,argument,sequence;
   vector<int> zref;
   int refele;
};


#endif
