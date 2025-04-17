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

class SeriesManager;

class LatticeParser : public StringProcessing {
 public:
   LatticeParser();
   virtual ~LatticeParser();
   bool parse(string,string, int, vector<Element *> &, SeriesManager *);

 private:
   Quadrupole *parseQuad(int,int,double);
   Drift *parseDrift(int,int,double);
   Marker *parseMarker(int,int,double);
   Chicane *parseChicane(int,int,double);
   Corrector *parseCorrector(int,int,double);
   ID *parseID(int,int,double, SeriesManager *);
   Phaseshifter *parsePhaseshifter(int,int,double, SeriesManager *);
   bool extractParameterValue(string, string, SeriesManager *, int, string, double *);

   int findIndex(vector<string> *,string);
   bool unroll(int, int,int);
   bool resolve(int, int,int);
   int checkMultiplier(string *);
   bool checkResetPosition(string *,double *);
   vector<string> label,type,argument,sequence;
   vector<int> zref;
   vector<double> zoff;
   int refele;

   void addSequence(const string& label, string argument, int rank, SeriesManager *sm);
};


#endif
