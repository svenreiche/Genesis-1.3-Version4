#ifndef __GENESIS_PARSELATTICE__
#define __GENESIS_PARSELATTICE__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>
#include <algorithm>
#include <stdlib.h>

#include "StringProcessing.h"

using namespace std;

struct LatticeLayout {
      string key;
      double z;
      double l;
      double zoff;
      int ref;
};



class ParseLattice : public StringProcessing {
public:
    ParseLattice();
    virtual ~ParseLattice();
    bool parse(string, string, int);
    void generateLattice(double dz);


private:
    void parseArguments(vector<string> &, string &);
    bool dereference(const string key, const string type, int recursion);
    bool unroll(const string, int recursion);
    int checkMultiplier(string &);
    double checkReference(string &);

    // the default values of beamline elements
    const map<string, string> cdrift = {{"l", "0"}};
    const map<string, string> cquad = {{"l",  "0"},
                                       {"k1", "0"},
                                       {"dx", "0"},
                                       {"dy", "0"}};
    const map<string, string> ccorr = {{"l",  "0"},
                                       {"cx", "0"},
                                       {"cy", "0"}};
    const map<string, string> cphase = {{"l",   "0"},
                                        {"phi", "0"}};
    const map<string, string> cdelay = {{"l",     "0"},
                                        {"delay", "0"},
                                        {"lb",    "0"},
                                        {"ld",    "0"}};
    const map<string, string> cmarker = {{"sort",      "0"},
                                         {"stop",      "0"},
                                         {"dumpfield", "0"},
                                         {"dumpbeam",  "0"}};
    const map<string, string> cund = {{"l",       "0"},
                                      {"lambdau", "0"},
                                      {"aw",      "0"},
                                      {"nwig",    "1"},
                                      {"kx",      "0"},
                                      {"ky",      "0"},
                                      {"ax",      "0"},
                                      {"ay",      "0"},
                                      {"gradx",   "0"},
                                      {"grady",   "0"},
                                      {"helical", "0"}};
    map<string, map<string, string>> cele = {{"drif", cdrift},
                                             {"quad", cquad},
                                             {"corr", ccorr},
                                             {"phas", cphase},
                                             {"chic", cdelay},
                                             {"mark", cmarker},
                                             {"undu", cund}};

    map<string, map<string, string>> elements;    // final database of dereferenced elements
    map<string, vector<string>> rawlat;          // intermediate database of element
    map<string, vector<string>> lines;           // lines definition
    vector<LatticeLayout> beamline;              // selected beamline for line definition
    int iref;                                    // index for reference
    int rank;


   /*
   int findIndex(vector<string> *,string);
   bool unroll(int, int,int);
   bool resolve(int, int,int);
   int checkMultiplier(string *);
   bool checkResetPosition(string *,double *);
   vector<string> label,type,argument,sequence;
   vector<int> zref;
   vector<double> zoff;
   int refele;
   */
};


#endif
