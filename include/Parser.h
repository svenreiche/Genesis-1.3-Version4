#ifndef __GENESIS_PARSER__
#define __GENESIS_PARSER__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>

#include "StringProcessing.h"


using namespace std;



class Parser : public StringProcessing {
 public:
   Parser();
   virtual ~Parser();
   bool open(string, int);
   bool parse(string *,map<string,string> *);

 private:

   int rank;
   ifstream fin;
   bool fillMap(string *,map<string,string> *);

};


#endif
