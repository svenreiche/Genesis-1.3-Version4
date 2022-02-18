#ifndef __GENESIS_PARSER__
#define __GENESIS_PARSER__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
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

   bool fail(void);

 private:

   int rank;
   istringstream input;
   bool fillMap(string *,map<string,string> *);

   bool is_parse_error;
};


#endif
