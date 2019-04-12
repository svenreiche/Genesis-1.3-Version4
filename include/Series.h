#ifndef __GENESIS__SERIES__
#define __GENESIS__SERIES__

#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <stdlib.h>

#include "mpi.h"
#include "StringProcessing.h"


using namespace std;


class SeriesBase{
 public:
  SeriesBase(){};
  ~SeriesBase(){};
  virtual double value()=0;
  virtual string init(int, map<string,string> *)=0;
  virtual void usage()=0;
};

class SeriesConst : public SeriesBase
{
 public:
  SeriesConst(){};
  ~SeriesConst(){};
  double value();
  string init(int, map<string,string> *);
  void usage();
 private:
  double c0;
};


//-------------------------

class Series{
 public:
   Series();
   ~Series();
   bool init(int, map<string,string> *, string); 
   bool check(string);
   double value(double, string);
 private:
   map<string,SeriesBase *> serie;
   
};



#endif
