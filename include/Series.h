#ifndef __GENESIS__SERIES__
#define __GENESIS__SERIES__

#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <stdlib.h>

#include <mpi.h>
#include "StringProcessing.h"
#include "Sequence.h"
#include "RandomU.h"
#include "Inverfc.h"

using namespace std;


class SeriesBase: public StringProcessing{
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


class SeriesPower : public SeriesBase
{
 public:
  SeriesPower(){};
  ~SeriesPower(){};
  double value();
  string init(int, map<string,string> *);
  void usage();
 private:
  double c0,dc,alpha;
  int n0,icount;
};


class SeriesRandom : public SeriesBase
{
 public:
  SeriesRandom(){};
  ~SeriesRandom(){};
  double value();
  string init(int, map<string,string> *);
  void usage();
 private:
  double c0,dc,seed;
  bool gauss;
  RandomU *seq;
  Inverfc erf;
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
