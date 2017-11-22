#ifndef __GENESIS__PROFILE__
#define __GENESIS__PROFILE__

#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <stdlib.h>

#include "mpi.h"
#include "StringProcessing.h"
#include "HDF5base.h"

using namespace std;


class ProfileBase{
 public:
  ProfileBase(){};
  ~ProfileBase(){};
  virtual double value(double)=0;
  virtual string init(int, map<string,string> *)=0;
  virtual void usage()=0;
};

class ProfileConst : public ProfileBase
{
 public:
  ProfileConst(){};
  ~ProfileConst(){};
  double value(double);
  string init(int, map<string,string> *);
  void usage();
 private:
  double c0;
};


class ProfilePolynom : public ProfileBase
{
 public:
  ProfilePolynom(){};
  ~ProfilePolynom(){};
  string init(int, map<string,string> *);
  double value(double);
  void usage();
 private:
  vector<double> c ; 
};


class ProfileStep : public ProfileBase
{
 public:
  ProfileStep(){};
  ~ProfileStep(){};
  string init(int, map<string,string> *);
  double value(double);
  void usage();
 private:
  double c0,send,sstart;
};

class ProfileGauss : public ProfileBase
{
 public:
  ProfileGauss(){};
  ~ProfileGauss(){};
  string init(int, map<string,string> *);
  double value(double);
  void usage();
 private:
  double c0,s0,sig;   
};



class ProfileFile : public ProfileBase, StringProcessing, HDF5Base
{
 public:
  ProfileFile(){};
  ~ProfileFile(){};
  string init(int, map<string,string> *);
  double value(double);
  void usage();
 private:
  string xdataset,ydataset;   
  bool isTime,revert;
  vector<double> xdat,ydat;
};


//-------------------------

class Profile{
 public:
   Profile();
   ~Profile();
   bool init(int, map<string,string> *, string); 
   double value(double,double,string);
 private:
   map<string,ProfileBase *> prof;
   
};



#endif
