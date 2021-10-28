#ifndef __GENESIS__PROFILE__
#define __GENESIS__PROFILE__

#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <stdlib.h>

#include <mpi.h>
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
  vector<string> names;
  string xdataset,ydataset;   
  bool isTime,revert,autoassign;
 private:
  vector<double> xdat,ydat;
};

class ProfileFileMulti: public StringProcessing, HDF5Base /* cannot get interpolated value from this class: ProfileBase is not a base class */
{
public:
	ProfileFileMulti(){};
	~ProfileFileMulti(){};
	bool setup(int, map<string, string> *, map<string, ProfileBase *> *);
	void usage();
};

class ProfileInterpolator: public ProfileBase
{
public:
	ProfileInterpolator(string, vector<double> *, vector<double> *);
	~ProfileInterpolator();
	double value(double);
	string init(int, map<string,string> *);
	void usage();

private:
	vector<double> xdat_, ydat_;
	string label_;
};


//-------------------------

class Profile{
 public:
   Profile();
   ~Profile();
   bool init(int, map<string,string> *, string); 
   bool check(string);
   double value(double,double,string);
 private:
   map<string,ProfileBase *> prof;
   
};



#endif
