#ifndef __GENESIS_SORTING__
#define __GENESIS_SORTING__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <cctype>

#include "Particle.h"
#include <mpi.h>

using namespace std;


class Sorting{
 public:
  Sorting();
  virtual ~Sorting();

  void init(int,int, bool, bool);
  void configure(double,double,double,double,double,double, bool);
  void globalSort(vector< vector< Particle > > *);

  int sort(vector< vector< Particle > > *);


 private:
  void update_stats(unsigned long long &, unsigned long long &);
  void globalSort_completion_msg(void);
  void fillPushVectors(vector< vector< Particle > > *);
  void localSort(vector< vector< Particle > > *);
  int centerShift(vector< vector< Particle > > *);
  void send(int, vector<double> *);
  void recv(int, vector<vector< Particle > >*, vector<double> *);
  
  int rank,size;

  double s0,slen,sendmin,sendmax,keepmin,keepmax;


  double reflen;
  bool doshift,dosort,globalframe;
  vector<double> pushforward,pushbackward;

  unsigned long long stat_max_xfer_size;
  int sort_round; // was introduced for generation of filenames for debug info
};



#endif
