#ifndef __GENESIS_TIME__
#define __GENESIS_TIME__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <cctype>

#include "StringProcessing.h"
#include "Setup.h"



using namespace std;



class Time: public StringProcessing{
 public:
   Time();
   virtual ~Time();
   bool init(int, int, map<string,string> *, Setup *);
   void finishInit(Setup *);

   void set(double,double,int,bool);

   int getPosition(vector<double> *);
   double getSampleRate();
   double getTimeWindowStart();
   double getTimeWindowLength();
   bool isTime();
   bool isScan();
   int getNodeOffset();
   int getNodeNSlice();
   bool isInitialized();

   void setSampleRate(double);
   void setTimeWindowStart(double);
   void setTimeWindowLength(double);


 private:
   void usage();
   bool dotime, doscan;
   double s0,slen,ds;
   int rank, size,sample;
   int ns_node,noff_node,nslice;
   bool initialized;
};

inline void Time::setSampleRate(double in){ sample=static_cast<int>(round(in)); }
inline void Time::setTimeWindowStart(double in) { s0=in; }
inline void Time::setTimeWindowLength(double in) { slen=in; }


inline double Time::getTimeWindowStart(){ return s0; }
inline double Time::getTimeWindowLength(){ return slen; }
inline bool Time::isTime(){ return dotime; }
inline bool Time::isScan(){ return doscan; }
inline int Time::getNodeOffset(){ return noff_node; }
inline int Time::getNodeNSlice(){ return ns_node; }
inline bool Time::isInitialized(){ return initialized; }

inline double Time::getSampleRate()
{
  if (dotime) { return sample;} else { return 1.; }
  //return sample;
}


#endif
