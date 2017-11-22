#ifdef MPE
#include "mpe.h"
#endif


#ifndef __GEN_MPEPROFILING__
#define __GEN_MPEPROFILING__


#include <string>



using namespace std;


class MPEProfiling{
 public:
  MPEProfiling();
  virtual ~MPEProfiling();
  
  void init(int);
  void finilize();
  void logIO(bool,bool,string);
  void logCalc(bool,bool,string);
  void logComm(bool,string);
  void logLoading(bool,string);

  void logEvent(string);

 private:
  int event1a,event1b;
  int event2a,event2b;
  int event3a,event3b;
  int event4a,event4b;
  int event5a,event5b;
  int event6a,event6b;
  int event;
};

extern MPEProfiling mpe;

#endif


  
