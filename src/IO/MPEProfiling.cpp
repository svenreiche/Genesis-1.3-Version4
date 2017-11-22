#include "MPEProfiling.h"


MPEProfiling mpe;

MPEProfiling::MPEProfiling(){}

MPEProfiling::~MPEProfiling(){}

void MPEProfiling::init(int rank)
{
#ifdef MPE
        MPE_Init_log();

        event1a = MPE_Log_get_event_number(); 
        event1b = MPE_Log_get_event_number(); 
        event2a = MPE_Log_get_event_number(); 
        event2b = MPE_Log_get_event_number(); 
        event3a = MPE_Log_get_event_number(); 
        event3b = MPE_Log_get_event_number(); 
        event4a = MPE_Log_get_event_number(); 
        event4b = MPE_Log_get_event_number(); 
        event5a = MPE_Log_get_event_number(); 
        event5b = MPE_Log_get_event_number(); 
        event6a = MPE_Log_get_event_number(); 
        event6b = MPE_Log_get_event_number(); 

        event = MPE_Log_get_event_number(); 


        if (rank == 0) {
	   MPE_Describe_state(event1a, event1b, "Calc-Beam", "red");
	   MPE_Describe_state(event6a, event6b, "Calc-Field", "orange");
	   MPE_Describe_state(event2a, event2b, "MPI-Comm",   "blue");
	   MPE_Describe_state(event3a, event3b, "IO-Coll.",    "green");
	   MPE_Describe_state(event5a, event5b, "IO-Ind.",    "dark green");
	   MPE_Describe_state(event4a, event4b, "Loading",    "brown");
          
           MPE_Describe_event(event,"Milestone","white");
        }
#endif
        return;

}

void MPEProfiling::logCalc(bool end, bool beam, string text){

#ifdef MPE
  int eventa=event1a;
  int eventb=event1b;
  if (beam==false){
    eventa=event6a;
    eventb=event6b;
  }

  if (end==false){ 
        MPE_Log_event(eventa, 0, text.c_str());
  } else {
        MPE_Log_event(eventb, 0, text.c_str());
  }
#endif
  return;
}

void MPEProfiling::logComm(bool end, string text){

#ifdef MPE
  if (end==false){ 
      MPE_Log_event(event2a, 0, text.c_str());
  } else {
      MPE_Log_event(event2b, 0, text.c_str());
  }
#endif
  return;
}

void MPEProfiling::logIO(bool end, bool collect, string text){

#ifdef MPE
  int eventa=event5a;
  int eventb=event5b;
  if (collect){
    eventa=event3a;
    eventb=event3b;
  }

  if (end==false){ 
      MPE_Log_event(eventa, 0, text.c_str());
  } else {
      MPE_Log_event(eventb, 0, text.c_str());
  }
#endif
  return;
}




void MPEProfiling::logLoading(bool end, string text){
#ifdef MPE
  if (end==false){ 
      MPE_Log_event(event4a, 0, text.c_str());
  } else {
      MPE_Log_event(event4b, 0, text.c_str());
  }
#endif
  return;
}


void MPEProfiling::logEvent(string text)
{
#ifdef MPE
  MPE_Log_event(event,0,text.c_str());
#endif
  return;
}


void MPEProfiling::finilize()
{
#ifdef MPE
  string file="MPElog";
  MPE_Finish_log(file.c_str());
#endif
  return;
}
