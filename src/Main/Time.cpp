#include "Time.h"

Time::Time()
{
  dotime=false;
  doscan=false;
  s0=0;
  slen=0;
  ds=1;
  sample=1;
  nslice=1;
  ns_node=1;
  noff_node=0;
  initialized=false;
}

Time::~Time(){}

void Time::usage(){

  cout << "List of keywords for TIME" << endl;
  cout << "&time" << endl;
  cout << " double s0 = 0" << endl;
  cout << " double slen   = 0" << endl;
  cout << " int sample = 1" << endl;
  cout << " bool time = true" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool Time::init(int inrank, int insize, map<string,string> *arg, Setup *setup)
{

  rank=inrank;
  size=insize;
  dotime=true;

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("s0")!=end)  {s0   = atof(arg->at("s0").c_str());  arg->erase(arg->find("s0"));}
  if (arg->find("slen")!=end){slen = atof(arg->at("slen").c_str()); arg->erase(arg->find("slen"));}
  if (arg->find("sample")!=end){sample = atoi(arg->at("sample").c_str()); arg->erase(arg->find("sample"));}
  if (arg->find("time")!=end){dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}
  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &time" << endl; this->usage();}
    return false;
  }

  initialized=true;  
  this->finishInit(setup);
  return true;


}


void Time::finishInit(Setup *setup)
{
  if (!initialized){return;}

  doscan=!dotime;

  ds=setup->getReferenceLength()*sample;
  nslice=static_cast<int> (round(slen/ds));
  if (nslice < size) { nslice = size ; }

  ns_node=nslice/size;
  if ((nslice % size)!=0) {ns_node++;}
  nslice=ns_node*size;
  noff_node=ns_node*rank;
  // needs to adjust the time window for harmonic conversions
  slen=ds*nslice;
  
  if (rank==0){
    cout << "Setting up time window of " << slen*1e6 << " microns with " << nslice <<" sample points..." << endl;
  }
}




int Time::getPosition(vector<double> *s)
{
  if (nslice<1){nslice=1;} 
  s->resize(nslice);
  for (int i=0;i <nslice; i++){
    s->at(i)=s0+i*ds;
  }
  return nslice;
}
