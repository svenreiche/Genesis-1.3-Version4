#include "Wake.h"
#include "Beam.h"

Wake::Wake()
{
  radius=2.5e-3;
  conductivity=0;
  relaxation=0;
  roundpipe=true;
}
Wake::~Wake(){}


void Wake::usage(){

  cout << "List of keywords for Wake" << endl;
  cout << "&wake" << endl;
  cout << " double radius = 2.5e-3" << endl;
  cout << " bool   roundpipe   = true" << endl;
  cout << " string material  = <empty>" << endl;
  cout << " double conductivity = 0e-6" << endl;
  cout << " double relaxation  = 0e-6" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool SponRad::init(int rank, int size, map<string,string> *arg,  Beam *beam)
{

  string material='';
  map<string,string>::iterator end=arg->end();

  if (arg->find("radius")!=end) {radius = atof(arg->at("radius").c_str());  arg->erase(arg->find("radius"));}
  if (arg->find("conductivity")!=end) {conductivity= atof(arg->at("conductivity").c_str());  arg->erase(arg->find("conductivity"));}
  if (arg->find("relaxation")!=end) {relaxation = atof(arg->at("relaxation").c_str());  arg->erase(arg->find("relaxation"));}
  if (arg->find("roundpipe")!=end)    {roundpipe    = atob(arg->at("roundpipe").c_str());  arg->erase(arg->find("roundpipe"));}
  if (arg->find("material")!=end){material = arg->at("material"); arg->erase(arg->find("material"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &wake" << endl; this->usage();}
    return false;
  }
  
  if ((material=="CU") || (material=="Cu") || (material=="cu")){
    conductivity=5.813e7;
    relaxation=8.1e-6;
  }

  if ((material=="AL") || (material=="Al") || (material=="al")){
    conductivity=3.571e7;
    relaxation=2.4e-6;
  }

  return true;
}
 
