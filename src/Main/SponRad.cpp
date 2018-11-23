#include "SponRad.h"
#include "Beam.h"

SponRad::SponRad(){}
SponRad::~SponRad(){}

void SponRad::usage(){

  cout << "List of keywords for SponRad" << endl;
  cout << "&sponrad" << endl;
  cout << " int seed = 1234 " << endl;
  cout << " bool doLoss   = false" << endl;
  cout << " bool doSpread = false" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool SponRad::init(int rank, int size, map<string,string> *arg,  Beam *beam)
{


  map<string,string>::iterator end=arg->end();

  if (arg->find("seed")!=end) {seed = atoi(arg->at("seed").c_str());  arg->erase(arg->find("seed"));}
  if (arg->find("doLoss")!=end)    {doLoss    = atob(arg->at("doLoss").c_str());  arg->erase(arg->find("doLoss"));}
  if (arg->find("doSpread")!=end)  {doSpread  = atob(arg->at("doSpread").c_str());  arg->erase(arg->find("doSpread"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &sponrad" << endl; this->usage();}
    return false;
  }

  beam->initIncoherent(seed,rank,doLoss,doSpread);

  return true;
}
 


