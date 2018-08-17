#include <sstream>
#include "Setup.h"

Setup::Setup()
{
  rootname="";
  lattice="";
  beamline="";
  partfile="";
  fieldfile="";
  one4one=false;
  shotnoise=true;
  nbins=4;
  npart=8192;
  gamma0=5800/0.511;
  lambda0=1e-10;
  delz=0.015; 
  seed=123456789;
  runcount = 0 ;  // count of runs in conjunction of calls of altersetup 
}

Setup::~Setup(){}

void Setup::usage(){

  cout << "List of keywords for SETUP" << endl;
  cout << "&setup" << endl;
  cout << " string rootname = <empty>" << endl;
  cout << " string lattice = <empty>" << endl;
  cout << " string beamline = <empty>" << endl;
  cout << " string partfile = <empty>" << endl;
  cout << " string fieldfile = <empty>" << endl;
  cout << " double gamma0 = 5800/0.511" << endl;
  cout << " double lambda0 = 1e-10" << endl;
  cout << " double delz = 0.015" << endl;
  cout << " int seed = 123456789" << endl;
  cout << " int npart = 8192" << endl;
  cout << " int nbins = 4" << endl;
  cout << " bool one4one = false" << endl;
  cout << " bool shotnoise = true" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool Setup::init(int inrank, map<string,string> *arg, Lattice *lat,string latstring,bool streaming)
{

  rank=inrank;
  map<string,string>::iterator end=arg->end();

  if (arg->find("rootname")!=end){rootname = arg->at("rootname"); arg->erase(arg->find("rootname"));}
  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice");  arg->erase(arg->find("lattice"));}
  if (arg->find("beamline")!=end){beamline = arg->at("beamline"); arg->erase(arg->find("beamline"));}
  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice");  arg->erase(arg->find("lattice"));}
  if (arg->find("partfile")!=end){partfile   = arg->at("partfile");  arg->erase(arg->find("partfile"));}
  if (arg->find("fieldfile")!=end){fieldfile = arg->at("fieldfile");  arg->erase(arg->find("fieldfile"));}
  if (arg->find("gamma0")!=end)  {gamma0   = atof(arg->at("gamma0").c_str());  arg->erase(arg->find("gamma0"));}
  if (arg->find("lambda0")!=end) {lambda0  = atof(arg->at("lambda0").c_str()); arg->erase(arg->find("lambda0"));}
  if (arg->find("delz")!=end)    {delz     = atof(arg->at("delz").c_str());  arg->erase(arg->find("delz"));}
  if (arg->find("seed")!=end)    {seed     = atoi(arg->at("seed").c_str());  arg->erase(arg->find("seed"));}
  if (arg->find("one4one")!=end) {one4one  = atob(arg->at("one4one"));  arg->erase(arg->find("one4one"));}
  if (arg->find("npart")!=end)    {npart  = atoi(arg->at("npart").c_str());  arg->erase(arg->find("npart"));}
  if (arg->find("nbins")!=end)    {nbins  = atoi(arg->at("nbins").c_str());  arg->erase(arg->find("nbins"));}
  if (arg->find("shotnoise")!=end){shotnoise  = atob(arg->at("shotnoise"));  arg->erase(arg->find("shotnoise"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &setup" << endl; this->usage();}
    return false;
  }

  if (streaming) {
    lattice="streaming";
    lat->parse(latstring,beamline,rank,streaming);
  }else{
    lat->parse(lattice,beamline,rank,streaming);
  }
  return true;
}



bool Setup::getRootName(string *filename)
{
  if (rootname.size()<1){
    return false;
  }
  *filename=rootname;
  if (runcount>0) {
    stringstream ss;
    ss << ".Run" << (runcount+1) ;
    *filename+=ss.str();
  }
  return true; 


}
