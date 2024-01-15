#include "Series.h"
#include <cmath>

// cross check if Series.cpp is used at all?


Series::Series()
{
}

Series::~Series()
{
  serie.clear();
}

bool Series::init(int rank, map<string,string> *arg,string element)
{

  SeriesBase *p;
  string label;

  if (element.compare("&sequence_const")==0){
    p=(SeriesBase *)new SeriesConst();
    label=p->init(rank,arg);
  } 
  if (element.compare("&sequence_power")==0){
    p=(SeriesBase *)new SeriesPower();
    label=p->init(rank,arg);
  } 

  if (element.compare("&sequence_random")==0){
    p=(SeriesBase *)new SeriesRandom();
    label=p->init(rank,arg);
  } 

  if (label.size()<1){
    return false;
  } else {
    serie[label]=p;
  }
  if (rank==0) {cout << "Adding sequence with label: " << label << endl;}
  return true;
}

bool Series::check(string label){
  if (label.size() < 1) { return true; }
  if (serie.find(label)!=serie.end()){  
    return true;
  }
  return false;
}


double Series::value(double val, string label)
{
  if ((label.size()<1)||(serie.find(label)==serie.end())){  
    return val;
  } else {
    return  serie[label]->value();
  }
}




//------------------------------------
// individual series/sequences


string SeriesConst::init(int rank, map<string,string>*arg)
{

  string label="";
  c0=0;
  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &sequence_const" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &sequence_const" << endl; this->usage();
  }
  return label;
}

double SeriesConst::value()
{
  return c0;
}

void SeriesConst::usage(){
  cout << "List of keywords for SEQUENCE_CONST" << endl;
  cout << "&sequence_const" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}


// sequence mostly for tapering

string SeriesPower::init(int rank, map<string,string>*arg)
{

  string label="";

  c0=0;
  dc=0;
  alpha=0;
  n0=1;
  icount=0;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("dc")!=end)   {dc    = atof(arg->at("dc").c_str());  arg->erase(arg->find("dc"));}
  if (arg->find("alpha")!=end){alpha = atof(arg->at("alpha").c_str());  arg->erase(arg->find("alpha"));}
  if (arg->find("n0")!=end)   {n0    = atoi(arg->at("n0").c_str());  arg->erase(arg->find("n0"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &sequence_power" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &sequence_power" << endl; this->usage();
  }
  return label;
}

double SeriesPower::value()
{
  icount++;
  if (icount <= n0) {
    return c0;
  }
  return c0+dc*pow(static_cast<double>(icount-n0),alpha);
}

void SeriesPower::usage(){
  cout << "List of keywords for SEQUENCE_POWER" << endl;
  cout << "&sequence_power" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double dc = 0" << endl;
  cout << " double alpha = 0" << endl;
  cout << " int n0 = 1" << endl;
  cout << "&end" << endl << endl;
  return;
}


// sequence mostly for undulator error

string SeriesRandom::init(int rank, map<string,string>*arg)
{

  string label="";

  c0=0;
  dc=0;
  seed=100;
  gauss=true;

  map<string,string>::iterator end=arg->end();

  if (arg->find("label")!=end){label = arg->at("label");  arg->erase(arg->find("label"));}
  if (arg->find("c0")!=end)   {c0    = atof(arg->at("c0").c_str());  arg->erase(arg->find("c0"));}
  if (arg->find("dc")!=end)   {dc    = atof(arg->at("dc").c_str());  arg->erase(arg->find("dc"));}
  if (arg->find("seed")!=end) {seed  = atof(arg->at("seed").c_str());  arg->erase(arg->find("seed"));}
  if (arg->find("normal")!=end){gauss= atob(arg->at("normal").c_str());  arg->erase(arg->find("normal"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &sequence_random" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &sequence_random" << endl; this->usage();
  }

  seq  = new RandomU (seed);
  return label;
}

double SeriesRandom::value()
{
  
  if (gauss){
    return c0+dc*(erf.value(2*seq->getElement()));
  }
  return c0+dc*(2*seq->getElement()-1);
}

void SeriesRandom::usage(){
  cout << "List of keywords for SEQUENCE_RANDOM" << endl;
  cout << "&sequence_random" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << " double dc = 0" << endl;
  cout << " int seed = 100" << endl;  
  cout << " bool normal = true" << endl;
  cout << "&end" << endl << endl;
  return;
}

