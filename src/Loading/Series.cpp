#include "Series.h"


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
    if (rank==0){ cout << "*** Error: Unknown elements in &profile_const" << endl; this->usage();}
    return "";
  }
  if ((label.size()<1)&&(rank==0)){
    cout << "*** Error: Label not defined in &profile_const" << endl; this->usage();
  }
  return label;
}

double SeriesConst::value()
{
  return c0;
}

void SeriesConst::usage(){
  cout << "List of keywords for SEQUENCE_CONST" << endl;
  cout << "&series_const" << endl;
  cout << " string label = <empty>" << endl;
  cout << " double c0 = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

