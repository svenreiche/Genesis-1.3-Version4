#include "EField.h"
#include "Beam.h"

EField::EField()= default;
EField::~EField()= default;

void EField::usage(){

  cout << "List of keywords for EFIELD" << endl;
  cout << "&efield" << endl;
  cout << " double rmax = 0 " << endl;
  cout << " int nz   = 0" << endl;
  cout << " int nphi = 0" << endl;
  cout << " int ngrid = 100" << endl;
  cout << " bool longrange = False" << endl;
  cout << "&end" << endl << endl;
}

bool EField::init(int rank, map<string,string> *arg,  Beam *beam, Setup *setup)
{

  double lambda=setup->getReferenceLength();
  bool longrange=false;

  auto end=arg->end();

  if (arg->find("rmax")!=end)  {rmax  = strtod(arg->at("rmax").c_str(), nullptr); arg->erase(arg->find("rmax"));}
  if (arg->find("ngrid")!=end) {ngrid = strtol(arg->at("ngrid").c_str(), nullptr, 10);  arg->erase(arg->find("ngrid"));}
  if (arg->find("nz")!=end)    {nz    = strtol(arg->at("nz").c_str(), nullptr,10);  arg->erase(arg->find("nz"));}
  if (arg->find("nphi")!=end)  {nphi  = strtol(arg->at("nphi").c_str(),nullptr,10);  arg->erase(arg->find("nphi"));}
  if (arg->find("longrange")!=end)  {longrange = atob(arg->at("longrange").c_str()); arg->erase(arg->find("longrange"));}

  if (!arg->empty()){
    if (rank==0){ cout << "*** Error: Unknown elements in &efield" << endl; this->usage();}
    return false;
  }

  beam->initEField(rmax,ngrid,nz,nphi,lambda,longrange);
  return true;
}
 

 
