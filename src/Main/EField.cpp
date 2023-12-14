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

bool EField::init(int rank, int size, map<string,string> *arg,  Beam *beam, Setup *setup, Time *time)
{

  bool dotime=time->isTime();                  // check for time simulation
  lambda=setup->getReferenceLength();
  longrange=false;
  redLorentz=true;

  auto end=arg->end();

  if (arg->find("rmax")!=end)  {rmax  = atof(arg->at("rmax").c_str()); arg->erase(arg->find("rmax"));}
  if (arg->find("ngrid")!=end) {ngrid = atoi(arg->at("ngrid").c_str());  arg->erase(arg->find("ngrid"));}
  if (arg->find("nz")!=end)    {nz    = atoi(arg->at("nz").c_str());  arg->erase(arg->find("nz"));}
  if (arg->find("nphi")!=end)  {nphi  = atoi(arg->at("nphi").c_str());  arg->erase(arg->find("nphi"));}
  if (arg->find("longrange")!=end)  {longrange = atob(arg->at("longrange").c_str()); arg->erase(arg->find("longrange"));}

  if (!arg->empty()){
    if (rank==0){ cout << "*** Error: Unknown elements in &efield" << endl; this->usage();}
    return false;
  }

  beam->initEField(rmax,ngrid,nz,nphi,lambda,longrange);
  return true;
}
 

 
