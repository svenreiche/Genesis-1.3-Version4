#include "AlterLattice.h"

AlterLattice::AlterLattice()
{
  zmatch=0;
  err_aw=0;
  err_ax=0;
  err_ay=0;
  err_qx=0;
  err_qy=0;
  nlat=1;
}

AlterLattice::~AlterLattice(){}

void AlterLattice::usage(){

  cout << "List of keywords for Lattice" << endl;
  cout << "&lattice" << endl;
  cout << " double err_aw = 0" << endl;
  cout << " double err_ax = 0" << endl;
  cout << " double err_ay = 0" << endl;
  cout << " double err_qx = 0" << endl;
  cout << " double err_qy = 0" << endl;
  cout << " double zmatch = 0" << endl;
  cout << " int nlat = 1" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool AlterLattice::init(int inrank, int insize, map<string,string> *arg, Lattice *lat, Setup *setup)
{

  rank=inrank;
  size=insize;

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("zmatch")!=end)  {zmatch= atof(arg->at("zmatch").c_str());  arg->erase(arg->find("zmatch"));}
  if (arg->find("err_aw")!=end)  {err_aw= atof(arg->at("err_aw").c_str());  arg->erase(arg->find("err_aw"));}
  if (arg->find("err_ax")!=end)  {err_aw= atof(arg->at("err_ax").c_str());  arg->erase(arg->find("err_ax"));}
  if (arg->find("err_ay")!=end)  {err_aw= atof(arg->at("err_ay").c_str());  arg->erase(arg->find("err_ay"));}
  if (arg->find("err_qx")!=end)  {err_aw= atof(arg->at("err_qx").c_str());  arg->erase(arg->find("err_qx"));}
  if (arg->find("err_qy")!=end)  {err_aw= atof(arg->at("err_qy").c_str());  arg->erase(arg->find("err_qy"));}
  if (arg->find("nlat")!=end)    {nlat  = atof(arg->at("nlat").c_str());    arg->erase(arg->find("nlat"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &lattice" << endl; this->usage();}
    return false;
  }
  
  if (zmatch>0) {
    lat->match(rank, zmatch, setup->getReferenceEnergy());
  }

  return true;


}

