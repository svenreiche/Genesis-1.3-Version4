#include "dump.h"
#include "writeFieldHDF5.h"
#include "writeBeamHDF5.h"

Dump::Dump()
{
  dumpfield="";
  dumpbeam ="";
}

Dump::~Dump(){}

void Dump::usage(){

  cout << "List of keywords for WRITE" << endl;
  cout << "&write" << endl;
  cout << " string field = <empty>" << endl;
  cout << " string beam  = <empty>" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool Dump::init(int inrank, int insize, map<string,string> *arg, Setup *setup, Beam *beam, vector<Field *> *field)
{

  rank=inrank;
  size=insize;

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("field")!=end){dumpfield=arg->at("field"); arg->erase(arg->find("field"));}
  if (arg->find("beam")!=end) {dumpbeam =arg->at("beam");  arg->erase(arg->find("beam"));}
  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &write" << endl; this->usage();}
    return false;
  }

  
  if (dumpfield.size()>0){
   WriteFieldHDF5 dump;
   dump.write(dumpfield,field);
  }
  if (dumpbeam.size()>0){
   WriteBeamHDF5 dump;
   dump.write(dumpbeam,beam);
  }

  return true;


}


