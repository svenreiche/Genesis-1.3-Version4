#include "dump.h"
#include "writeFieldHDF5.h"
#include "writeBeamHDF5.h"

Dump::Dump(){}
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
  string dumpfield, dumpbeam;
  map<string,string>::iterator end=arg->end();

  if (arg->find("field")!=end){dumpfield=arg->at("field"); arg->erase(arg->find("field"));}
  if (arg->find("beam")!=end) {dumpbeam =arg->at("beam");  arg->erase(arg->find("beam"));}
  if (arg->size()!=0){
    if (inrank==0){ cout << "*** Error: Unknown elements in &write" << endl; this->usage();}
    return(false);
  }

  
  if (dumpfield.size()>0){
   WriteFieldHDF5 dump;
   string completefn;
   setup->RootName_to_FileName(&completefn, &dumpfield);
   if(!dump.write(completefn,field)) {
     if(inrank==0) {
       cout << "   write operation was not successful!" << endl;
     }
     return(false);
   }
  }
  if (dumpbeam.size()>0){
   WriteBeamHDF5 dump;

   // propagate beam dump settings before initiating writing procedure
   beam->setWriteFilter(setup->BWF_get_enabled(), setup->BWF_get_from(), setup->BWF_get_to(), setup->BWF_get_inc());

   string completefn;
   setup->RootName_to_FileName(&completefn, &dumpbeam);
   dump.write(completefn,beam);
  }

  return(true);
}
