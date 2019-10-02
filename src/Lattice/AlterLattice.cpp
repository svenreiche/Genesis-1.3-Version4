#include "AlterLattice.h"

AlterLattice::AlterLattice()
{
  zmatch=0;
  element = "";
  field   = "";
  value   = 0;
  instance= 0;
  resolve = false;
  add     = true;
}

AlterLattice::~AlterLattice(){}

void AlterLattice::usage(){

  cout << "List of keywords for Lattice" << endl;
  cout << "&lattice" << endl;
  cout << " double zmatch = 0" << endl;
  cout << " string element = <empty string>" << endl;
  cout << " string field   = <empty string>" << endl;
  cout << " double value   = 0 / reference"  << endl;
  cout << " int instance   = 0 " << endl;
  cout << " bool resolvePeriod = false" << endl;
  cout << " bool add = true" << endl;
   cout << "&end" << endl << endl;
  return;
}


bool AlterLattice::init(int inrank, int insize, map<string,string> *arg, Lattice *lat, Setup *setup)
{

  rank=inrank;
  size=insize;
  zmatch=0;
  element = "";
  field   = "";
  value   = 0;
  instance= 0;
  resolve = false;
  add     = true;
 
  map<string,string>::iterator end=arg->end();

  if (arg->find("zmatch")!=end)  {zmatch= atof(arg->at("zmatch").c_str());  arg->erase(arg->find("zmatch"));}
  if (arg->find("value")!=end)   {value = atof(arg->at("value").c_str());  arg->erase(arg->find("value"));}
  if (arg->find("element")!=end) {element= arg->at("element");  arg->erase(arg->find("element"));}
  if (arg->find("field")!=end)   {field= arg->at("field");  arg->erase(arg->find("field"));}
  if (arg->find("instance")!=end){instance= atoi(arg->at("instance").c_str());  arg->erase(arg->find("instance"));}
  if (arg->find("resolvePeriod")!=end)  {resolve= atob(arg->at("resolvePeriod"));  arg->erase(arg->find("resolvePeriod"));}
  if (arg->find("add")!=end)     {add  = atob(arg->at("add"));    arg->erase(arg->find("add"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &lattice" << endl; this->usage();}
    return false;
  }
  
  if (zmatch>0) {
    lat->match(rank, zmatch, setup->getReferenceEnergy());
  }

  if (element!=""){
    bool val= lat->alterElement(element,field,value,instance,add);
    if (!val) {
      if (rank == 0 ) {
	cout << "*** Error: Input element and field does not match any supported elements" << endl;
      }
      return false;
    }
  }

  return true;


}

