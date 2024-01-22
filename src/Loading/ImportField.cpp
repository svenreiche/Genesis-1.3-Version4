#include "ImportField.h"
#include "readFieldHDF5.h"


ImportField::ImportField() {}
ImportField::~ImportField(){}

void ImportField::usage(){

  cout << "List of keywords for IMPORTFIELD" << endl;
  cout << "&importfield" << endl;
  cout << " string file = <empty>" << endl;
  cout << " int harmonic = 1" << endl;
  cout << " bool time = true" << endl;
  cout << " bool replace = false" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool ImportField::init(int rank, int size, map<string,string> *arg, vector<Field *> *fieldin, Setup *setup, Time *time)
{
  // parameters and their defaults
  string file;
  int harm=1;
  bool dotime=true;
  bool force_replace=false;
  double offset=0.;

  double lambda=setup->getReferenceLength();   // reference length for theta
  double gamma=setup->getReferenceEnergy();    // get default energy from setup input deck


  map<string,string>::iterator end=arg->end();
  if (arg->find("file")!=end    ){file=arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("time")!=end)    {dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}
  if (arg->find("harmonic")!=end){harm = atoi(arg->at("harmonic").c_str()); arg->erase(arg->find("harmonic"));}
  if (arg->find("replace")!=end) {force_replace = atob(arg->at("replace").c_str()); arg->erase(arg->find("replace"));}
  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &importfield" << endl; this->usage();}
    return false;
  }


  ReadFieldHDF5 import;
  bool check=import.readGlobal(rank, size, file, setup, time, harm, dotime);
  if (!check) { 
    import.close();
    return check; 
  }
  

  // sample rate and time dependent run could have changed when taken by externaldistribution in readGlobal
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  dotime=time->isTime();                                            // check for time simulation

  vector<double> s;
  int nslice=time->getPosition(&s);

  // check for already existing field at specified harmonic
  int idx=-1;
  Field *field=nullptr, *old_field=nullptr;
  for (int i=0; i<fieldin->size();i++){
    if (fieldin->at(i)->harm==harm){
      field=fieldin->at(i);
      idx=i;
    }
  }
  if (idx<0)
  {
    if (rank==0) {cout << "Importing radiation field distribution from file: " << file << " ..." << endl; }
    field=new Field;
    fieldin->push_back(field);
    idx=fieldin->size()-1;
  }
  else
  {
    // we have a field at this harmonic...
    if(!force_replace) {
      if (rank==0) {cout << "*** Error: Cannot import field, because field is already defined" << endl; }
      return false;
    }

    if (rank==0) {
      cout << "Importing radiation field distribution from file: " << file << " (replacing already existing field)..." << endl;
      cout << "   !Replacing of fields is a new function, use with care!" << endl;
    }
    field = new Field;
    // replace already existing field by new one to be filled w/ data (no need to update the index)
    old_field = fieldin->at(idx);
    fieldin->at(idx) = field;
  }

  field->init(time->getNodeNSlice(),import.getNGrid(),import.getDGrid(),lambda,sample*lambda,s[0],harm);

  // read the field, slice by slice
  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    double sloc=s[i];
    import.readSlice(s[i],&field->field[j]);
  }
  import.close();
  
  delete old_field; // ok to delete nullptr -> no action

  return true;
}



