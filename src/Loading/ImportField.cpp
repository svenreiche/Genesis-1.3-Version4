#include "ImportField.h"
#include "readFieldHDF5.h"


ImportField::ImportField()
{
  offset=0;
  harm=1;
  dotime=true;
}

ImportField::~ImportField(){}

void ImportField::usage(){

  cout << "List of keywords for IMPORTFIELD" << endl;
  cout << "&importfield" << endl;
  cout << " string file = <empty>" << endl;
  cout << " int harmonic = 1" << endl;
  cout << " bool time = true" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool ImportField::init(int rank, int size, map<string,string> *arg, vector<Field *> *fieldin, Setup *setup, Time *time)
{

 
    
  double lambda=setup->getReferenceLength();   // reference length for theta
  double gamma=setup->getReferenceEnergy();           // get default energy from setup input deck


  map<string,string>::iterator end=arg->end();

  if (arg->find("file")!=end    ){file=arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("time")!=end)    {dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}
  if (arg->find("harmonic")!=end){harm = atoi(arg->at("harmonic").c_str()); arg->erase(arg->find("harmonic"));}


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


  int idx=-1;
  Field *field;
  for (int i=0; i<fieldin->size();i++){
    if (fieldin->at(i)->harm==harm){
      field=fieldin->at(i);
      idx=i;
    }
  }
  if (idx<0){
    field=new Field;
    if (rank==0) {cout << "Importing radiation field distribution from file: " << file << " ..." << endl; }
  } else {
    if (rank==0) {cout << "*** Error: Cannot import field, because field is already defined" << endl; }
    return false;
  }



  field->init(time->getNodeNSlice(),import.getNGrid(),import.getDGrid(),lambda,sample*lambda,s[0],harm);
  
  if (idx<0){
    fieldin->push_back(field);
    idx=fieldin->size()-1;
  }


 

  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    double sloc=s[i];
    import.readSlice(s[i],&field->field[j]);
  }
  import.close();
  

  return true;

}



