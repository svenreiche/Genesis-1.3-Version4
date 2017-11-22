#include "ImportField.h"
#include "readFieldHDF5.h"
#include "MPEProfiling.h"


ImportField::ImportField()
{
  offset=0;
  dotime=true;
}

ImportField::~ImportField(){}

void ImportField::usage(){

  cout << "List of keywords for IMPORTFIELD" << endl;
  cout << "&importfield" << endl;
  cout << " string file = <empty>" << endl;
  cout << " double offset = 0" << endl;
  cout << " bool time = true" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool ImportField::init(int rank, int size, map<string,string> *arg, vector<Field *> *beam, Setup *setup, Time *time)
{

  
  if (beam->size()>0){
    if (rank==0) {cout << "*** Error: Cannot import field, because field is already defined" << endl; }
    return false;
  }
    
  double lambda=setup->getReferenceLength();   // reference length for theta
  double gamma=setup->getReferenceEnergy();           // get default energy from setup input deck


  map<string,string>::iterator end=arg->end();

  if (arg->find("file")!=end    ){file=arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("offset")!=end  ){offset=atof(arg->at("offset").c_str());  arg->erase(arg->find("offset"));}
  if (arg->find("time")!=end)    {dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &importfield" << endl; this->usage();}
    return false;
  }

  

  if (rank==0){
    cout << "Analysing field distribution..." << endl; 
  }

  // not yet working!!!!!!!!!!!!!!!

  /*

  ReadBeamHDF5 import;

  bool check=import.readGlobal(rank, size, file, setup, time, offset,dotime);
  if (!check) { 
    import.close();
    return check; 
  }

  // sample rate and time dependent run could have changed when taken by externaldistribution in readGlobal
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  dotime=time->isTime();                                            // check for time simulation

  vector<double> s;
  int nslice=time->getPosition(&s);
  beam->init(time->getNodeNSlice(),nbins,lambda,sample*lambda,s[0],one4one);


  mpe.logLoading(false,"Importing Particle Distribution");

  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    double sloc=s[i];
    import.readSlice(s[i],&beam->beam[j],&beam->current[j]);
  }
  import.close();
  
  mpe.logLoading(true,"Importing Particle Distribution");
  */ 
  mpe.logEvent("End: ImportField::init");

  return true;

}



