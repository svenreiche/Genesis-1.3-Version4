#include "ImportBeam.h"
#include "readBeamHDF5.h"


ImportBeam::ImportBeam()
{
  offset=0;
  dotime=true;
}

ImportBeam::~ImportBeam(){}

void ImportBeam::usage(){

  cout << "List of keywords for IMPORTBEAM" << endl;
  cout << "&importbeam" << endl;
  cout << " string file = <empty>" << endl;
  cout << " bool time = true" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool ImportBeam::init(int rank, int size, map<string,string> *arg, Beam *beam, Setup *setup, Time *time)
{


  if (beam->beam.size()>0){
    if (rank==0) {cout << "*** Error: Cannot import beam, because beam is already defined" << endl; }
    return false;
  }
    
  double lambda=setup->getReferenceLength();   // reference length for theta
  bool one4one=setup->getOne4One();            // check for one4one simulations
  double gamma=setup->getReferenceEnergy();           // get default energy from setup input deck
  int nbins=setup->getNbins();


  map<string,string>::iterator end=arg->end();

  if (arg->find("file")!=end    ){file=arg->at("file"); arg->erase(arg->find("file"));}
  if (arg->find("time")!=end)    {dotime = atob(arg->at("time").c_str()); arg->erase(arg->find("time"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &importbeam" << endl; this->usage();}
    return false;
  }

  

  if (rank==0){
    cout << "Importing particle distribution from file: " << file << " ..." << endl; 
  }

  ReadBeamHDF5 import;

  bool check=import.readGlobal(rank, size, file, setup, time, dotime);
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



  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    double sloc=s[i];
    import.readSlice(s[i],&beam->beam[j],&beam->current[j],one4one);
  }
  import.close();

  return true;

}


