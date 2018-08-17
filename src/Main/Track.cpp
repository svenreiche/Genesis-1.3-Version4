#include "Track.h"
#include "Gencore.h"

Track::Track()
{
  zstop=1e9;
  output_step=1;
  sort_step=0;
  dumpFieldStep=0;
  dumpBeamStep=0;
  bunchharm=1;
}

Track::~Track(){}

void Track::usage(){

  cout << "List of keywords for TRACK" << endl;
  cout << "&track" << endl;
  cout << " double zstop = 1e9" << endl;
  cout << " double s0    = <taken from TIME module>" << endl;
  cout << " double slen  = <taken from TIME module>" << endl;
  cout << " int output_step  = 1" << endl;
  cout << " int field_dump_step  = 0" << endl;
  cout << " int beam_dump_step  = 0" << endl;
  cout << " int sort_step = 0" << endl;
  cout << " int bunchharm = 1" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool Track::init(int inrank, int insize, map<string,string> *arg, Beam *beam, vector<Field *> *field,Setup *setup, Lattice *lat, AlterLattice *alt,Time *time,bool supressOutput)
{
 
  rank=inrank;
  size=insize;
  bunchharm=1; //reset to default for each tracking
  
  bool isTime=time->isTime();
  bool isScan=time->isScan();
  double sample=time->getSampleRate();
  s0=time->getTimeWindowStart();
  slen=time->getTimeWindowLength();

  
  map<string,string>::iterator end=arg->end();

  if (arg->find("zstop")!=end)  {zstop= atof(arg->at("zstop").c_str());  arg->erase(arg->find("zstop"));}
  if (arg->find("s0")!=end)     {s0= atof(arg->at("s0").c_str());  arg->erase(arg->find("s0"));}
  if (arg->find("slen")!=end)   {slen= atof(arg->at("slen").c_str());  arg->erase(arg->find("slen"));}
  if (arg->find("output_step")!=end)   {output_step= atoi(arg->at("output_step").c_str());  arg->erase(arg->find("output_step"));}
  if (arg->find("field_dump_step")!=end)  {dumpFieldStep= atoi(arg->at("field_dump_step").c_str()); arg->erase(arg->find("field_dump_step"));}
  if (arg->find("beam_dump_step")!=end)   {dumpBeamStep = atoi(arg->at("beam_dump_step").c_str());  arg->erase(arg->find("beam_dump_step"));}
  if (arg->find("sort_step")!=end)   {sort_step= atoi(arg->at("sort_step").c_str());  arg->erase(arg->find("sort_step"));}
  if (arg->find("bunchharm")!=end)   {bunchharm= atoi(arg->at("bunchharm").c_str());  arg->erase(arg->find("bunchharm"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &track" << endl; this->usage();}
    return false;
  }


  if (output_step < 1) { output_step=1; }

  string file;
  setup->getRootName(&file);
  file.append(".out.h5");


 
  Undulator *und = new Undulator;

  lat->generateLattice(setup->getStepLength(),setup->getReferenceLength(),setup->getReferenceEnergy(),alt, und);  
  und->updateOutput(zstop,output_step);
  und->updateMarker(dumpFieldStep,dumpBeamStep,sort_step,zstop);

  beam->setBunchingHarmonicOutput(bunchharm);
  // call to gencore to do the actual tracking.  

  Gencore core;
  core.run(file.c_str(),beam,field,und,isTime,isScan,supressOutput);


  delete und;
   
  if  (rank==0) { cout << "End of Track" << endl;}
 
  return true;

}
