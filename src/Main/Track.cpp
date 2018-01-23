#include "Track.h"
#include "Gencore.h"

Track::Track()
{
  zstop=1e9;
  output_step=1;
  dumpFieldStep=0;
  dumpBeamStep=0;
  sortStep=0;
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
  cout << "&end" << endl << endl;
  return;
}

bool Track::init(int inrank, int insize, map<string,string> *arg, Beam *beam, vector<Field *> *field,Setup *setup, Lattice *lat, AlterLattice *alt,Time *time)
{
 
  rank=inrank;
  size=insize;
  
  bool isTime=time->isTime();
  bool isScan=time->isScan();
  double sample=time->getSampleRate();
  s0=time->getTimeWindowStart();
  slen=time->getTimeWindowLength();

  
  map<string,string>::iterator end=arg->end();

  if (arg->find("zstop")!=end)  {zstop= atof(arg->at("zstop").c_str());  arg->erase(arg->find("zstop"));}
  if (arg->find("s0")!=end)     {s0= atof(arg->at("s0").c_str());  arg->erase(arg->find("s0"));}
  if (arg->find("slen")!=end)   {slen= atof(arg->at("slen").c_str());  arg->erase(arg->find("slen"));}
  if (arg->find("output_step")!=end)   {output_step= atof(arg->at("output_step").c_str());  arg->erase(arg->find("output_step"));}
  if (arg->find("field_dump_step")!=end)  {dumpFieldStep= atof(arg->at("field_dump_step").c_str()); arg->erase(arg->find("field_dump_step"));}
  if (arg->find("beam_dump_step")!=end)   {dumpBeamStep = atof(arg->at("beam_dump_step").c_str());  arg->erase(arg->find("beam_dump_step"));}
  if (arg->find("sort_step")!=end)   {output_step= atof(arg->at("sort_step").c_str());  arg->erase(arg->find("sort_step"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &track" << endl; this->usage();}
    return false;
  }


  if (output_step < 1) { output_step=1; }

  string file;
  setup->getRootName(&file);
  file.append(".out.h5");




  if (rank==0){

    hid_t fid=H5Fcreate(file.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT); 



    lat->writeLattice(fid,setup->getStepLength(),setup->getReferenceLength(),setup->getReferenceEnergy(),alt);
    setup->writeGlobal(fid,zstop,output_step,dumpFieldStep,dumpBeamStep,sortStep,s0,slen,sample,isTime,isScan);

    H5Fclose(fid);  

  }

  MPI::COMM_WORLD.Barrier(); // synchronize all nodes till root has finish writing the output file

  // call to gencore to do the actual tracking.  

  Gencore core;
  core.run(file.c_str(),beam,field);

  if  (rank==0) { cout << "End of Track" << endl;}
 
  return true;

}
