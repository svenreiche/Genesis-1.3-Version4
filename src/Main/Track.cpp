#include "Track.h"
#include "Gencore.h"

#include "BeamDiag_Demo.h"

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
  cout << " bool field_dump_at_undexit = false" << endl;
  cout << " int beam_dump_step  = 0" << endl;
  cout << " int sort_step = 0" << endl;
  cout << " int bunchharm = 1" << endl;
  cout << "&end" << endl << endl;
  /* currently undocumented debugging option: dbg_report_lattice */

  return;
}

/* from StringProcessing.cpp */
bool Track::atob(string in)
{
	bool ret=false;
	if ((in.compare("1")==0)||(in.compare("true")==0)||(in.compare("t")==0)) { ret=true; }
	return ret;
}

bool Track::init(int inrank, int insize, map<string,string> *arg, Beam *beam, vector<Field *> *field,Setup *setup, Lattice *lat, AlterLattice *alt,Time *time)
{
  rank=inrank;
  size=insize;
  bunchharm=1; //reset to default for each tracking

  bool dbg_report_lattice=false;  
  bool dbg_report_moddiag=false;
  bool dumpFieldUE=false;
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
  if (arg->find("field_dump_at_undexit")!=end) {dumpFieldUE = atob(arg->at("field_dump_at_undexit")); arg->erase(arg->find("field_dump_at_undexit"));}
  if (arg->find("beam_dump_step")!=end)   {dumpBeamStep = atoi(arg->at("beam_dump_step").c_str());  arg->erase(arg->find("beam_dump_step"));}
  if (arg->find("sort_step")!=end)   {sort_step= atoi(arg->at("sort_step").c_str());  arg->erase(arg->find("sort_step"));}
  if (arg->find("bunchharm")!=end)   {bunchharm= atoi(arg->at("bunchharm").c_str());  arg->erase(arg->find("bunchharm"));}
  if (arg->find("dbg_report_lattice")!=end) {dbg_report_lattice = atob(arg->at("dbg_report_lattice")); arg->erase(arg->find("dbg_report_lattice"));}
  if (arg->find("dbg_report_moddiag")!=end) {dbg_report_moddiag = atob(arg->at("dbg_report_moddiag")); arg->erase(arg->find("dbg_report_moddiag"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &track" << endl; this->usage();}
    return false;
  }


  if (output_step < 1) { output_step=1; }

  string file;
  setup->getRootName(&file);
  file.append(".out.h5");


 
  Undulator *und = new Undulator;

  lat->generateLattice(setup, alt, und); /* !changes to 'lat' after this function call may have no effect, as this function also stores the generated lattice into the provided Undulator instance! */
  if(dumpFieldUE)
    und->markUndExits();

  und->updateOutput(zstop,output_step);
  und->updateMarker(dumpFieldStep,dumpBeamStep,sort_step,zstop);
  beam->setBunchingHarmonicOutput(bunchharm);

  // controling the output

  bool ssrun=true;    
  if ((size > 1) or (beam->beam.size()>1)){
    ssrun=false;              // no steady-state run if there is more than one core or more than one slice
  }
  if (ssrun){    // disable global output when there is only one slice calculated 
    beam->set_global_stat(false);
    for (int i=0; i<field->size();i++){
      field->at(i)->set_global_stat(false);
    }
  } else {
    beam->set_global_stat(setup->getBeamGlobalStat());
    for (int i=0; i<field->size();i++){
      field->at(i)->set_global_stat(setup->getFieldGlobalStat());
    }
  } 
  for (int i=0; i<field->size();i++){
    field->at(i)->setOutput(
       setup->outputFFT(),
       setup->outputSpatial(),
       setup->outputIntensity(),
       setup->outputFieldDump());
  }
  beam->setOutput(setup->outputCurrent(),setup->outputEnergy(),setup->outputSpatial(),setup->outputAux());

  if(dbg_report_lattice && (rank==0)) {
    stringstream ss;
    int rc = setup->getCount();

    /* generate report filename (same logic as in Setup::getRootName) */
    ss << "latreport_u"; // "_u" indicates lattice report from Undulator class (also Lattice class can write lattice report)
    if(rc>0) {
      ss << ".Run" << (rc+1);
    }
    ss << ".txt";
    und->reportLattice(ss.str());
  }


  // Setup beam diagnostics modules (their 'do_diag_ member function is 
  // called at every integration step).
  // !Note: Do not free these here, the call to 'clear_beam_diag' releases!
  // !the allocated memory.                                               !
// #define DO_BEAMDIAG_DEMO
#ifdef DO_BEAMDIAG_DEMO
  BeamDiag_Demo *bd_demo = new BeamDiag_Demo();
  // configure the demo diag module (every diag module provides
  // its specific configuration functions, if needed)
  bd_demo->config(12345);
  bd_demo->set_verbose(true);
  // Register demo module
  beam->register_beam_diag(bd_demo);
#endif

  if((0==rank) && dbg_report_moddiag)
    beam->beam_diag_list_registered();

  // call to gencore to do the actual tracking.  
  Gencore core;
  core.run(file.c_str(),beam,field,und,isTime,isScan);

  // Clear beam diagnostics (this instance of 'Beam' class could be re-used for the next &track command)
  // Note: This also destroys the beam diagnostics objects
  beam->clear_beam_diag();

  delete und;
   
  if  (rank==0) { cout << "End of Track" << endl;}
 
  return true;

}
