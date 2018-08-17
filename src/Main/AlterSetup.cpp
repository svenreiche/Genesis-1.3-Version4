#include "AlterSetup.h"


AlterSetup::AlterSetup()
{
  rootname="";
  lattice="";
  beamline="";
  harmonic=1;
  subharmonic=1;
  resample=false;
  disable=false;

}

AlterSetup::~AlterSetup(){}

void AlterSetup::usage(){

  cout << "List of keywords for ALTER_SETUP" << endl;
  cout << "&alter_setup" << endl;
  cout << " string rootname = <taken from SETUP + Increment>" << endl;
  //  cout << " string lattice = <taken from SETUP>" << endl;
  cout << " string beamline = <empty>" << endl;
  cout << " double delz = <taken from SETUP>" << endl;
  cout << " int harmonic = 1" << endl;
  cout << " int subharmonic = 1" << endl;
  cout << " bool resample = false" << endl;
  cout << " bool disable = false" << endl;
  cout << "&end" << endl << endl;
  return;
}

bool AlterSetup::init(int inrank, map<string,string> *arg, Setup *setup, Lattice *lat, Time *time, Beam *beam, vector<Field *> *field)
{

  beamline="";  // clear beamline name for multiple calls of altersetup
  lattice=setup->getLattice();
  delz=setup->getStepLength();
  rank=inrank;
  
  map<string,string>::iterator end=arg->end();

  if (arg->find("rootname")!=end){rootname = arg->at("rootname"); arg->erase(arg->find("rootname"));}
  //  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice"); arg->erase(arg->find("lattice"));}
  if (arg->find("beamline")!=end){beamline = arg->at("beamline"); arg->erase(arg->find("beamline"));}
  if (arg->find("delz")!=end)    {delz     = atof(arg->at("delz").c_str());  arg->erase(arg->find("delz"));}
  if (arg->find("subharmonic")!=end){subharmonic  = atoi(arg->at("subharmonic").c_str());  arg->erase(arg->find("subharmonic"));}
  if (arg->find("harmonic")!=end){harmonic  = atoi(arg->at("harmonic").c_str());  arg->erase(arg->find("harmonic"));}
  if (arg->find("resample")!=end){resample  = atob(arg->at("resample"));  arg->erase(arg->find("resample"));}
  if (arg->find("disable")!=end){disable  = atob(arg->at("disable"));  arg->erase(arg->find("disable"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &alter_setup" << endl; this->usage();}
    return false;
  }

  setup->setStepLength(delz);
  if (!time->isTime()){ resample=false; }  // no resampling in non-timedependent runs allowed

  // step one: Select new lattice if selected
  if (beamline!="") {
    bool status = lat->parse(lattice,beamline,rank, false);
     if (status==false) { return status ; }
  }

  // step two: Add increment for rootname or assign new root name
  if (rootname!=""){
    setup->setRootName(&rootname); 
  } else {
    setup->incrementCount();
  }

  //step three: do subharmonic conversion
  // step 3.1 - time window adjustment
  if (subharmonic>1){
    if (rank==0){ cout << "Converting to subharmonic: " << subharmonic << "..." << endl;}
    double lam0=setup->getReferenceLength()*static_cast<double>(subharmonic);
    setup->setReferenceLength(lam0);
    if (!resample) {
      double samp=time->getSampleRate()/static_cast<double>(subharmonic);
      time->setSampleRate(samp);
    }
    time->finishInit(setup);
    // step 3.2 - beam
    if (!beam->subharmonicConversion(subharmonic,resample)){ 
      if (rank==0) {cout << "*** Error: Cannot convert beam distribution to lower harmonic" << endl;}
      return false;
    }
    // step 3.3 - field
    for (int i=0; i<field->size();i++){
      if (field->at(i)->getHarm()!=1){
	if (rank==0) {cout << "Deleting higher radiation harmonic: " << field->at(i)->getHarm() << " in subharmonic conversion" << endl;}
	field->erase(field->begin()+i);
	i--; // reseting the count
      }else{
	if (rank==0) {cout << "Converting radiation fundamental to harmonic: " << subharmonic << " and keep it in memory" << endl;}
	 if (!(field->at(i)->subharmonicConversion(subharmonic,resample))){
	    if (rank==0) {cout << "*** Error: Cannot convert field distribution to higher harmonic" << endl;}
	    return false;	    
	 }
      }
    }
  }

  // step four: do harmonic conversion
  // step 4.1 - time window
  if (harmonic>1){
    if (rank==0){ cout << "Converting to harmonic: " << harmonic << "..." << endl;}
    double lam0=setup->getReferenceLength()/static_cast<double>(harmonic);
    setup->setReferenceLength(lam0);
    if (!resample) {
      double sam=time->getSampleRate()*static_cast<double>(harmonic);
      time->setSampleRate(sam);
    }
    time->finishInit(setup);

    // step 4.2 - beam
    if (!beam->harmonicConversion(harmonic,resample)){ 
      if (rank==0) {cout << "*** Error: Cannot convert beam distribution to higher harmonic" << endl;}
      return false;
    }
    
    // step 4.3 - field
    for (int i=0; i<field->size();i++){
      if (field->at(i)->getHarm()!=harmonic){
        if (disable) {
	    if (rank==0) {cout << "Disabling non-matching radiation harmonic: " << field->at(i)->getHarm() << " in harmonic conversion" << endl;}
	    field->at(i)->disable(1./static_cast<double>(harmonic));

	} else {
	    if (rank==0) {cout << "Deleting non-matching radiation harmonic: " << field->at(i)->getHarm() << " in harmonic conversion" << endl;}
	    field->erase(field->begin()+i);
	    i--; // reseting the count
        }
      }else{
	if (rank==0) {cout << "Converting radiation harmonic: " << harmonic << " to fundamental and keep it in memory" << endl;}
	 if (!(field->at(i)->harmonicConversion(harmonic,resample))){
	    if (rank==0) {cout << "*** Error: Cannot convert field distribution to higher harmonic" << endl;}
	    return false;	    
	 }
      }
  
    }
  }

  return true;


}


