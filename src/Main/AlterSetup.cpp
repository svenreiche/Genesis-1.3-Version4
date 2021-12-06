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

  bool my_beam_write_filter=false;
  bool update_beam_write_filter=false;

  if (arg->find("rootname")!=end){rootname = arg->at("rootname"); arg->erase(arg->find("rootname"));}
  //  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice"); arg->erase(arg->find("lattice"));}
  if (arg->find("beamline")!=end){beamline = arg->at("beamline"); arg->erase(arg->find("beamline"));}
  if (arg->find("delz")!=end)    {delz     = atof(arg->at("delz").c_str());  arg->erase(arg->find("delz"));}
  if (arg->find("subharmonic")!=end){subharmonic  = atoi(arg->at("subharmonic").c_str());  arg->erase(arg->find("subharmonic"));}
  if (arg->find("harmonic")!=end){harmonic  = atoi(arg->at("harmonic").c_str());  arg->erase(arg->find("harmonic"));}
  if (arg->find("resample")!=end){resample  = atob(arg->at("resample"));  arg->erase(arg->find("resample"));}
  if (arg->find("disable")!=end){disable  = atob(arg->at("disable"));  arg->erase(arg->find("disable"));}


  /* electron beam slice downselector parameters */
  // If about to do harmonic/subharmonic conversion:
  // (1) load defaults
  //     --> do it before parsing the corresponding user input so that updated values before effective
  // (2) disable write selector (unless requested otherwise)
  //     --> done one we have parsed the corresponding user input
  if ((harmonic>1) || (subharmonic>1)) {
    if(setup->BWF_get_enabled()) {
      if(rank==0) {
        cout << "(sub-)harmonic conversion to be performed: loading default parameters for beam write slice selector" << endl;
      }
      setup->BWF_load_defaults(); // load defaults (but do not change 'enabled' flag)
    }
  }

  // same code as in Setup.cpp
  if (arg->find("beam_write_slices_from")!=end) {
    int t = atoi(arg->at("beam_write_slices_from").c_str());
    arg->erase(arg->find("beam_write_slices_from"));
    setup->BWF_set_from(t);
    my_beam_write_filter=true;     // user can override this if needed
    update_beam_write_filter=true; // this is updated only after the harmonic conversion steps
  }
  if (arg->find("beam_write_slices_to")!=end) {
    int t = atoi(arg->at("beam_write_slices_to").c_str());
    arg->erase(arg->find("beam_write_slices_to"));
    setup->BWF_set_to(t);
    my_beam_write_filter=true;     // user can override this if needed
    update_beam_write_filter=true; // this is updated only after the harmonic conversion steps
  }
  if (arg->find("beam_write_slices_inc")!=end) {
    int t = atoi(arg->at("beam_write_slices_inc").c_str());
    arg->erase(arg->find("beam_write_slices_inc"));
    setup->BWF_set_inc(t);
    my_beam_write_filter=true;     // user can override this if needed
    update_beam_write_filter=true; // this is updated only after the harmonic conversion steps
  }
  if (arg->find("beam_write_slices_filter")!=end) {
    my_beam_write_filter = atob(arg->at("beam_write_slices_filter"));
    update_beam_write_filter=true; // this is updated only after the harmonic conversion steps
    arg->erase(arg->find("beam_write_slices_filter"));
  }

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &alter_setup" << endl; this->usage();}
    return false;
  }

  // if (sub-)harmonic conversion was requested and slice selector is enabled: disable it, unless
  // one the parameters was touched by the user input
  if ((harmonic>1) || (subharmonic>1)) {
    if(!update_beam_write_filter) {
      if(rank==0) {
        cout << "(sub-)harmonic conversion: disabling beam write slice selector" << endl;
      }
      setup->BWF_set_enabled(false);
    }
  }

  setup->setStepLength(delz);
  if (!time->isTime()){ resample=false; }  // no resampling in non-timedependent runs allowed

  // step one: Select new lattice if selected
  if (beamline!="") {
    bool status = lat->parse(lattice,beamline,rank);
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

    // this process changed the slice numbering: slice downselector (for ebeam dumping was already disabled above)
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
	    delete field->at(i);
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

    // this process changed the slice numbering: slice downselector (for ebeam dumping was already disabled above)
  }

  // enable beam write filter again if requested (logic is identical to that in &setup)
  if(update_beam_write_filter)
    setup->BWF_set_enabled(my_beam_write_filter);

  return true;
}


