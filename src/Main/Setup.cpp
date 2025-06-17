#include <sstream>
#include "Setup.h"
#ifdef USE_DPI
  #include "DiagnosticHookS.h"
#endif

Setup::Setup()
{
  rootname="";
  outputdir="";
  lattice="";
  beamline="";
  partfile="";
  fieldfile="";
  one4one=false;
  shotnoise=true;
  nbins=4;
  npart=8192;
  gamma0=5800/0.511;
  lambda0=1e-10;
  delz=0.015; 
  seed=123456789;


  beam_global_stat=false;
  field_global_stat=false;

  exclude_spatial_output=false;
  exclude_fft_output=false;
  exclude_intensity_output=false;
  exclude_energy_output=false;
  exclude_aux_output=false;
  exclude_current_output=true;
  exclude_twiss_output=true;
  exclude_field_dump=false;
  do_write_outfile=true;

  // filtering of beam slices during dump process (information is forwarded into active instance of Beam class when actually needed there)
  BWF_set_enabled(false);
  BWF_load_defaults();

  // count of runs in conjunction of calls of altersetup
  runcount = 0;
  write_meta_file=false;
  sema_file_enabled_start=false;
  sema_file_enabled_done=false;
}

Setup::~Setup()= default;

void Setup::usage(){

  cout << "List of keywords for SETUP" << endl;
  cout << "&setup" << endl;
  cout << " string rootname = <taken from command line>" << endl;
  cout << " string outputdir = <empty>" << endl;
  cout << " string lattice  = <taken from command line>" << endl;
  cout << " string beamline = <empty>" << endl;
  cout << " double gamma0 = 5800/0.511" << endl;
  cout << " double lambda0 = 1e-10" << endl;
  cout << " double delz = 0.015" << endl;
  cout << " int seed = 123456789" << endl;
  cout << " int npart = 8192" << endl;
  cout << " int nbins = 4" << endl;
  cout << " bool one4one = false" << endl;
  cout << " bool shotnoise = true" << endl;
  cout << " bool beam_global_stat = false" << endl;
  cout << " bool field_global_stat = false" << endl;
  cout << " bool exclude_spatial_output = false" << endl;
  cout << " bool exclude_fft_output = false" << endl;
  cout << " bool exclude_intensity_output = false" << endl;
  cout << " bool exclude_energy_output = false" << endl;
  cout << " bool exclude_aux_output = false" << endl;
  cout << " bool exclude_current_output = true" << endl;
  cout << " bool exclude_twiss_output = true" << endl;
  cout << " bool exclude_field_dump = false" << endl;
  cout << " bool write_meta_file = false" << endl;
  cout << " bool write_semaphore_file = false" << endl;
  cout << " bool write_semaphore_file_done = false" << endl;
  cout << " bool write_semaphore_file_started = false" << endl;
  cout << " string semaphore_file_name = <derived from 'rootname'>" << endl;
  cout << "&end" << endl << endl;
}

bool Setup::init(int inrank, map<string,string> *arg, Lattice *lat, SeriesManager *sm, FilterDiagnostics &filter){

  rank=inrank;
  auto end=arg->end();

  if (arg->find("rootname")!=end){rootname = arg->at("rootname"); arg->erase(arg->find("rootname"));}
  if (arg->find("outputdir")!=end){outputdir = arg->at("outputdir"); arg->erase(arg->find("outputdir"));}
  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice");  arg->erase(arg->find("lattice"));}
  if (arg->find("beamline")!=end){beamline = arg->at("beamline"); arg->erase(arg->find("beamline"));}
  if (arg->find("lattice")!=end) {lattice  = arg->at("lattice");  arg->erase(arg->find("lattice"));}
  if (arg->find("gamma0")!=end)  {gamma0   = atof(arg->at("gamma0").c_str());  arg->erase(arg->find("gamma0"));}
  if (arg->find("lambda0")!=end) {lambda0  = atof(arg->at("lambda0").c_str()); arg->erase(arg->find("lambda0"));}
  if (arg->find("delz")!=end)    {delz     = atof(arg->at("delz").c_str());  arg->erase(arg->find("delz"));}
  if (arg->find("seed")!=end)    {seed     = atoi(arg->at("seed").c_str());  arg->erase(arg->find("seed"));}
  if (arg->find("one4one")!=end) {one4one  = atob(arg->at("one4one"));  arg->erase(arg->find("one4one"));}
  if (arg->find("npart")!=end)    {npart  = atoi(arg->at("npart").c_str());  arg->erase(arg->find("npart"));}
  if (arg->find("nbins")!=end)    {nbins  = atoi(arg->at("nbins").c_str());  arg->erase(arg->find("nbins"));}
  if (arg->find("shotnoise")!=end){shotnoise  = atob(arg->at("shotnoise"));  arg->erase(arg->find("shotnoise"));}
  if (arg->find("beam_global_stat")!=end)  {beam_global_stat  = atob(arg->at("beam_global_stat"));   arg->erase(arg->find("beam_global_stat"));}
  if (arg->find("field_global_stat")!=end) {field_global_stat = atob(arg->at("field_global_stat"));  arg->erase(arg->find("field_global_stat"));}

  if (arg->find("exclude_spatial_output")!=end)   {exclude_spatial_output  = atob(arg->at("exclude_spatial_output"));   arg->erase(arg->find("exclude_spatial_output"));}
  if (arg->find("exclude_fft_output")!=end)       {exclude_fft_output      = atob(arg->at("exclude_fft_output"));       arg->erase(arg->find("exclude_fft_output"));}
  if (arg->find("exclude_intensity_output")!=end) {exclude_intensity_output= atob(arg->at("exclude_intensity_output")); arg->erase(arg->find("exclude_intensity_output"));}
  if (arg->find("exclude_energy_output")!=end)    {exclude_energy_output   = atob(arg->at("exclude_energy_output"));    arg->erase(arg->find("exclude_energy_output"));}
  if (arg->find("exclude_aux_output")!=end)       {exclude_aux_output      = atob(arg->at("exclude_aux_output"));       arg->erase(arg->find("exclude_aux_output"));}
  if (arg->find("exclude_current_output")!=end)   {exclude_current_output  = atob(arg->at("exclude_current_output"));   arg->erase(arg->find("exclude_current_output"));}
  if (arg->find("exclude_twiss_output")!=end)   {exclude_twiss_output  = atob(arg->at("exclude_twiss_output"));   arg->erase(arg->find("exclude_twiss_output"));}
  if (arg->find("exclude_field_dump")!=end)   {exclude_field_dump  = atob(arg->at("exclude_field_dump"));   arg->erase(arg->find("exclude_field_dump"));}

  if (arg->find("write_meta_file")!=end)   {write_meta_file = atob(arg->at("write_meta_file"));   arg->erase(arg->find("write_meta_file"));}
  if (arg->find("write_semaphore_file")!=end)   {sema_file_enabled_done  = atob(arg->at("write_semaphore_file"));   arg->erase(arg->find("write_semaphore_file"));}
  /* alias for write_semaphore_file */
  if (arg->find("write_semaphore_file_done")!=end)   {sema_file_enabled_done  = atob(arg->at("write_semaphore_file_done"));   arg->erase(arg->find("write_semaphore_file_done"));}
  if (arg->find("write_semaphore_file_started")!=end)   {sema_file_enabled_start  = atob(arg->at("write_semaphore_file_started"));   arg->erase(arg->find("write_semaphore_file_started"));}
  if (arg->find("semaphore_file_name")!=end) {
    // Providing a file name for the semaphore file always switches on writing the "done" semaphore file, overriding 'write_semaphore_file' flag.
    // This allows to switch on semaphore functionality just by specifying corresponding command line argument -- no modification of G4 input file needed.
    sema_file_enabled_done = true;
    setSemaFN(arg->at("semaphore_file_name")); arg->erase(arg->find("semaphore_file_name"));
  }
  /* Note: if requested, semaphore file of type "started" is generated in GenMain after successfully returning from Setup::init */

  // same code also in AlterSetup.cpp
  if (arg->find("beam_write_slices_from")!=end) {
    beam_write_slices_from = atoi(arg->at("beam_write_slices_from").c_str());
    arg->erase(arg->find("beam_write_slices_from"));
    beam_write_filter=true; // user can override this if needed
  }
  if (arg->find("beam_write_slices_to")!=end) {
    beam_write_slices_to = atoi(arg->at("beam_write_slices_to").c_str());
    arg->erase(arg->find("beam_write_slices_to"));
    beam_write_filter=true; // user can override this if needed
  }
  if (arg->find("beam_write_slices_inc")!=end) {
    BWF_set_inc(atoi(arg->at("beam_write_slices_inc").c_str())); // setter function does checks if valid value
    arg->erase(arg->find("beam_write_slices_inc"));
    beam_write_filter=true; // user can override this if needed
  }
  if (arg->find("beam_write_slices_filter")!=end) {
    beam_write_filter = atob(arg->at("beam_write_slices_filter"));
    arg->erase(arg->find("beam_write_slices_filter"));
  }

  if (!arg->empty()){
    if (rank==0){ cout << "*** Error: Unknown elements in &setup" << endl; this->usage();}
    return false;
  }

  // sort the filter flags
  filter.beam.global = beam_global_stat;
  filter.beam.spatial = !exclude_spatial_output;
  filter.beam.energy = !exclude_energy_output;
  filter.beam.current = !exclude_current_output;
  filter.beam.twiss = !exclude_twiss_output;
  filter.beam.auxiliar = !exclude_aux_output;
  filter.field.global = field_global_stat;
  filter.field.spatial = !exclude_spatial_output;
  filter.field.fft = !exclude_fft_output;
  filter.field.intensity = !exclude_intensity_output;

  if (one4one)
  {
    nbins = 1;
  }

  return lat->parse(lattice,beamline,rank,sm);
}



bool Setup::getRootName(string *filename)
{
  if (rootname.empty()){
    return false;
  }
  *filename=rootname;
  if (runcount>0) {
    stringstream ss;
    ss << ".Run" << (runcount+1) ;
    *filename+=ss.str();
  }
  return true; 
}

bool Setup::RootName_to_FileName(string *fn_out, string *filename)
{
    // replace placeholder symbol @ with rootname
    std::string placeholder("@");
    int pos;
    while ((pos = filename->find(placeholder)) != std::string::npos)
        filename->replace(pos, placeholder.length(),rootname);

  // If 'outputdir' parameter not specified
  // ==> do not change string

  if (outputdir.empty()) {
    *fn_out = *filename;
    return true;
  }

  string t;
  t = outputdir + "/" + *filename;
  *fn_out = t;

  return true;
}

void Setup::BWF_load_defaults()
{
  // note: only loading defaults, not disabling the filter (if desired, this has to be done in addition)
  beam_write_slices_from=-1;
  beam_write_slices_to=-1;
  beam_write_slices_inc=1;
}




/* returns true if filename for semaphore file was generated */
bool Setup::getSemaFN(string *fnout) {
	// user-defined filename overrides the logic to derive from rootname (also a user-provided 'outputdir' parameter does not have an effect here)
	if(! sema_file_name.empty()) {
		*fnout = sema_file_name;
		return(true);
	}

	if(rootname.empty()) {
		return(false);
	}
	string t;
	t = rootname;
	t += ".sema";
	RootName_to_FileName(fnout, &t);
	return(true);
}
void Setup::setSemaFN(string fn_in) {
	sema_file_name = fn_in;
}
