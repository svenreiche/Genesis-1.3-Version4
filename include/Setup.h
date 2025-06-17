#ifndef __GENESIS_SETUP__
#define __GENESIS_SETUP__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <cctype>
#include <mpi.h>

#include "StringProcessing.h"
#include "Lattice.h"
#include "Diagnostic.h"
#ifdef USE_DPI
  class DiagFieldPluginCfg;
  class DiagBeamPluginCfg;
#endif

using namespace std;

class Lattice;
class SeriesManager;

class Setup: public StringProcessing{
 public:
   Setup();
   virtual ~Setup();
   bool init(int, map<string,string> *, Lattice *, SeriesManager *sm, FilterDiagnostics &filter);
   //   bool writeGlobal(hid_t, double, int, int,int,int,double, double, double, bool, bool);
   double getReferenceLength();
   void   setReferenceLength(double);
   double getReferenceEnergy();
   double getStepLength();
   void   setStepLength(double);
   bool   getOne4One();
   bool   getShotNoise();
   bool   getBeamGlobalStat();
   bool   getFieldGlobalStat();
   bool   outputFFT();
   bool   outputSpatial();
   bool   outputIntensity();
   bool   outputCurrent();
   bool   outputEnergy();
   bool   outputAux();
   bool   outputFieldDump();
   bool   getSemaEnStart();
   bool   getSemaEnDone();
   void   setSemaFN(string);
   bool   getSemaFN(string *);
   bool   get_write_meta_file();
   void   set_do_write_outfile(bool);
   bool   get_do_write_outfile();

   int    getNpart();
   int    getNbins();
   int    getSeed();
   //   bool   getInputFileNameBeam(string *);
   //   bool   getInputFileNameField(string *);
   bool   getRootName(string *);
   void   setRootName(string *);
   void   getOutputdir(string *);
   bool   RootName_to_FileName(string *, string *);
   int getCount();
   void incrementCount();
   string getLattice();

   void   BWF_load_defaults();
   bool   BWF_get_enabled(); // BWF=beam write filter
   void   BWF_set_enabled(bool);
   int    BWF_get_from();
   void   BWF_set_from(int);
   int    BWF_get_to();
   void   BWF_set_to(int);
   int    BWF_get_inc();
   void   BWF_set_inc(int);

#ifdef USE_DPI
   // FIXME: make this 'private' and implement functions for access
   std::vector<DiagFieldPluginCfg> diagpluginfield_;
   std::vector<DiagBeamPluginCfg> diagpluginbeam_;
#endif

 private:
   static void usage();
   string rootname,outputdir,lattice,beamline,partfile,fieldfile;
   double gamma0,lambda0,delz;
   bool one4one,shotnoise;
   bool beam_global_stat, field_global_stat;
   bool exclude_spatial_output, exclude_fft_output, exclude_intensity_output, exclude_energy_output, exclude_aux_output, exclude_current_output, exclude_twiss_output, exclude_field_dump;
   bool do_write_outfile;

   bool beam_write_filter;
   int beam_write_slices_from, beam_write_slices_to, beam_write_slices_inc;

   int seed, rank,npart,nbins,runcount;

   bool write_meta_file;
   bool sema_file_enabled_start, sema_file_enabled_done;
   string sema_file_name; // user-defined name of semaphore file, if empty: file name will be derived in function getSemaFN
};

inline string Setup::getLattice(){return lattice;}
inline double Setup::getReferenceLength(){ return lambda0; }
inline void   Setup::setReferenceLength(double lam){ lambda0=lam; return; }
inline double Setup::getReferenceEnergy(){ return gamma0; }
inline bool   Setup::getOne4One(){return one4one;}
inline bool   Setup::getShotNoise(){return shotnoise;}
inline bool   Setup::getBeamGlobalStat(){return beam_global_stat;}
inline bool   Setup::getFieldGlobalStat(){return field_global_stat;}
inline int    Setup::getNpart(){return npart;}
inline int    Setup::getNbins(){return nbins;}
inline int    Setup::getSeed(){return seed;}
inline double Setup::getStepLength(){return delz;}
inline void   Setup::setStepLength(double din){delz=din;return;}
inline int    Setup::getCount(){return(runcount);}
inline void   Setup::incrementCount(){runcount++; return;}
inline void   Setup::setRootName(string *newname){rootname=*newname; runcount=0; return;}
inline bool   Setup::outputFFT(){ return exclude_fft_output;}
inline bool   Setup::outputSpatial(){ return exclude_spatial_output;}
inline bool   Setup::outputIntensity(){ return exclude_intensity_output;}
inline bool   Setup::outputEnergy(){ return exclude_energy_output;}
inline bool   Setup::outputCurrent(){ return exclude_current_output;}
inline bool   Setup::outputAux(){ return exclude_aux_output;}
inline bool   Setup::outputFieldDump() { return exclude_field_dump;}
inline void   Setup::set_do_write_outfile(bool v) {do_write_outfile=v;}
inline bool   Setup::get_do_write_outfile()       {return do_write_outfile;}

inline bool   Setup::BWF_get_enabled()    { return beam_write_filter; }
inline void   Setup::BWF_set_enabled(bool in) { beam_write_filter=in; }
inline int    Setup::BWF_get_from()       { return beam_write_slices_from; }
inline void   Setup::BWF_set_from(int in) { beam_write_slices_from=in; }
inline int    Setup::BWF_get_to()         { return beam_write_slices_to; }
inline void   Setup::BWF_set_to(int in)   { beam_write_slices_to=in; }
inline int    Setup::BWF_get_inc()        { return beam_write_slices_inc; }
inline void   Setup::BWF_set_inc(int in)
{
	if(in<=0) {
		int r;
		MPI_Comm_rank(MPI_COMM_WORLD, &r); // unclear if member variable 'rank' is valid at any time this setter function is called
		if(r==0) {
			cout << "*** beam_write_slices_inc: positive value expected, forcing to 1" << endl;
		}
		in=1;
	}
	beam_write_slices_inc=in;
}

inline bool   Setup::getSemaEnStart() { return sema_file_enabled_start; }
inline bool   Setup::getSemaEnDone()  { return sema_file_enabled_done; }
inline void   Setup::getOutputdir(string *q) {*q = outputdir;}
inline bool   Setup::get_write_meta_file(){ return write_meta_file;}
#endif
