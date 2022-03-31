#ifndef __GENESIS_FIELDSOLVER__
#define __GENESIS_FIELDSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>

#ifdef FFTW
#include <fftw3.h>
#endif

#include <hdf5.h>


class Field;
class Beam;

#include "Particle.h"
#include "Undulator.h"


using namespace std;

class HDF5_CollWriteCore;
struct dump_settings {
	bool do_dump;

	hid_t fid;

	// one instance per obj. in the file
	HDF5_CollWriteCore *pcwc;
	HDF5_CollWriteCore *pcwc_filt;

	// parameters describing the data layout along the first axis of the HDF5 data objs to be generated (=longitudinal axis)
	int curr_slice;
	int smin;
	int nstot;
};

class FieldSolver{
 public:
   FieldSolver();
   virtual ~FieldSolver();
   void getDiag(double,double,double,int);
   void advance(double, Field *, Beam *, Undulator *);
   void init(int);
   void initSourceFilter(bool, double, double, double);
   bool initSourceFilter_DbgDumpSettings(bool, int, string, string);

 private:
   int ngrid;
   double delz_save;
   complex<double> cstep;
   vector< complex< double > > r,c,cbet,cwet,crsource;

   bool   difffilter_;
   vector<double> sigmoid_;

#ifdef FFTW
   bool hasPlan;
   complex<double> *in, *out;
   fftw_plan p,pi;
#endif

   bool crsource_dump_en_;
   bool crsource_dump_filter_, crsource_dump_crsource_, crsource_dump_crsource_filt_;
   int crsource_dump_every_;
   string crsource_dump_rootname_;
   int call_cntr_adv_;
   void dump_file_open(struct dump_settings *, Field *);
   void dump_file_close(struct dump_settings *);
   void dump_crsource(struct dump_settings *, HDF5_CollWriteCore *);
   void dump_filter(struct dump_settings *, hid_t);

   void filterSourceTerm(struct dump_settings *);
   void ADI(vector<complex< double > > &);
   void tridagx(vector<complex< double > > &);
   void tridagy(vector<complex< double > > &);
};

#endif
