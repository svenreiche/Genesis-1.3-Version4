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

	// one instance per obj. in the file
	HDF5_CollWriteCore *pcwc;
	HDF5_CollWriteCore *pcwc_filt;

	// parameters describing the data layout along the first axis of the HDF5 data objs to be generated (=longitudinal axis)
	int curr_slice;
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
   void initSourceFilter_DbgDumpSettings(bool, int, string);

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
   int crsource_dump_every_;
   string crsource_dump_rootname_;
   int call_cntr_adv_;
   void dump_crsource(struct dump_settings *, HDF5_CollWriteCore *);

   void filterSourceTerm(struct dump_settings *);
   void ADI(vector<complex< double > > &);
   void tridagx(vector<complex< double > > &);
   void tridagy(vector<complex< double > > &);
};

#endif
