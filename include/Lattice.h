#ifndef __GENESIS_LATTICE__
#define __GENESIS_LATTICE__

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <map>


//#include "HDF5base.h"
#include "RandomU.h"
#include "LatticeParser.h"
#include "Optics.h"
#include "Setup.h"
#include "AlterLattice.h"
#include "Inverfc.h"
#include "Undulator.h"

using namespace std;

class AlterLattice;
class Series;
class SeriesManager;

class Lattice{
 public:
   Lattice();
   virtual ~Lattice();
   bool generateLattice(Setup *, AlterLattice *, Undulator *);

   //   bool writeLattice(hid_t,double, double, double,AlterLattice *);

   bool parse(string,string,int,SeriesManager *);
   void getMatchedOptics(double *, double *, double *, double *);
   void match(int, double, double);

   bool alterElement(string,string, double, string, SeriesManager *,int,bool);

   void report(string);

 private:

    map<double,int> layout;
    vector<Element *> lat;
    int nlat{};
    bool matched;
    double mbetax{},mbetay{},malphax{},malphay{};
    int findElement(double, double, string);
    int findMarker(double, string);
    void unrollLattice(double,double, double,bool,unsigned long);
    void calcSlippage(double, double);
    void generateFieldErrors(double aw0,vector<double> *awout, double awerr);
    void generateKickErrors( vector<double> *awout, vector<double> *kicks, double delz, double gammaref);
    Inverfc erf;
    RandomU *seq;
    vector<double> lat_aw,lat_qf, lat_kx, lat_ky, lat_ku, lat_slip, lat_phase,lat_ps;
    vector<double> lat_ax, lat_ay, lat_qx,lat_qy, lat_dz, lat_z;
    vector<double> lat_lb, lat_ld, lat_lt, lat_delay, lat_cx, lat_cy; // lattice for chicane and corrector elements
    vector<int>    lat_mk,lat_helical;    // array for markers and helicity
    vector<double> lat_gradx, lat_grady; // transverse gradient of undulator

};


#endif
