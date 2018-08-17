#ifndef __GENESIS_LATTICE__
#define __GENESIS_LATTICE__

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <map>


//#include "HDF5base.h"

#include "LatticeParser.h"
#include "Optics.h"
#include "AlterLattice.h"
#include "Undulator.h"

using namespace std;

class AlterLattice;

class Lattice{
 public:
   Lattice();
   virtual ~Lattice();
   bool generateLattice(double, double, double,AlterLattice *, Undulator *);

   //   bool writeLattice(hid_t,double, double, double,AlterLattice *);

   bool parse(string,string,int,bool);
   void getMatchedOptics(double *, double *, double *, double *);
   void match(int, double, double);
 private:

   map<double,int> layout;
   vector<Element *> lat;
   int nlat;
   bool matched;
   double mbetax,mbetay,malphax,malphay;
   int findElement(double, double, string);
   int findMarker(double, string);
   void unrollLattice(double);
   void calcSlippage(double, double);

   vector<double> lat_aw,lat_qf, lat_kx, lat_ky, lat_ku, lat_slip, lat_phase,lat_ps;
   vector<double> lat_ax, lat_ay, lat_qx,lat_qy, lat_dz, lat_z;
   vector<double> lat_lb, lat_ld, lat_lt, lat_delay, lat_cx, lat_cy; // lattice for chicane and corrector elements
   vector<int>    lat_mk,lat_helical;    // array for markers and helicity
   vector<double> lat_gradx, lat_grady; // transverse gradient of undulator

};


#endif
