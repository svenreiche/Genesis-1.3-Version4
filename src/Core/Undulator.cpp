#include "Undulator.h"

Undulator::~Undulator(){}

Undulator::Undulator()
{
  istepz=-1;
  zstop=1e9;
  zfrac=1;
}


void Undulator::updateOutput(int nzout)
{


// calculate the size of the output record

  istepz=-1;
  nout=1;
  out.resize(nstepz+1);
  out[0]=true;  // first is always output
  for (int i=1;i<(nstepz+1);i++){   
    out[i]=false;
    if (((i % nzout)==0)&&( (marker[i-1] & 8) == 0)){
      out[i]=true;
      nout++;
    } 
  }

  return;
}

void Undulator::updateMarker(int nfld, int npar, int nsort, double zstop_in)
{
    nstepz=aw.size();
    zstop = zstop_in;
    zfrac = 1. ;
    if (zstop < z[nstepz-1]) {
        zfrac = zstop / z[nstepz - 1];
    }

    for (int i=0; i<marker.size();i++) {
      if (nfld > 0) {  // field dump
          if ((i % nfld) == 0) {
              marker.at(i) |= 1;
          }
      }
      if (npar > 0) {    // particle dump
          if ((i % npar) == 0) {
              marker.at(i) |= 2;
          }
      }
      if (nsort > 0) {    // sorting
          if ((i % nsort) == 0) {
              marker.at(i) |= 4;
          }
      }
    }
    for (int i=0; i<z.size();i++) {
        if (z.at(i) > zstop) {  // stop calculation
            marker.at(i)|=8;
        }
    }
    return;
}

// Field dumps at the exit of the undulator (one dump for each undulator
// in the expanded lattice)?
void Undulator::markUndExits(void)
{
  int nz=marker.size();
  for (int i=0;i<nz;i++)
  {
    // NOTE: Setting flag for next integration step (the first one with aw==0),
    // as the field is dumped before any other work done for this step

    // !!!        vvvvvvvvvvvvv=== NOTE: as of git commit 1a9d191 (2021-08-11): vector 'aw' is one element shorter than 'marker', see Lattice::generateLattice !!!
    if ( (i>0) && (i<aw.size()) ) {
      if ((aw.at(i-1)!=0) && (aw.at(i)==0)) {
        marker.at(i) |= 1; /* request field dump */
      }
    }
  }
}

void Undulator::reportLattice(string fn_report)
{
  ofstream fo;
  int nz=aw.size(); /* aw is one element shorter than marker (as of git commit 1a9d191), see Lattice::generateLattice */

  fo.open(fn_report.c_str());
  fo << "i,z,aw,qf,marker,marker_decoded" << endl;
  for (int i=0;i<nz;i++)
  {
    int m  = marker[i];
    string decoded;

    decoded = "";
    if(m&1)
      decoded+="F"; // Field dump
    if(m&2)
      decoded+="P"; // Particle dump
    if(m&4)
      decoded+="S"; // Sort
    if(m&8)
      decoded+="X"; // eXit

    fo << i << "," << z[i] << "," << aw[i] << "," << qf[i] << "," << marker[i] << "," << decoded << endl;
  }
  fo.close();
}


/*

bool Undulator::init(hid_t fid)
{

  nzout=1; 
  zstop=1e9;

  readDataDouble(fid,(char *)"/Global/gamma0",&gammaref,1);
  readDataDouble(fid,(char *)"/Global/zstop",&zstop,1);
  readDataInt(fid,(char *)"/Global/nzout",&nzout,1);
  

  int nwork=getDatasetSize(fid,(char *)"/Lattice/aw");

  aw.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/aw",&aw[0],nwork);
  helical.resize(nwork);
  readDataInt(fid,(char *)"/Lattice/helical",&helical[0],nwork);
  ax.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ax",&ax[0],nwork);
  ay.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ay",&ay[0],nwork);
  ku.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ku",&ku[0],nwork);
  kx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/kx",&kx[0],nwork);
  ky.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ky",&ky[0],nwork);
  gradx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/gradx",&gradx[0],nwork);
  grady.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/grady",&grady[0],nwork);
  qf.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/qf",&qf[0],nwork);
  qx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/qx",&qx[0],nwork);
  qy.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/qy",&qy[0],nwork);
  z.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/z",&z[0],nwork);
  dz.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/dz",&dz[0],nwork);
  slip.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/slippage",&slip[0],nwork);
  phaseshift.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/phaseshift",&phaseshift[0],nwork);
  cx.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/cx",&cx[0],nwork);
  cy.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/cy",&cy[0],nwork);
  chic_lb.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/lb",&chic_lb[0],nwork);
  chic_ld.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/ld",&chic_ld[0],nwork);
  chic_lt.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/lt",&chic_lt[0],nwork);
  chic_angle.resize(nwork);
  readDataDouble(fid,(char *)"/Lattice/delay",&chic_angle[0],nwork);
  marker.resize(nwork);
  readDataInt(fid,(char *)"/Lattice/marker",&marker[0],nwork);

  // convert delay into bending angle
  

  for (int i=0; i<aw.size(); i++){
    if (chic_angle[i]!=0){
       double delay=fabs(chic_angle[i]);
       double tmin=0;
       double tmax=asin(1)-0.001;
       bool converged=false;
       double theta,d;
       while (!converged){
         theta=0.5*(tmax+tmin);
         d=4*chic_lb[i]*(theta/sin(theta)-1)+2*chic_ld[i]*(1/cos(theta)-1);
         if (d>delay) {
           tmax=theta;
         } else {
          tmin=theta;
	 }
         if (fabs(delay-d)<1e-15) { converged=true; }
       }
       chic_angle[i]=theta;
    }
  }
 






// calculate the size of the output record

  istepz=-1;
  nstepz=aw.size();
  nout=1;
  out.resize(nstepz+1);
  out[0]=true;  // first is always output
  for (int i=1;i<(nstepz+1);i++){   
    out[i]=false;
    if (((i % nzout)==0)&&(z[i]<zstop)){
      out[i]=true;
      nout++;
    } 
  }

  return true;

}
*/


bool Undulator::advance(int rank)
{
  istepz++;
  
  if (istepz >= nstepz){  // end reached?
    return false;
  }

  if ((marker[istepz]&8)>0){ // check for 3rd bit set in marker value for stoping calculation
    if (rank==0){
      cout << "Calculation terminated due to requested stop." << endl;
    }
    return false; 
  }

  int dstepz=static_cast<int> (round(nstepz/10.*zfrac));
  if (dstepz<1){dstepz=1;}

  // Lechner, 2023-07-27: Changed the condition to display also progress of 0%
  // (to inform user that simulation is now running)
  // Remark: First execution of this code block with istepz==0
  if (((istepz % dstepz) == 0) /* && (istepz >0) */ ) {
    if (rank==0) {
      cout << "  Calculation: " <<10*istepz/dstepz << "% done" << endl;
    }
  }
  return true;
}

double Undulator::fc(int h)
{
  BesselJ bessj;

  double coup=aw[istepz];
  if (this->isHelical()){
    if (h==1){
      return coup;
    } else {
      return 0;
    }
  } else {
    double xi = aw[istepz]*aw[istepz];
    xi=0.5 * xi /(1+xi)*static_cast<double>(h);
    int h0,h1;
    if ((h % 2) == 1){
      h0=(h-1)/2;
      h1 = h0+1;
      return coup*(bessj.value(h0,xi)-bessj.value(h1,xi))*pow(-1.,h0);
    } else {
      //      h0=(h-2)/2;
      //      h1 = h0+2;
      //      return coup*0.5*(bessj.value(h0,xi)-bessj.value(h1,xi))*pow(-1.,h0); // this needs polarization enebaled - not implemented yet
      return 0;
    }
  }
}


double Undulator::faw2(double x, double y){  // square of the transverse dependence of the undulator field.
  double dx=x-ax[istepz];
  double dy=y-ay[istepz]; 
  return (1+kx[istepz]*dx*dx+ky[istepz]*dy*dy+2*(gradx[istepz]*dx+grady[istepz]*dy)); // note kx is scaled as XKX*ku*ku in Lattice.cpp, gradx as ku*GRADX.
}



double Undulator::faw(double x, double y){  // transverse dependence of the undulator field.
  double dx=x-ax[istepz];
  double dy=y-ay[istepz]; 
  return (1+0.5*(kx[istepz]*dx*dx+ky[istepz]*dy*dy)+gradx[istepz]*dx+grady[istepz]*dy); // note kx is scaled as XKX*ku*ku in Lattice.cpp, gradx as ku*GRADX.
}
