// include all routine which are arising from a collective contribution from the beam. This are:
// Wake
// Longrange Space Charge
// CSR

#include "Collective.h"
#include "Beam.h"
#include "Undulator.h"

Collective::Collective()
{
  hasWake=false;
  transient=false;
  ztrans=0;
  radius=1.0;
}

Collective::~Collective()
{
}


void Collective::initWake(unsigned int ns_in, unsigned int nsNode, double ds_in, double *wakeext_in, double *wakeres_in, double *wakegeo_in, double * wakerou_in, double ztrans_in, double radius_in, bool transient_in)
{ 


  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node
  MPI_Comm_size(MPI_COMM_WORLD, &size); // assign ranksize to node

  if (MPISingle){
      rank=0;
      size=1;
  }

#if 0
  if(0==rank) {
    cout << "in Collective::initWake"<<endl;
  }
#endif

  
  ns=ns_in;  // full number of slices with highest resolution
  ds=ds_in;
  ncur=size*nsNode;   // number of slices in simulations (should be ns/sample)
  dscur=ds*ns/ncur;
  ztrans=ztrans_in;
  radius=radius_in;
  transient=transient_in;

 

  wakeext = new double[nsNode];  // global wake, explicityle defined in input deck (e.g. constant energy loss)
  wakeint = new double[nsNode];  // wake, internally calculated when updating the wake potential (e.g. sorting)
  count   = new int [nsNode];

  // array to hold current profile with simulation resolution.
  // cur     = new double[ncur+1];
  cur.resize(ncur+1);
  for(int kk=0; kk<ncur+1; kk++)
    cur[kk]=0;


  wake    = new double[ns];
  wakegeo = new double[ns];
  wakeres = new double[ns];
  wakerou = new double[ns];
  current = new double[ns];
  dcurrent= new double[ns];

  // fill the wakes or single particle wakes
  for (int i=0; i <nsNode; i++){
    wakeext[i]=wakeext_in[i];
    wakeint[i]=0;
  }


  for (int i=0; i <ns; i++){
    wakegeo[i]=wakegeo_in[i];
    wakeres[i]=wakeres_in[i];
    wakerou[i]=wakerou_in[i];
  }
  wakegeo[0]*=0.5;  // self-loading theorem
  wakeres[0]*=0.5;  // self-loading theorem
  wakerou[0]*=0.5;

  hasWake=true;
  needsUpdate=true;
  return;
}


 

void Collective::apply(Beam *beam, Undulator *und, double delz)
{

  if (!hasWake){
    return;
  }

  if (needsUpdate || transient) { this->update(beam, und->getz()); }


  // apply wakes

  for (int ic = 0; ic <beam->current.size(); ic++){
    beam->eloss[ic]=wakeext[ic]+wakeint[ic];
    double dg=beam->eloss[ic]*delz/511000;    // actuall beam loss per integration step
    int npart=beam->beam.at(ic).size();
    for (int ip=0; ip<npart; ip++){
      beam->beam.at(ic).at(ip).gamma+=dg;
    }
  }

  return;



}



void Collective::update(Beam *beam, double zpos)
{  
  // ---------------
  // step 0 - checks
  int nsNode=beam->current.size();
  // the following access with 'at' deterministically crashes the program before the MPI_Allgather operation if destination array is too small
  cur.at(ncur)=0; // used for interpolation

  // check if the total number of slices used to prepare global current buffer has changed
  if(nsNode*size!=ncur)
    abort();
  
  // ---------------------------
  // step 1 - gather current profile from all nodes
  MPI_Allgather(
     &beam->current[0],nsNode,MPI_DOUBLE,
     cur.data(),       nsNode,MPI_DOUBLE,
     MPI_COMM_WORLD);
  cur[ncur]=0;   // used for interpolation


 
  //-------------------------------------------------
  // step 2 - interpolate current profile to high resolution and initialize the total wake function
   
  for (int is=0; is <ns ; is++){
    double s=ds*static_cast<double> (is);
    unsigned int idx=static_cast<int> (floor(s/dscur));
    double wei=1-(s-idx*dscur)/dscur;  
    current[is]=wei*cur[idx]+(1-wei)*cur[idx+1];
    current[is]*=ds/ce;   // convert current to number of electrons
    dcurrent[is]=-(cur[idx+1]-cur[idx])*ds/ce/dscur;
    wake[is]=0;
  }
 
  //----------------------------------------   
  // step 3 - calculate the startng position for transient

  int icut=0;
  double z=zpos+ztrans;  // effective length from first source point
  double delta=0.5*radius*radius;
  if (z <=0) {      // if value is negative no wakefields are calculated
      icut=ns;
  } else {
      icut=static_cast<int>(floor(delta/z/ds));
  }
  if (!transient){
     icut=0;
  }
 

  //--------------------------------------------
  // step 4 - initialize the loss vector for the given time-window slice of the given node


  for (int ic = 0; ic <nsNode; ic++){
    count[ic]=0;
    wakeint[ic]=0;
  }



  double sc0=dscur*(nsNode*rank);     // s positions of the timewindow slice for a given rank.
  double sc1=dscur*(nsNode*(rank+1));
  int is0=static_cast<int>(round(sc0/ds));
  int is1=static_cast<int>(round(sc1/ds));

  for (int is=is0; is< is1; is++){  // loop the evaluation point from back to front
    double s=is*ds;
    double wakeloc=0;
    for (int i=icut; i < ns-is; i++){  // loob from evaluation point till the head of the bunch
      wakeloc+=current[is+i]*(wakeres[i]+wakerou[i]);
      wakeloc+=dcurrent[is+i]*wakegeo[i];
    }
    int idx = floor((s-sc0)/dscur);
    count[idx]++;
    wakeint[idx]+=wakeloc;
  }

  for (int ic = 0; ic <nsNode; ic++){
    if (count[ic] > 0) {
      wakeint[ic]/=static_cast<double>(count[ic]);
    } else {
      wakeint[ic]=0;
    }
  }

  needsUpdate=false;
  return;
   
}














