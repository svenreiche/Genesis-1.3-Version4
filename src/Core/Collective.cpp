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
    this->clearWake();
}

void Collective::clearWake(){
    /*
    if (hasWake){
        delete [] wake;
        delete [] wakegeo;
        delete [] wakeres;
        delete [] wakerou;
        delete [] current;
        delete [] dcurrent;
        delete [] wakeext;
        delete [] wakeint;
        delete [] count;
    }*/
    hasWake = false;
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

  
  ns=ns_in;  // full number of slices with highest resolution (as for sample = 1)
  ds=ds_in;  // length of a slice at highest resolution
  ncur=static_cast<int>(size)*nsNode;   // number of slices in simulations (should be ns/sample)
  dscur=ds*ns/ncur;                     // length of slice in simulation
  ztrans=ztrans_in;
  radius=radius_in;
  transient=transient_in;

 
  //  these wakes are only defined for the slices on a single node!!!
  wakeext.resize(nsNode);  // global wake, explicityle defined in input deck (e.g. constant energy loss)
  wakeint.resize(nsNode);  // wake, internally calculated when updating the wake potential (e.g. sorting)
  count.resize(nsNode);

  // array to hold current profile with simulation resolution.
  // cur     = new double[ncur+1];
  cur.resize(ncur+1);
  for(int kk=0; kk<ncur+1; kk++)
    cur[kk]=0;

  // full single particle wakes with highest reoslution.
  wake.resize(ns);     // total wake at high resolution
  wakegeo.resize(ns);    // gemoetric wake
  wakeres.resize(ns);    // resistive wall
  wakerou.resize(ns);    // roughness wake
  current.resize(ns);    // current profile in high resolution
  dcurrent.resize(ns);    // current differential

  // fill the wakes or single particle wakes
  for (int i=0; i <nsNode; i++){
    wakeext.at(i)=wakeext_in[i];   // external wake - user defined
    wakeint.at(i)=0;               // internal wake - calculated from single particle wakes
  }


  for (int i=0; i <ns; i++){
    wakegeo.at(i)=wakegeo_in[i];
    wakeres.at(i)=wakeres_in[i];
    wakerou.at(i)=wakerou_in[i];
  }
  wakegeo.at(0)*=0.5;  // self-loading theorem
  wakeres.at(0)*=0.5;  // self-loading theorem
  wakerou.at(0)*=0.5;

  hasWake=true;
  needsUpdate=true;
}


 

void Collective::apply(Beam *beam, Undulator *und, double delz)
{

  if (!hasWake){
    return;
  }

  if (needsUpdate || transient) { this->update(beam, und->getz()); }


  // apply wakes

  for (int ic = 0; ic <beam->current.size(); ic++){
    beam->eloss[ic]=wakeext.at(ic)+wakeint.at(ic);
    double dg=beam->eloss[ic]*delz/511000;    // actual beam loss per integration step
    unsigned long npart=beam->beam.at(ic).size();
    for (unsigned long ip=0; ip<npart; ip++){
      beam->beam.at(ic).at(ip).gamma+=dg;
    }
  }
}



void Collective::update(Beam *beam, double zpos)
{  
  // ---------------
  // step 0 - checks
  int nsNode=static_cast<int>(beam->current.size());
  // the following access with 'at' deterministically crashes the program before the MPI_Allgather operation if destination array is too small
  cur.at(ncur)=0; // used for interpolation

  // check if the total number of slices used to prepare global current buffer has changed
  if(nsNode*size!=ncur)
    abort();
  
  // ---------------------------
  // step 1 - gather current profile from all nodes and stored in the vector cur

  MPI_Allgather(
     &beam->current[0],nsNode,MPI_DOUBLE,
     cur.data(),       nsNode,MPI_DOUBLE,
     MPI_COMM_WORLD);
  cur.at(ncur)=0;   // used for interpolation


 
  //-------------------------------------------------
  // step 2 - interpolate current profile to high resolution and initialize the total wake function
   
  for (int is=0; is <ns ; is++){
    double s=ds*static_cast<double> (is);
    unsigned int idx=static_cast<int> (floor(s/dscur));
    double wei=1-(s-idx*dscur)/dscur;  
    current.at(is)=wei*cur.at(idx)+(1-wei)*cur.at(idx+1);
    current.at(is)*=ds/ce;   // convert current to number of electrons
    dcurrent.at(is)=-(cur.at(idx+1)-cur.at(idx))*ds/ce/dscur;
    wake.at(is)=0;
  }
 
  //----------------------------------------   
  // step 3 - calculate the starting position for transient

  unsigned long icut;
  double z=zpos+ztrans;  // effective length from first source point
  double delta=0.5*radius*radius;
  if (z <=0) {      // if value is negative no wakefields are calculated
      icut=ns;
  } else {
      icut=static_cast<unsigned long>(floor(delta/z/ds));
  }
  if (!transient){
     icut=0;
  }
 

  //--------------------------------------------
  // step 4 - initialize the loss vector for the given time-window slice of the given node


  // --- new code to avoid some nasty rounding errors
  auto sample = static_cast<int>(dscur/ds);
  if (sample < 1) {sample = 1;} // steps per slice
  for (int ic=0; ic <nsNode;ic++) { // looping over the internal slices
      wakeint.at(ic)=0;
       auto is0 = (nsNode*rank+ic)*sample;
       for (int j=0; j < sample; j++){
            auto is = is0+j;
            for (unsigned long i=0; i < ns-is; i++){  // loob from evaluation point till the head of the bunch
               wakeint.at(ic)+=current.at(is+i)*(wakeres.at(i)+wakerou.at(i));
               wakeint.at(ic)+=dcurrent.at(is+i)*wakegeo.at(i);
            }
       }
       wakeint.at(ic)/=static_cast<double>(sample);  // since the current is weighted multiple times when looping over sample
  }



  /*
  for (int ic = 0; ic <nsNode; ic++){
    count.at(ic)=0;
    wakeint.at(ic)=0;
  }

  double sc0=dscur*(nsNode*rank);     // s positions of the timewindow slice for a given rank.
  double sc1=dscur*(nsNode*(rank+1));
  int is0=static_cast<int>(round(sc0/ds));
  int is1=static_cast<int>(round(sc1/ds));

  for (int is=is0; is< is1; is++){  // loop the evaluation point from back to front
    double s=is*ds;
    double wakeloc=0;
    for (unsigned long i=icut; i < ns-is; i++){  // loob from evaluation point till the head of the bunch
      wakeloc+=current.at(is+i)*(wakeres.at(i)+wakerou.at(i));
      wakeloc+=dcurrent.at(is+i)*wakegeo.at(i);
    }
    int idx = floor((s-sc0)/dscur);
    if ((idx < 0) or (idx >= count.size())){
        cout << "Idx out of range " << idx << " rank: " << rank << endl;
        cout << "s: " << s << endl;
        cout << "sc0: " << sc0 << endl;
        cout << "dscur: " << dscur << endl;
    }
    count.at(idx)++;
    wakeint.at(idx)+=wakeloc;
  }

  for (int ic = 0; ic <nsNode; ic++){
    if (count.at(ic) > 0) {
      wakeint.at(ic)/=static_cast<double>(count.at(ic));
    } else {
      wakeint.at(ic)=0;
    }
  }
  */

  needsUpdate=false;
}














