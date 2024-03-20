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

  loc_count_workaround=0;
}

Collective::~Collective()
{
    this->clearWake();
}

void Collective::clearWake(){
    if (hasWake){
        wake.clear();
        wakegeo.clear();
        wakeres.clear();
        wakerou.clear();
        current.clear();
        dcurrent.clear();
        wakeext.clear();
        wakeint.clear();
        count.clear();
    }
    hasWake = false;
}

void Collective::resize_and_zero(vector<double>& v, size_t n)
{
  v.resize(n);
  for(size_t j=0; j<n; j++)
    v.at(j)=0;
}
void Collective::resize_and_zero_i(vector<int>& v, size_t n)
{
  v.resize(n);
  for(size_t j=0; j<n; j++)
    v.at(j)=0;
}

void Collective::initWake(unsigned int ns_in, unsigned int nsNode, double ds_in, double *wakeext_in, double *wakeres_in, double *wakegeo_in, double * wakerou_in, double ztrans_in, double radius_in, bool transient_in)
{ 
  // obtain MPI rank/size
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
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
  ncur=static_cast<int>(size)*nsNode;   // number of slices in simulations (should be ns/sample)
  dscur=ds*ns/ncur;
  ztrans=ztrans_in;
  radius=radius_in;
  transient=transient_in;


  // array to hold current profile with simulation resolution.
  resize_and_zero(cur,ncur+1);

  resize_and_zero(wakeext, nsNode);  // global wake, explicityle defined in input deck (e.g. constant energy loss)
  resize_and_zero(wakeint, nsNode);  // wake, internally calculated when updating the wake potential (e.g. sorting)
  resize_and_zero_i(count, nsNode);

  resize_and_zero(wake,    ns);
  resize_and_zero(wakegeo, ns);
  resize_and_zero(wakeres, ns);
  resize_and_zero(wakerou, ns);
  resize_and_zero(current, ns);
  resize_and_zero(dcurrent,ns);

  // fill the wakes or single particle wakes
  for (int i=0; i <nsNode; i++){
    wakeext[i]=wakeext_in[i];
    wakeint.at(i)=0;
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
}


 

void Collective::apply(Beam *beam, Undulator *und, double delz)
{

  if (!hasWake){
    return;
  }

  if (needsUpdate || transient) { this->update(beam, und->getz()); }


  // apply wakes

  for (int ic = 0; ic <beam->current.size(); ic++){
    beam->eloss[ic]=wakeext[ic]+wakeint.at(ic);
    double dg=beam->eloss[ic]*delz/511000;    // actuall beam loss per integration step
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
    /* FIXME: range checking needed here (to avoid possible issues with floor operation resulting in idx=-1? */
    double wei=1-(s-idx*dscur)/dscur;  
    current[is]=wei*cur.at(idx)+(1-wei)*cur.at(idx+1);
    current[is]*=ds/ce;   // convert current to number of electrons
    dcurrent[is]=-(cur.at(idx+1)-cur.at(idx))*ds/ce/dscur;
    wake[is]=0;
  }
 
  //----------------------------------------   
  // step 3 - calculate the startng position for transient

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
    for (unsigned long i=icut; i < ns-is; i++){  // loop from evaluation point till the head of the bunch
      wakeloc+=current[is+i]*(wakeres[i]+wakerou[i]);
      wakeloc+=dcurrent[is+i]*wakegeo[i];
    }
    
    double floorarg = (s-sc0)/dscur;
    int idx = floor(floorarg);
    /* Lechner, 2024-Jan: workaround (range checking) from commit id 28dd1fb, use only until issue is fixed */
    if(idx<0) {
        const double eps=1.0e-10;
        loc_count_workaround++;
        if(floorarg >= -eps) {
            idx=0;
            if(1==loc_count_workaround) {
                cout << "workaround/hack in Collective.cpp: set idx=" << idx << ", floorarg=" << floorarg << endl;
            } else if (2==loc_count_workaround) {
                cout << "workaround/hack in Collective.cpp: set idx=" << idx << ", floorarg=" << floorarg << " (not reporting further interventions)" << endl;
            } else {
                /* only reporting the first two interventions */
            }
        } else {
            cout << "!!! BADSTUFF in Collective.cpp: idx=" << idx << ", floorarg=" << floorarg << endl;
            abort();
        }
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
  needsUpdate=false;
}














