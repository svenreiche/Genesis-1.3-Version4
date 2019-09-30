// include all routine which are arising from a collective contribution from the beam. This are:
// Wake
// Longrange Space Charge
// CSR

#include "Collective.h"
#include "Beam.h"

Collective::Collective()
{
  hasWake=false;
}

Collective::~Collective()
{
}


void Collective::initWake(unsigned int ns_in, double ds_in, double *wakeres_in, double *wakegeo_in, double *wakerou_in, double ztrans_in, double radius_in, bool trans_in){
  ztrans=ztrans_in;
  transient=trans_in;
  radius=radius_in;
  ns=ns_in;
  ds=ds_in;
  wake    =new double[ns];
  wakeres =new double[ns];    // holds the pre-calculated single wake potential
  wakegeo =new double[ns];    // holds the pre-calculated single wake potential
  wakerou =new double[ns];    // holds the pre-calculated single wake potential
  current =new double[ns];
  dcurrent=new double[ns];    // is the differential in the current
  curwork =new double[ns+1];
  
  for (int i=0; i<ns;i++){
    wakeres[i]=wakeres_in[i];
    wakegeo[i]=wakegeo_in[i];
    wakerou[i]=wakerou_in[i];
  }
  wakeres[0]*=0.5;  // self-loading theorem
  wakegeo[0]*=0.5;  // self-loading theorem
  wakerou[0]*=0.5;  // self-loading theorem

  hasWake=true;
  return;
}



void Collective::apply(Beam *beam, Undulator *und, double delz)
{  
  if (!hasWake){
    return;
  }

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node
  MPI_Comm_size(MPI_COMM_WORLD, &size); // assign rank to node

  if (MPISingle){
      rank=0;
      size=1;
  }
  
  if (rank==0){
    cout << "applying wakes" << endl;
  }
   
  // get full current profile
  int ncur=beam->current.size();
  double dscur=beam->slicelength;
  cout << "Rank: " << rank << " ncur:" << ncur << endl;

  double *cur=new double [ncur*size+1];
  for (int i=0; i < ncur*size;i++){cur[i]=0;}

  //  MPI::COMM_WORLD.Allgather(&beam->current[0],ncur,MPI::DOUBLE, &cur[0],ncur,MPI::DOUBLE);
  MPI_Allgather(&beam->current[0],ncur,MPI_DOUBLE,&cur[0],ncur,MPI_DOUBLE,MPI_COMM_WORLD);
  ncur=ncur*size;
  cur[ncur]=0;  // needed for interpolation

  double sam=dscur/ds;
  if (rank == 0) { cout << "Test: " << sam << endl;}



  // interpolate current profile to high resolution and initialize the totak wake function
   
  for (int is=0; is <ns ; is++){
    double s=ds*static_cast<double> (is);
    unsigned int idx=static_cast<int> (floor(s/dscur));
    double wei=1-(s-idx*dscur)/dscur;  
    current[is]=wei*cur[idx]+(1-wei)*cur[idx+1];
    current[is]*=ds/ce;   // convert current to number of electrons
    dcurrent[is]=-(cur[idx+1]-cur[idx])/ce/sam;
    wake[is]=0;
  }

  if (rank==0) { cout << "Current is interpolated" << endl ;}
  

  // calculate the startng position for transient
  int icut=0;
  double z=und->getz()+ztrans;  // effective length from first source point
  double delta=0.5*radius*radius;
  if (z <=0) {      // if value is negative no wakefields are calculated
    icut=ns;
  } else {
    icut=static_cast<int>(floor(delta/z/ds));
  }
  if (!transient){
    icut=0;
  }
 
   if (rank==0) { cout << "Transient cutoff done" << endl ;}


  // do the convolution
  for (int is=0; is< ns; is++){  // loop the evaluation point from back to front
    for (int i=icut; i < ns-is; i++){  // loob from evaluation point till the head of the bunch
      wake[is]+=current[is+i]*(wakeres[i]+wakerou[i]);
      wake[is]+=dcurrent[is+i]*wakegeo[i];
    }
  }
 
  if (rank==0) { cout << "Convolution calculated" << endl ;}



  ncur=beam->current.size();
  int ioff = rank*ncur;


  double sc0=dscur*(ioff);     // s positions of the timewindow slice for a given rank.
  double sc1=dscur*(ncur+ioff);
  int *count=new int [ncur];

  for (int ic = 0; ic <ncur; ic++){
    count[ic]=0;
    beam->eloss[ic]=0;
  }

  if (rank==0) { cout << "Loss vector initialized" << endl ;}


  for (int is=0;is < ns; is++){
    double s=is*ds;
    if ((s >= sc0) and (s < sc1)){
      int idx = floor((s-sc0)/dscur);
      count[idx]++;
      beam->eloss[idx]+=wake[is];
    }
  }
  if (rank==0) { cout << "Loss calculated" << endl ;}

  
  // save in beam->eloss for output file

  for (int ic = 0; ic <ncur; ic++){
    if (count[ic]>0) {
      beam->eloss[ic]/=static_cast<double>(count[ic]);
    } else {
      beam->eloss[ic]=0;
    }
    double dg=beam->eloss[ic]*delz/511000;    // actuall beam loss per integration step
    int npart=beam->beam.at(ic).size();
    for (int ip=0; ip<npart; ip++){
      beam->beam.at(ic).at(ip).gamma+=dg;
    }
  }

  if (rank==0) { cout << "Loss applied" << endl ;}

  delete[] cur;
  delete[] count;
  return;

}
