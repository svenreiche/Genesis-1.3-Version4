#include "Beam.h"
#include "Field.h"
#include "Sorting.h"

#include <mpi.h>

extern bool MPISingle;

Beam::~Beam()= default;
Beam::Beam(){
      do_global_stat=false;
      doCurrent=false;
      doSpatial=true;
      doEnergy=true;
      doAux=true;
      beam_write_filter=false;
      beam_write_slices_from=-1;
      beam_write_slices_to=-1;
      beam_write_slices_inc=1;
}

void Beam::init(int nsize, int nbins_in, double reflen_in, double slicelen_in, double s0_in, bool one4one_in ) {

    nbins = nbins_in;
    reflength = reflen_in;  // the length corresponding to 2pi in ponderomotive phase.
    slicelength = slicelen_in;  // reflength times samplerate.
    s0 = s0_in;
    one4one = one4one_in;
    do_global_stat = false;

    current.resize(nsize);
    eloss.resize(nsize);
    longESC.resize(nsize);  // array to hold long range space charge field
    for (int i = 0; i < nsize; i++) {
        eloss[i] = 0;
        longESC[i] = 0;
    }
    beam.resize(nsize);
}


double Beam::getSize(int is) {  // a calculation of the rms sizes needed for space charge calculation
    if (this->current[is] <= 0){
        return 1;    // dummy value if there is no current
    }

    double x1=0;
    double x2=0;
    double y1=0;
    double y2=0;
    int count = 0;
    for (auto const &par: this->beam.at(is)){
            x1 += par.x;
            x2 += par.x*par.x;
            y1 += par.y;
            y2 += par.y*par.y;
            count++;
    }
    if (count == 0) {
        return 1.;
    }
    x1 /=static_cast<double>(count);
    x2 /=static_cast<double>(count);
    y1 /=static_cast<double>(count);
    y2 /=static_cast<double>(count);
    return sqrt(std::abs(x2-x1*x1))*sqrt(std::abs(y2-y1*y1));
}



// initialize the sorting routine
// reference position is in ponderomotive phase. The valid bucket size is from 0 to 2 pi*sample
void Beam::initSorting(int rank,int size,bool doshift,bool dosort)
{
  auto isz=beam.size();
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.init(rank,size,doshift,dosort);
  sorting.configure(0,sl,0,sl*isz,0,sl*isz,false); 
}

int Beam::sort()
{
  int shift=0;
  double dQ=ce/slicelength;
  if (one4one){
    shift= sorting.sort(&beam);
    for (int i=0; i<beam.size();i++){    // correct the local current
      int np=beam.at(i).size();
      current.at(i)=static_cast<double>(np)*ce/slicelength;
    }
    col.forceUpdate();
  }  
  return shift;
}


void Beam::track(double delz,vector<Field *> *field, Undulator *und){

  for (auto & ifld : *field){
    ifld->setStepsize(delz);
  }

  solver.track(delz*0.5,this,und,false);   // track transverse coordinates first half of integration step

  solver.advance(delz,this,field,und);     // advance longitudinal variables 

  incoherent.apply(this,und,delz);         // apply effect of incoherent synchrotron
  col.apply(this,und,delz);         // apply effect of collective effects

  solver.applyR56(this,und,reflength);    // apply the longitudinal phase shift due to R56 if a chicane is selected.

  solver.track(delz*0.5,this,und,true);      }



void Beam::report_storage(string infotxt)
{
  bool do_report=false;
  const int report_rank=30;
  int rank;

  if(!do_report)
    return;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* display infos only on a single rank */
  if(rank==report_rank) {
    cout << "*** Reporting beam storage on rank " << report_rank << " (" << infotxt << ") ***" << endl;
    for(int k=0; k<beam.size(); k++) {
      cout << k << " size=" << beam[k].size() << " capacity=" << beam[k].capacity() << endl;
    }
    cout << "End of report on rank " << report_rank << endl;
  }

  unsigned long long glbl_size=0, glbl_capacity=0;
  unsigned long long my_tot_size=0, my_tot_capacity=0;
  for(int k=0; k<beam.size(); k++) {
    my_tot_size    +=beam[k].size();
    my_tot_capacity+=beam[k].capacity();
  }
  MPI_Allreduce(&my_tot_size, &glbl_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&my_tot_capacity, &glbl_capacity, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(rank==report_rank) {
    cout << "Totalized over all ranks: size=" << glbl_size << ", capacity=" << glbl_capacity << endl;
    cout << "*** End of report ***" << endl;
  }
}
bool Beam::dbg_skip_shrink()
{
  const char *q = getenv("G4_TEST_NOSHRINK");
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(q!=nullptr) {
    if(rank==0) {
      cout << "class Beam: skipping shrinking of vectors" << endl;
    }
    return(true);
  }
  return(false);
}

void Beam::make_compact()
{
  // Controls amount of additional memory to be reserved.
  // This ensures that we don't fall back to STL's "geometric resizing"
  // when the first additional element is inserted (for instance when
  // sorting is on during the subsequent '&track').
  // For libstdc++, this would immediately double the capacity.
  const double extra=0.05;

  if(dbg_skip_shrink())
    return;

  for(int k=0; k<beam.size(); k++) {
    beam[k].shrink_to_fit();

    unsigned long long newcap = beam[k].capacity();
    newcap += (1 + extra*newcap);
    //         ^ ensures that we get a 'ceil'
    beam[k].reserve(newcap);
  }
}

bool Beam::harmonicConversion(int harmonic, bool resample)
{
  if (resample && !one4one) { return false;}  // resampling requires one4one
 
  reflength=reflength/static_cast<double>(harmonic);
  if (resample){ 
    slicelength=slicelength/static_cast<double>(harmonic);
  }
  for (int i=0; i<beam.size(); i++){
    for (int j=0; j<beam[i].size();j++){
      beam[i].at(j).theta*=static_cast<double>(harmonic); 
    }
  }
  if (!resample) { return true; }

  col.clearWake();  // clear the wake definitions. Needs an explicit wake commando in input deck
  report_storage("before harmonic upconversion");

  // blowing up the slice number
  int nsize=beam.size();

  beam.resize(harmonic*nsize);
  current.resize(harmonic*nsize);
  eloss.resize(harmonic*nsize);
  longESC.resize(harmonic*nsize);
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.configure(0,sl,0,sl*nsize*harmonic,0,sl*nsize*harmonic,false);   


  // first make all new slices zero length
  for (int i=nsize; i < harmonic*nsize; i++){
    beam.at(i).resize(0);    
    current.at(i)=0;
  }


  // second move the old slices 0,1,.. n-1 to h*0,h*1,... h*(n-1)
  for (int i=nsize-1; i>0; i--){  // runs down to 1 only because there is no need to copy from 0 to h*0 = 0

    // By ensuring that the destination slice is empty, we can be
    // sure that the source slice is empty after the swap operation
    if(!beam.at(i*harmonic).empty()) {
      abort();
    }
    std::swap(beam.at(i*harmonic), beam.at(i));
  }

  // updating the sorting algorithm
  int shift=this->sort();  // sort the particles and update current

  report_storage("after harmonic upconversion");
  make_compact();
  report_storage("after shrinking operation");

  return true;
}

bool Beam::subharmonicConversion(int harmonic, bool resample)
{
  if (resample && !one4one) { return false;}
   
 
  reflength=reflength*static_cast<double>(harmonic);
  if (resample){ // needs a lot of working here............
    slicelength=slicelength*static_cast<double>(harmonic);
  }


  for (int i=0; i<beam.size(); i++){
    for (int j=0; j<beam[i].size();j++){
      beam[i].at(j).theta/=static_cast<double>(harmonic);   // preparing to push everything into first slice
    }
  }
  if (!resample) { return true; }


  col.clearWake();  // clear the wake definitions. Needs an explicit wake commando in input deck

// prepare to copy everyting into the first slice
  int nsize=beam.size();
  Particle p;

  if ((nsize % harmonic) !=0) { return false;}  // check whether the number of slices cannot merged into a smaller number

  //return true;

  double dtheta=4.*asin(1)*slicelength/reflength/static_cast<double>(harmonic);

  for (int i=1; i<nsize;i++){
      for (int k=0; k<beam.at(i).size();k++){
	p.gamma=beam[i].at(k).gamma;
	p.theta=beam[i].at(k).theta+i*dtheta;
	p.x    =beam[i].at(k).x;
	p.y    =beam[i].at(k).y;
	p.px   =beam[i].at(k).px;
	p.py   =beam[i].at(k).py;
	beam[0].push_back(p);
      }
      beam[i].clear(); 
  }
  beam.resize(nsize/harmonic);
  current.resize(nsize/harmonic);
  eloss.resize(nsize/harmonic);
  longESC.resize(nsize/harmonic);
  // updating the sorting algorithm
  int isz=beam.size();
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.configure(0,sl,0,sl*isz,0,sl*isz,false); 


  int shift=this->sort();  // sort the particles and update current
  return true;
}





