#include "Beam.h"
#include "Field.h"
#include "Sorting.h"

#include "BeamDiag.h"

Beam::~Beam(){}
Beam::Beam(){
  bd_std=NULL;
  can_change_diaghooks=true;
}

void Beam::init(int nsize, int nbins_in, double reflen_in, double slicelen_in, double s0_in, bool one4one_in )
{  

  nbins=nbins_in;
  reflength=reflen_in;  // the length corresponding to 2pi in ponderomotive phase.
  slicelength=slicelen_in;  // reflength times samplerate.
  s0=s0_in;
  one4one=one4one_in;

  current.resize(nsize);
  eloss.resize(nsize);
  for (int i=0; i<nsize; i++){
    eloss[i]=0;
  }

  beam.resize(nsize);

  return;
}

/*** Start: Code for modular beam diagnostics ***/
void Beam::register_beam_diag(BeamDiag *bd) {
	if(can_change_diaghooks)
		diaghooks.push_back(bd);
	else {
		cout << "*** Error in Beam::register_beam_diag: cannot hook additional diagnostics after calling initDiagnostics ***" << endl;
	}
}
void Beam::clear_beam_diag(void) {
	for (unsigned int k=0; k<diaghooks.size(); k++)
		delete diaghooks.at(k);

	diaghooks.clear();

	can_change_diaghooks=true;
}
void Beam::beam_diag_store_results(hid_t parentobj) {
	for(unsigned int k=0; k<diaghooks.size(); k++) {
		diaghooks.at(k)->output(parentobj);
	}
}
void Beam::beam_diag_list_registered(void) {
	cout << "Listing diagnostics currently registered in this Beam instance" << endl;
	for(unsigned int k=0; k<diaghooks.size(); k++) {
		cout << "   element " << k << ": " << *(diaghooks.at(k)) << endl;
	}
	cout << "   (end of list)" << endl;
}
/*** End: Code for modular beam diagnostics ***/

void Beam::initDiagnostics(int nz)
{
  idx=0;
  int ns=current.size();

  can_change_diaghooks=false;
  for(unsigned int k=0; k<diaghooks.size(); k++) {
    diaghooks.at(k)->init(nz,ns);
  }
}

// initialize the sorting routine
// reference position is in ponderomotive phase. The valid bucket size is from 0 to 2 pi*sample
void Beam::initSorting(int rank,int size,bool doshift,bool dosort)
{
  int isz=beam.size();
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.init(rank,size,doshift,dosort);
  sorting.configure(0,sl,0,sl*isz,0,sl*isz,false); 
  return;
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

  for (int i=0; i<field->size();i++){
    field->at(i)->setStepsize(delz);
  }

  solver.track(delz*0.5,this,und,false);   // track transverse coordinates first half of integration step

  solver.advance(delz,this,field,und);     // advance longitudinal variables 

  incoherent.apply(this,und,delz);         // apply effect of incoherent synchrotron
  col.apply(this,und,delz);         // apply effect of collective effects
  
  solver.applyR56(this,und,reflength);    // apply the longitudinalphase shift due to R56 if a chicane is selected.

  solver.track(delz*0.5,this,und,true);      // apply corrector settings and track second half for transverse coordinate
  return;


}

bool Beam::harmonicConversion(int harmonic, bool resample)
{
  if ((resample==true) && (one4one==false)) { return false;}  // resampling requires one4one
 
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

  // blowing up the slice number
  int nsize=beam.size();

  beam.resize(harmonic*nsize);
  current.resize(harmonic*nsize);
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.configure(0,sl,0,sl*nsize*harmonic,0,sl*nsize*harmonic,false);   


  // first make all new slices zero length
  for (int i=nsize; i < harmonic*nsize; i++){
    beam.at(i).resize(0);    
    current.at(i)=0;
  }
  // second copy the old slices 0,1,.. n-1 to h*0,h*1,... h*(n-1)

  Particle p;
  for (int i=nsize-1; i>0; i--){  // runs down to 1 only because there is no need to copy from 0 to h*0 = 0
    int nloc=beam[i].size();
    for (int j=0; j<nloc; j++){
      p.gamma=beam[i].at(j).gamma;
      p.theta=beam[i].at(j).theta;
      p.x    =beam[i].at(j).x;
      p.y    =beam[i].at(j).y;
      p.px   =beam[i].at(j).px;
      p.py   =beam[i].at(j).py;
      beam[i*harmonic].push_back(p);
    }
    beam[i].clear(); 
  }

  // updating the sorting algorithm

  int shift=this->sort();  // sort the particles and update current
  return true;
}

bool Beam::subharmonicConversion(int harmonic, bool resample)
{
  if ((resample==true) && (one4one==false)) { return false;}
   
 
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

  // updating the sorting algorithm
  int isz=beam.size();
  double sl=4*asin(1.)*slicelength/reflength;
  sorting.configure(0,sl,0,sl*isz,0,sl*isz,false); 


  int shift=this->sort();  // sort the particles and update current
  return true;
}


void Beam::diagnostics(bool output, double z)
{
  if (!output) { return; }


  for(unsigned int k=0; k<diaghooks.size(); k++) {
    diaghooks.at(k)->do_diag(this, z);
  }

  // Beam diagnostics complete: increment index into arrays holding the diagnostic data
  idx++;
}

