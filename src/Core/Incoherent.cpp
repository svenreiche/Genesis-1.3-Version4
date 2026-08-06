#include "Incoherent.h"
#include "Beam.h"

Incoherent::Incoherent(){
  base=0;
  noff=-1;
  doLoss=false;
  doSpread=false;
}

Incoherent::~Incoherent(){}

void Incoherent::init(int base_in, bool doLoss_in,bool doSpread_in)
{

  doLoss=doLoss_in;
  doSpread=doSpread_in;
  base=static_cast<unsigned int>(base_in);

  sran.clear();
  noff=-1;
  return;
}


// One generator per slice, keyed on the global slice index. The alternative,
// a single generator per core drawing for one slice after another, ties the
// numbers a slice receives to the core layout: both the seed and the position
// in the sequence depend on which core owns the slice and on how many slices
// preceded it there.
void Incoherent::ensureStreams(const Beam *beam)
{
  size_t nslice=beam->beam.size();
  if ((sran.size()==nslice)&&(noff==beam->noff)) { return; }

  noff=beam->noff;
  sran.resize(nslice);
  for (size_t i=0; i<nslice; i++){
    sran[i].set(seedFromIndex(base,static_cast<unsigned long>(noff)+i,SeedStream::incoherent));
  }
}



void Incoherent::apply(Beam *beam, Undulator *und, double delz)
{  

  if (!und->inUndulator()) { return; }
  if ((!doLoss) && (!doSpread)) { return; }

  double gam0=und->getGammaRef();
  double awz=und->getaw();
  double xkw0=und->getku();


  double dgamsig=1.015e-27* xkw0 * xkw0 * awz * awz;

  if (und->isHelical()){
    dgamsig*= 1.42 *awz + 1./(1.+1.5*awz+0.95*awz*awz);
  } else {
    dgamsig*= 1.697*awz + 1./(1.+1.88*awz+0.8*awz*awz); 
  }

  if (!doSpread){ dgamsig=0;}

  dgamsig=sqrt(dgamsig*gam0*gam0*gam0*gam0*xkw0*delz)*sqrt(3.);


  double dgamavg=xkw0*gam0*awz;
  if(!doLoss) { dgamavg=0;}

  dgamavg=dgamavg*dgamavg*1.88e-15*delz;


  // apply energy change to electorn bunch
  this->ensureStreams(beam);

  int nbins=beam->nbins;
  if (beam->one4one){ nbins=1;}
  double dg=0;

  for (int islice=0;islice< beam->beam.size();islice++){
    int npart=beam->beam.at(islice).size();
    RandomU &seq=sran[islice];
    for (int ip=0; ip<npart; ip++){
      if ((ip % nbins) == 0){
         dg=-dgamavg+dgamsig*(2*seq.getElement()-1);
      }
      beam->beam.at(islice).at(ip).gamma+=dg;
    }
  }


  return;

}

