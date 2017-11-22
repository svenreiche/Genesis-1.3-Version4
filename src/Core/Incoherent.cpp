#include "Incoherent.h"
#include "Beam.h"

Incoherent::Incoherent(){
  sran=NULL;
  doLoss=false;
  doSpread=false;
}

Incoherent::~Incoherent(){}

void Incoherent::init(int base, int rank, bool doLoss_in,bool doSpread_in)
{

  doLoss=doLoss_in;
  doSpread=doSpread_in;


  RandomU rseed(base);
  double val;
  for (int i=0; i<=rank;i++){
    val=rseed.getElement();
  }
  val*=1e9;
  int locseed=static_cast<int> (round(val));
  if (sran !=NULL) { delete sran; }
  sran  = new RandomU (locseed);
  return;
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
  int nbins=beam->nbins;
  if (beam->one4one){ nbins=1;}
  double dg=0;

  for (int islice=0;islice< beam->beam.size();islice++){
    int npart=beam->beam.at(islice).size();
    for (int ip=0; ip<npart; ip++){
      if ((ip % nbins) == 0){
         dg=-dgamavg+dgamsig*(2*sran->getElement()-1);
      }
      beam->beam.at(islice).at(ip).gamma+=dg;
    }
  }


  return;

}

