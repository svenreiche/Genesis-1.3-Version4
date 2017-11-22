#include "ShotNoise.h"

ShotNoise::ShotNoise(){
  sran=NULL;
  nwork=1000;
  work=new double [nwork];
}

ShotNoise::~ShotNoise()
{
  delete [] work;
}



void ShotNoise::init(int base,int rank)
{
  RandomU rseed(base);
  double val;
  for (int i=0; i<=rank;i++){
    val=rseed.getElement();
  }
  val*=1e9;
  int locseed=static_cast<int> (round(val));
  if (sran !=NULL) { delete sran; }
  sran  = new RandomU (locseed);
}

void ShotNoise::applyShotNoise(Particle *beam, int npart, int nbins, double ne)
{


  if (npart>nwork){
    delete [] work;
    nwork=npart;
    work = new double [nwork];
  }

  int mpart=npart/nbins;

  double nbl=ne/static_cast<double>(mpart);  // number of simulated electrons per beamlet 

  for (int i=0; i< npart; i++){
    work[i]=0;
  }
  for (int ih=0; ih< (nbins-1)/2; ih++){
    for (int i1=0;i1<mpart;i1++){
      double phi=4*asin(1)*sran->getElement();
      double an=sqrt(-log(sran->getElement())/nbl)*2/static_cast<double>(ih+1);
      for (int i2=0; i2<nbins; i2++){
        int idx=i1*nbins+i2;
        work[idx]-=an*sin(beam[idx].theta*(ih+1)+phi);
      }
    }
  }
  
  for (int i=0; i< npart; i++){
    beam[i].theta+=work[i];
  }

  return;
}
