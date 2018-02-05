#include "BeamSolver.h"
#include "Field.h"
#include "Beam.h"

BeamSolver::BeamSolver()
{
  onlyFundamental=false;
}

BeamSolver::~BeamSolver(){}


void BeamSolver::advance(double delz, Beam *beam, vector< Field *> *field, Undulator *und)
{
   
  // here the harmonics needs to be taken into account

  vector<int> nfld;
  vector<double> rtmp;
  rpart.clear();
  rharm.clear();
  xks=1;  // default value in the case that no field is defined

   for (int i=0; i < field->size(); i++){
    int harm=field->at(i)->getHarm();
    if ((harm==1) || !onlyFundamental){
      xks=field->at(i)->xks/static_cast<double>(harm);    // fundamental field wavenumber used in ODE below
      nfld.push_back(i);
      rtmp.push_back(und->fc(harm)/field->at(i)->xks);      // here the harmonics have to be taken care
      rpart.push_back(0);
      rharm.push_back(static_cast<double>(harm));
    }
  }  


  xku=und->getku();
  if (xku==0){   // in the case of drifts - the beam stays in phase if it has the reference energy // this requires that the phase slippage is not applied
    xku=xks*0.5/und->getGammaRef()/und->getGammaRef();
  }
	    
  double aw=und->getaw();

  double autophase=und->autophase();


  // Runge Kutta solver to advance particle

  

  for (int is=0; is<beam->beam.size(); is++){    


      for (int ip=0; ip<beam->beam.at(is).size();ip++){
        gamma=beam->beam.at(is).at(ip).gamma;
        theta=beam->beam.at(is).at(ip).theta+autophase; // add autophase here
        double x =beam->beam.at(is).at(ip).x;
        double y =beam->beam.at(is).at(ip).y;
        double px=beam->beam.at(is).at(ip).px;
        double py=beam->beam.at(is).at(ip).py;
	double awloc=und->faw(x,y);                 // get the transverse dependence of the undulator field
        btpar=1+px*px+py*py+aw*aw*awloc*awloc;	  

	ez=0;

	cpart=0;
	double wx,wy;
	int idx;
        for (int ifld=0;ifld<nfld.size();ifld++){

	  int islice=(is+field->at(nfld[ifld])->first) % field->at(nfld[ifld])->field.size(); 

	  if (field->at(nfld[ifld])->getLLGridpoint(x,y,&wx,&wy,&idx)){ // check whether particle is on grid
           cpart=field->at(nfld[ifld])->field[islice].at(idx)*wx*wy;
           idx++;
           cpart+=field->at(nfld[ifld])->field[islice].at(idx)*(1-wx)*wy;
           idx+=field->at(nfld[ifld])->ngrid-1;
           cpart+=field->at(nfld[ifld])->field[islice].at(idx)*wx*(1-wy);
           idx++;
           cpart+=field->at(nfld[ifld])->field[islice].at(idx)*(1-wx)*(1-wy);
           rpart[ifld]=rtmp[ifld]*awloc*conj(cpart);
	  }
	}
	this->RungeKutta(delz);

        beam->beam.at(is).at(ip).gamma=gamma;
        beam->beam.at(is).at(ip).theta=theta; 
      }
    
  }
  return;
}

void BeamSolver::RungeKutta(double delz)
{
  // Runge Kutta Solver 4th order - taken from pushp from the old Fortran source


  // first step
  k2gg=0;
  k2pp=0;

  this->ODE(gamma,theta);

  // second step
  double stpz=0.5*delz;

  gamma+=stpz*k2gg;
  theta+=stpz*k2pp;
  
  k3gg=k2gg;
  k3pp=k2pp;

  k2gg=0;
  k2pp=0;

  this->ODE(gamma,theta);

  // third step
  gamma+=stpz*(k2gg-k3gg);
  theta+=stpz*(k2pp-k3pp);

  k3gg/=6;
  k3pp/=6;

  k2gg*=-0.5;
  k2pp*=-0.5;

  this->ODE(gamma,theta);

  // fourth step
  stpz=delz;

  gamma+=stpz*k2gg;
  theta+=stpz*k2pp;

  k3gg-=k2gg;
  k3pp-=k2pp;

  k2gg*=2;
  k2pp*=2;

  this->ODE(gamma,theta);
  gamma+=stpz*(k3gg+k2gg/6.0);
  theta+=stpz*(k3pp+k2pp/6.0);

  return;
}


void BeamSolver::ODE(double tgam,double tthet)
{

  // differential equation for longitudinal motion
  double ztemp1=-2./xks;
  complex<double> ctmp=0;
  for (int i=0; i<rpart.size();i++){
    ctmp+=rpart[i]*complex<double> (cos(rharm[i]*tthet), -sin(rharm[i]*tthet));
  }
  double btper0=btpar+ztemp1*ctmp.real();   //perpendicular velocity
  double btpar0=sqrt(1.-btper0/(tgam*tgam));     //parallel velocity
  k2pp+=xks*(1.-1./btpar0)+xku;             //dtheta/dz
  k2gg+=ctmp.imag()/btpar0/tgam-ez;         //dgamma/dz

  return; 
}

