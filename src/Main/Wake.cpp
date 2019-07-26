#include "Wake.h"
#include "Beam.h"

Wake::Wake()
{
  radius=2.5e-3;
  conductivity=0;
  relaxation=0;
  roundpipe=true;
  transient=false;
  ztrans=0;
  gap=0;
  lgap=1;
  hrough=0;
  lrough=1;
  ns=0;
}

Wake::~Wake(){}


void Wake::usage(){

  cout << "List of keywords for Wake" << endl;
  cout << "&wake" << endl;
  cout << " double radius = 2.5e-3" << endl;
  cout << " bool   roundpipe   = true" << endl;
  cout << " string material  = <empty>" << endl;
  cout << " double conductivity = 0e-6" << endl;
  cout << " double relaxation  = 0e-6" << endl;
  cout << " double gap = 0e-6" << endl;
  cout << " double lgap  = 1.0" << endl;
  cout << " double hrough = 0.0" << endl;
  cout << " double lrough = 1.0" << endl;
  cout << " bool transient = false" << endl;
  cout << " double ztrans = 0" << endl;
  cout << "&end" << endl << endl;
  return;
}

// input parameter

bool Wake::init(int rank, int size, map<string,string> *arg,  Time *time, Setup *setup, Beam *beam)
{

  string material="";
  map<string,string>::iterator end=arg->end();
  map<string,string>::iterator iter=arg->begin();
  

  if (arg->find("radius")!=end) {radius = atof(arg->at("radius").c_str());  arg->erase(arg->find("radius"));}
  if (arg->find("conductivity")!=end) {conductivity= atof(arg->at("conductivity").c_str());  arg->erase(arg->find("conductivity"));}
  if (arg->find("relaxation")!=end) {relaxation = atof(arg->at("relaxation").c_str());  arg->erase(arg->find("relaxation"));}
  if (arg->find("roundpipe")!=end)    {roundpipe    = atob(arg->at("roundpipe").c_str());  arg->erase(arg->find("roundpipe"));}
  if (arg->find("material")!=end){material = arg->at("material"); arg->erase(arg->find("material"));}
  if (arg->find("gap")!=end)  {gap = atof(arg->at("gap").c_str());  arg->erase(arg->find("gap"));}
  if (arg->find("lgap")!=end) {lgap = atof(arg->at("lgap").c_str());  arg->erase(arg->find("lgap"));}
  if (arg->find("hrough")!=end) {hrough = atof(arg->at("hrough").c_str());  arg->erase(arg->find("hrough"));}
  if (arg->find("lrough")!=end) {lrough = atof(arg->at("lrough").c_str());  arg->erase(arg->find("lrough"));}
  if (arg->find("transient")!=end)    {transient    = atob(arg->at("transient").c_str());  arg->erase(arg->find("transient"));}
  if (arg->find("ztrans")!=end) {ztrans = atof(arg->at("ztrans").c_str());  arg->erase(arg->find("ztrans"));}


  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &wake" << endl; this->usage();}
    return false;
  }
  
  if ((material=="CU") || (material=="Cu") || (material=="cu")){
    conductivity=5.813e7;
    relaxation=8.1e-6;
  }

  if ((material=="AL") || (material=="Al") || (material=="al")){
    conductivity=3.571e7;
    relaxation=2.4e-6;
  }

  // getting the time parameter from the time window definition
  unsigned int nslen=time->getNodeNSlice()*size;
  unsigned int nsample=time->getSampleRate();

  slen=time->getTimeWindowLength();
  ds=slen/static_cast<double>(nslen*nsample);

  // clear single particle wakes
  ns=nslen*nsample;    // sample every wavelength
  wakeres = new double [ns];
  wakegeo = new double [ns];
  wakerou = new double [ns];

  // calculate the single particle wakes
  this->singleWakeResistive(rank);
  this->singleWakeGeometric(rank);
  this->singleWakeRoughness(rank);


  // transfer wakes into beam class
  beam->initWake(ns,ds,wakeres,wakegeo,wakerou, ztrans,radius,transient);

  delete[] wakeres,wakegeo,wakerou;

  return true;

}



void Wake::KernelRoughness(vector< complex<double> > *K, complex<double> q1, complex<double> q2)
{
  int N=K->size();
  complex<double> dq=(q2-q1)/static_cast<double>(N-1);
  complex<double> i=complex<double> (0,1);

  for (int j=0; j< N;j++){
    complex<double> q=q1+static_cast<double>(j)*dq;
    complex<double> S=(sqrt(2.*q+1.)-i*sqrt(2.*q-1.))*q/sqrt(4.*q*q-1.);
    K->at(j)=(S+1.)/(1.-i*rrough*q*S)/(1.+i*rrough*q);
  }
  K->at(0)*=0.5;
  K->at(N-1)*=0.5;   // for trapazoidal integration
}



double Wake::TrapIntegrateRoughness(vector< complex<double> > *K, complex<double> q1, complex<double> q2, double tau)
{
  int N=K->size();
  complex<double> dq=(q2-q1)/static_cast<double>(N-1);
  complex<double> i=complex<double> (0,1);  
  complex<double> val=0;

  for (int j=0; j< N;j++){
    complex<double> q=q1+static_cast<double>(j)*dq;
    val+=exp(-i*q*tau)*K->at(j);
  }
  val*=dq;
  return val.real();
}





void Wake::singleWakeRoughness(int rank)
{

  if (hrough <=0){
    for (int i=0; i< ns; i++){
      wakerou[i]=0;
    }
    return;
  }


  double pi=2.*asin(1);
  rrough = pi*pi*pi/lrough/lrough/lrough*hrough*hrough*radius; // aspect ratio parameter. If much smaller than one an inductive mode
             // is expected, if unity or larger it is rather a synchronous mode.

  double tau=0;  // example 
  complex<double> q1 = complex<double> (0,0);
  complex<double> q2 = complex<double> (0,2e-3);
  complex<double> q3 = complex<double> (1,2e-3);
  complex<double> q4 = complex<double> (1,0);
  complex<double> q5 = complex<double> (100,0);

  int N=128;
  vector< complex<double> > K1,K2,K3,K4;
  K1.resize(N);
  K2.resize(8*N);
  K3.resize(N);
  K4.resize(8*N);

  this->KernelRoughness(&K1,q1,q2);
  this->KernelRoughness(&K2,q2,q3);
  this->KernelRoughness(&K3,q3,q4);
  this->KernelRoughness(&K4,q4,q5);


  double res;
  double coef=rrough/pi*4/radius/radius*1.6e-19/4/pi/8.854e-12;   

  for (int i=0; i<ns; i++){
    tau=2*pi*ds*i/lrough;
    res =this->TrapIntegrateRoughness(&K1,q1,q2,tau);
    res+=this->TrapIntegrateRoughness(&K2,q2,q3,tau);
    res+=this->TrapIntegrateRoughness(&K3,q3,q4,tau);
    res+=this->TrapIntegrateRoughness(&K4,q4,q5,tau);
    wakerou[i]=coef*res;
  }
  
  if (rank==0){
       cout << "Roughness Wake calculated..." << endl;   
  }

 
}


void Wake::singleWakeGeometric(int rank)
{

  if (gap <=0){
    for (int i=0; i< ns; i++){
      wakegeo[i]=0;
    }
    return;
  }

  double pi=2.*asin(1.);
  double coef=-vacimp*ce/(pi*pi*radius*lgap)*2*sqrt(0.5*gap); // scaling coefficient
  if (!roundpipe) { coef*=0.956; }     //
  for (int i = 0;i<ns;i++){
    wakegeo[i]=coef*sqrt(ds*i);     
  }

  if (rank==0){
       cout << "Geometric Wake calculated..." << endl;   
  }

}




 
void Wake::singleWakeResistive(int rank)
{

 if (conductivity <=0){
    for (int i=0; i< ns; i++){
      wakeres[i]=0;
    }
    return;
  }
   double s0=pow(2*radius*radius/vacimp/conductivity,1./3.); // characteristic length in SI units
   double gamma=relaxation/s0;
   double coef = radius/(s0*s0);
   double pi=2.*asin(1.);

   double kappamax=100;  // empirical cut-off in impedance spectrum

   unsigned int nk=1000;
   unsigned int nq=10000;

   double Zre[1000],Zim[1000]; 
   if (roundpipe) {
     for (int i =0; i<nk; i++){
        double kappa=(i+1.)*kappamax/nk;  // value of kappa	 
 	double t = kappa/sqrt(1+kappa*kappa*gamma*gamma);
        double lambdaRe=coef*sqrt(t)*sqrt(1.-t*gamma);
  	double lambdaIm=coef*sqrt(t)*sqrt(1.+t*gamma)-kappa*kappa*radius*0.5/s0/s0;
	double nomi=2.*kappa/(3e8*radius*s0)/(lambdaRe*lambdaRe+lambdaIm*lambdaIm);
	Zre[i]=lambdaRe*nomi;
	Zim[i]=-lambdaIm*nomi;
     }
   } else {
      double coh[10000],sih[10000];
      double dq=15/static_cast<double>(nq-1);
      for (int i =1; i<10000;i++){
	   coh[i]=0.5*(exp(dq*i)+exp(-dq*i));
	   sih[i]=coh[i]-exp(-dq*i);
	   sih[i]/=dq*i;
      }
      coh[0]=1.;
      sih[0]=1.;
        
      for (int i =0; i<nk; i++){       
           double kappa=(i+1.)*kappamax/nk;  // value of kappa	 
 	   double t = kappa/sqrt(1+kappa*kappa*gamma*gamma);
           double scale=2.*15.*kappa/(3e8*radius*s0*(2*nq-1));
	   Zre[i]=0;
	   Zim[i]=0;	 
       // integrate over q-> infty which is actually exp(30)
           for (int j=1;j<nq;j++){
               double lambdaRe=coef*sqrt(t)*sqrt(1.-t*gamma)*coh[j]*coh[j];
    	       double lambdaIm=coef*sqrt(t)*sqrt(1.+t*gamma)*coh[j]*coh[j]-kappa*kappa*radius*0.5/s0/s0*sih[j]*coh[j];
	       double nomi=scale/(lambdaRe*lambdaRe+lambdaIm*lambdaIm);
	       Zre[i]+=lambdaRe*nomi;
	       Zim[i]+=-lambdaIm*nomi;	 
	   }
      }
   }	 	 
   
 
   // reconstructing the wakepotential
   coef=ds*kappamax/nk/s0;
   for (int i=0;i<ns;i++){
      wakeres[i]=0;
      for (int j=0;j<nk;j++){
	  double phi=i*coef*(1.e0+j);
	  wakeres[i]+=Zre[j]*cos(phi)+Zim[j]*sin(phi);
      }
   }
   coef=-kappamax/nk/s0*3e8/pi*(vacimp*ce/4/pi); // to scale to SI units
   for (int i = 0; i < ns; i++){
     wakeres[i]*=coef;
   }


  /*
   ofstream myfile;

   myfile.open ("Zre.txt");
   for (int i = 0; i < nk ; i++){
	myfile << Zre[i] << endl;
   }
   myfile.close();   

   myfile.open ("Zim.txt");
   for (int i = 0; i < nk ; i++){
	myfile << Zim[i] << endl;
   }
   myfile.close(); 

   myfile.open ("wake.txt");
   for (int i = 0; i < ns ; i++){
	myfile << wakeres[i] << endl;
   }
   myfile.close();
   */



   if (rank==0){
       cout << "Resistive Wake calculated (s0 = " << s0 << ")..." << endl;   
   }

}

