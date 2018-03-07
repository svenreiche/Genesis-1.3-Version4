#include "LoadBeam.h"

LoadBeam::LoadBeam()
{
  gamma=5800/0.511; gammaref="";
  delgam=0; delgamref="";
  current=0; currentref="";
  ex=0.3e-6; exref="";
  ey=0.3e-6; eyref="";
  betax=15; betaxref="";
  betay=15; betayref="";
  alphax=0; alphaxref="";
  alphay=0; alphayref="";
  xcen=0; xcenref="";
  ycen=0; ycenref="";
  pxcen=0; pxcenref="";
  pycen=0; pycenref="";
  bunch=0; bunchref="";
  bunchphase=0; bunchphaseref="";
  emod=0; emodref="";
  emodphase=0; emodphaseref="";
}

LoadBeam::~LoadBeam(){}

void LoadBeam::usage(){

  cout << "List of keywords for BEAM" << endl;
  cout << "&beam" << endl;
  cout << " double gamma = gammaref / reference" << endl;
  cout << " double delgam = 0 / reference" << endl;
  cout << " double current = 1000 / reference" << endl;
  cout << " double ex = 0.3e-6 / reference" << endl;
  cout << " double ey = 0.3e-6 / reference" << endl;
  cout << " double betax = 15 / reference / matched" << endl;
  cout << " double betay = 15 / reference / matched" << endl;
  cout << " double alphax = 0 / reference / matched" << endl;
  cout << " double alphay = 0 / reference / matched" << endl;
  cout << " double xcenter = 0 / reference" << endl;
  cout << " double ycenter = 0 / reference" << endl;
  cout << " double pxcenter = 0 / reference" << endl;
  cout << " double pycenter = 0 / reference" << endl;
  cout << " double bunch = 0 / reference" << endl;
  cout << " double bunchphase = 0 / reference" << endl;
  cout << " double emod = 0 / reference" << endl;
  cout << " double emodphase = 0 / reference" << endl;
  cout << "&end" << endl << endl;
  return;
}


bool LoadBeam::init(int rank, int size, map<string,string> *arg, Beam *beam, Setup *setup, Time *time, Profile *prof, Lattice *lat)
{


  if (beam->beam.size()>0){
    if (rank==0) {cout << "*** Error: Cannot generat beam, because beam is already defined" << endl; }
    return false;
  }


  double lambda=setup->getReferenceLength();   // reference length for theta
  double sample=static_cast<double>(time->getSampleRate());         // check slice length
  bool one4one=setup->getOne4One();            // check for one4one simulations
  bool shotnoise=setup->getShotNoise();
  int npart=setup->getNpart();
  int nbins=setup->getNbins();
  bool dotime=time->isTime();                  // check for time simulation

  gamma=setup->getReferenceEnergy();           // get default energy from setup input deck
  lat->getMatchedOptics(&betax,&alphax,&betay,&alphay);  // use matched value if calculated

 
  map<string,string>::iterator end=arg->end();

  if (arg->find("gamma")!=end    ){this->reference(arg->at("gamma"),&gamma,&gammaref); arg->erase(arg->find("gamma"));}
  if (arg->find("delgam")!=end   ){this->reference(arg->at("delgam"),&delgam,&delgamref); arg->erase(arg->find("delgam"));}
  if (arg->find("current")!=end  ){this->reference(arg->at("current"),&current,&currentref); arg->erase(arg->find("current"));}
  if (arg->find("ex")!=end       ){this->reference(arg->at("ex"),&ex,&exref); arg->erase(arg->find("ex"));}
  if (arg->find("ey")!=end       ){this->reference(arg->at("ey"),&ey,&eyref); arg->erase(arg->find("ey"));}
  if (arg->find("betax")!=end    ){this->reference(arg->at("betax"),&betax,&betaxref); arg->erase(arg->find("betax"));}
  if (arg->find("betay")!=end    ){this->reference(arg->at("betay"),&betay,&betayref); arg->erase(arg->find("betay"));}
  if (arg->find("alphax")!=end   ){this->reference(arg->at("alphax"),&alphax,&alphaxref); arg->erase(arg->find("alphax"));}
  if (arg->find("alphay")!=end   ){this->reference(arg->at("alphay"),&alphay,&alphayref); arg->erase(arg->find("alphay"));}
  if (arg->find("xcenter")!=end  ){this->reference(arg->at("xcenter"),&xcen,&xcenref); arg->erase(arg->find("xcenter"));}
  if (arg->find("ycenter")!=end  ){this->reference(arg->at("ycenter"),&ycen,&ycenref); arg->erase(arg->find("ycenter"));}
  if (arg->find("pxcenter")!=end ){this->reference(arg->at("pxcenter"),&pxcen,&pxcenref); arg->erase(arg->find("pxcenter"));}
  if (arg->find("pycenter")!=end ){this->reference(arg->at("pycenter"),&pycen,&pycenref); arg->erase(arg->find("pycenter"));}
  if (arg->find("bunch")!=end    ){this->reference(arg->at("bunch"),&bunch,&bunchref); arg->erase(arg->find("bunch"));}
  if (arg->find("bunchphase")!=end ){this->reference(arg->at("bunchphase"),&bunchphase,&bunchphaseref); arg->erase(arg->find("bunchphase"));}
  if (arg->find("emod")!=end     ){this->reference(arg->at("emod"),&emod,&emodref); arg->erase(arg->find("emod"));}
  if (arg->find("emodphase")!=end ){this->reference(arg->at("emodphase"),&emodphase,&emodphaseref); arg->erase(arg->find("emodphase"));}

  if (arg->size()!=0){
    if (rank==0){ cout << "*** Error: Unknown elements in &beam" << endl; this->usage();}
    return false;
  }



  if (rank==0){cout << "Generating input particle distribution..." << endl; }


  double theta0=4.*asin(1.);

  if (one4one){
    nbins=1;
    theta0*=sample;
  }  
  
  if ( (npart % nbins) != 0){
    if (rank==0) { cout << "*** Error: NPART is not a multiple of NBINS" << endl; } 
    return false;
  }
 

  vector<double> s;
  int nslice=time->getPosition(&s);

  beam->init(time->getNodeNSlice(),nbins,lambda,sample*lambda,s[0],one4one);
  beam->initSorting(rank,size,false,one4one);  // sorting routine is initialized, with default values to suppress field slippage but do sorting if one4one is enabled

  int nbeam=1024;
  Particle *beamslice = new Particle [nbeam]; 
  BeamSlice slice;


  QuietLoading ql;
  ShotNoise sn;
  

  // choice of distribution and seed needs to be done here
  
  if (one4one){
    int base[2]={setup->getSeed(),rank};
    ql.init(one4one,&base[0]);
  } else {
    int base[6] = { 0,1,2,3,4,5};
    ql.init(one4one,&base[0]);
  }

  
  sn.init(setup->getSeed(),rank);


  for (int j=0; j<time->getNodeNSlice(); j++){
    int i=j+time->getNodeOffset();
    slice.current=prof->value(s[i],current,currentref);
    double ne=slice.current*lambda*sample/ce;
    slice.gamma =prof->value(s[i],gamma,gammaref);
    slice.delgam=prof->value(s[i],delgam,delgamref);
    slice.ex    =prof->value(s[i],ex,exref);
    slice.ey    =prof->value(s[i],ey,eyref);
    slice.betax =prof->value(s[i],betax,betaxref);
    slice.betay =prof->value(s[i],betay,betayref);
    slice.alphax=prof->value(s[i],alphax,alphaxref);
    slice.alphay=prof->value(s[i],alphay,alphayref);
    slice.xcen  =prof->value(s[i],xcen,xcenref);
    slice.ycen  =prof->value(s[i],ycen,ycenref);
    slice.pxcen =prof->value(s[i],pxcen,pxcenref);
    slice.pycen =prof->value(s[i],pycen,pycenref);
    slice.bunch =prof->value(s[i],bunch,bunchref);
    slice.emod  =prof->value(s[i],emod,emodref);
    slice.bunchphase =prof->value(s[i],bunchphase,bunchphaseref);
    slice.emodphase  =prof->value(s[i],emodphase,emodphaseref);
    
    int npartloc=npart;
    if (one4one) { npartloc=static_cast<int>(round(ne)); }
 
    if (npartloc>nbeam){
        delete [] beamslice;
        nbeam=npartloc;
        beamslice = new Particle [nbeam];
    }

    ql.loadQuiet(beamslice, &slice, npartloc, nbins,theta0,i);
    if ((shotnoise)&&(!one4one)&&(dotime)){ sn.applyShotNoise(beamslice,npartloc,nbins,ne); }

    beam->beam[j].resize(npartloc);
    beam->current[j]=slice.current;
    for (int k=0;k<npartloc;k++){
      beam->beam[j].at(k).gamma=beamslice[k].gamma;
      beam->beam[j].at(k).theta=beamslice[k].theta;
      beam->beam[j].at(k).x    =beamslice[k].x;
      beam->beam[j].at(k).y    =beamslice[k].y;
      beam->beam[j].at(k).px   =beamslice[k].px;
      beam->beam[j].at(k).py   =beamslice[k].py;
      }
  }



  delete [] beamslice;  



 
  return true;

}

