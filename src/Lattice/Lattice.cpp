#include "Lattice.h"

Lattice::Lattice()
{
  matched=false;
}


Lattice::~Lattice()
{
  for (int i=0;i<lat.size();i++){
    delete lat[i];
  }
  lat.clear();
}


bool Lattice::parse(string filename, string beamline, int rank, bool streaming)
{

  // release old lattice

  for (int i=0;i<lat.size();i++){
    delete lat[i];
  }
  lat.clear();

  LatticeParser parser;
  matched=false;

  if (rank == 0) { cout << "Parsing lattice file..." << endl; }
  bool err=parser.parse(filename,beamline,rank, lat,streaming);
  if (err==false) { 
    return err; 
  }
  
  layout.clear();  // layout holds start and end point of each value
  
  for(int i=0; i<lat.size();i++){
    double z0=lat[i]->z;    
    layout[z0]=1;
    layout[z0+lat[i]->l]=1;
  }

  return true;
}


// generate from the input lattice the explicit representation, including some modification with 'alt'

bool Lattice::generateLattice(double delz, double lambda, double gamma, AlterLattice *alt,Undulator *und)
{

  this->unrollLattice(delz);  
  this->calcSlippage(lambda,gamma);

  und->setGammaRef(gamma);

  int ndata=lat_aw.size();

  und->aw.resize(ndata);
  und->ax.resize(ndata);
  und->ay.resize(ndata);
  und->ku.resize(ndata);
  und->kx.resize(ndata);
  und->ky.resize(ndata);
  und->gradx.resize(ndata);
  und->grady.resize(ndata);
  und->qf.resize(ndata);
  und->qx.resize(ndata);
  und->qy.resize(ndata);
  und->cx.resize(ndata);
  und->cy.resize(ndata);
  und->chic_angle.resize(ndata);
  und->chic_lb.resize(ndata);
  und->chic_ld.resize(ndata);
  und->chic_lt.resize(ndata);
  und->slip.resize(ndata);
  und->phaseshift.resize(ndata);
  und->z.resize(ndata);
  und->dz.resize(ndata);
  und->helical.resize(ndata);
  und->marker.resize(ndata+1);

  for (int i=0; i<ndata;i++){
      und->aw[i]=lat_aw[i];
      und->ax[i]=lat_ax[i];
      und->ay[i]=lat_ay[i];
      und->ku[i]=lat_ku[i];
      und->kx[i]=lat_kx[i];
      und->ky[i]=lat_ky[i];
      und->gradx[i]=lat_gradx[i];
      und->grady[i]=lat_grady[i];
      und->qf[i]=lat_qf[i]; 
      und->qx[i]=lat_qx[i];
      und->qy[i]=lat_qx[i];
      und->cx[i]=lat_cx[i];
      und->cy[i]=lat_cy[i];

      und->chic_angle[i]=lat_delay[i];  // here it is the delay but will converted to angle
      und->chic_lb[i]=lat_lb[i];
      und->chic_ld[i]=lat_ld[i];
      und->chic_lt[i]=lat_lt[i];       // here it is the total length but it will change to the time delay.
      
      if (und->chic_angle[i]!=0){
       double delay=fabs(und->chic_angle[i]);
       double tmin=0;
       double tmax=asin(1)-0.001;
       bool converged=false;
       double theta,d;
       while (!converged){
         theta=0.5*(tmax+tmin);
         d=4*und->chic_lb[i]*(theta/sin(theta)-1)+2*und->chic_ld[i]*(1/cos(theta)-1);
         if (d>delay) {
           tmax=theta;
         } else {
          tmin=theta;
	 }
         if (fabs(delay-d)<1e-15) { converged=true; }
       }
       und->chic_angle[i]=theta;
      }


      und->slip[i]=lat_slip[i];
      und->phaseshift[i]=lat_phase[i];
      und->z[i]=lat_z[i];
      und->dz[i]=lat_dz[i];
      und->helical[i]=lat_helical[i];
      und->marker[i]=lat_mk[i];

  }
  int lastmark=this->findMarker(und->z[ndata-1]+und->dz[ndata-1],"Marker");
  if (lastmark < 0){
    und->marker[ndata]=0;
  } else {
        Marker *mark=(Marker *)lat[lastmark];
        und->marker[ndata]=mark->action; 
  }


  return true;




}



void Lattice::calcSlippage(double lambda, double gamma)
{

  int nz=lat_aw.size();
  lat_slip.resize(nz);                                  // this needs to be improved for chicanes
  lat_phase.resize(nz);
  
  // calc the path length for a chicane

  double Lz=0;    // projected path
  double tmp;

  for (int i=0; i< nz;i++){
    if (lat_aw[i]>0){ // within undulator
      tmp=2*gamma*gamma*lambda/(1+lat_aw[i]*lat_aw[i]);
      lat_slip[i]=lat_dz[i]/tmp;
      lat_phase[i]=0;

      if ((Lz>0) && (i>0)){ // apply previous drift section
        tmp=Lz/2/gamma/gamma/lambda;  
        lat_slip[i-1]+=floor(tmp)+1;  //auto phasing would always add some slippage   
	//        tmp-=floor(tmp);
	//        lat_phase[i-1]+=tmp*4*asin(1);   // auto phasing
	Lz=0; // clear last intra beam section 
      } 
      // reset drift parameters
    } else {
      Lz   +=lat_dz[i];
      //      if (lat_delay[i]>0) { cout << "Delay: " << lat_delay[i] << " Lambda: " << lambda << endl; } 
      tmp=lat_delay[i]/lambda;  // affect of chicane is always autophasing!!!
      lat_slip[i]=floor(tmp);
      lat_phase[i]=lat_ps[i];     // phase shifter goes here
    }
  }
  // correct for the end of the lattice that autophasing is applied in case for second, suceeding run
  tmp=Lz/2/gamma/gamma/lambda;  
  lat_slip[nz-1]+=floor(tmp)+1;   
  //  tmp-=floor(tmp);
  // lat_phase[nz-1]+=tmp*4*asin(1);

}



int Lattice::findElement(double z0, double z1,string type)
{
  double zmid=0.5*(z0+z1);
  for(int i=0;i<lat.size();i++){
    double zz0 =lat[i]->z;
    double zz1 =zz0+lat[i]->l;
    if ((zmid>zz0)&&(zmid<zz1)&&(lat[i]->type.compare(type)==0)){
      return i;
    }
  }
  return -1; 
}

int Lattice::findMarker(double z0,string type)
{

  for(int i=0;i<lat.size();i++){
    double zz0 =lat[i]->z-z0;
    if ((zz0*zz0<1e-6) && (lat[i]->type.compare(type)==0)){  
      return i;
    }
  }
  return -1; 
}




void Lattice::getMatchedOptics(double *betax, double *alphax, double *betay, double *alphay)
{

  if (matched){
    *betax =mbetax;
    *alphax=malphax;
    *betay =mbetay;
    *alphay=malphay;
  }

  return;

}

void Lattice::match(int rank, double z0, double gammaref)
{

  double phix,phiy;

  this->unrollLattice(20); 
  Optics opt;
  opt.init();

  int i=0;

  while((lat_z[i]<z0)&&(i<lat_z.size())){
 
    double qf=lat_qf[i];
    double qx=lat_aw[i]*lat_aw[i]/gammaref/gammaref;
    double qy=qx*lat_ky[i];
    qx*=lat_kx[i];
    double dz=lat_dz[i];
    if ((lat_z[i]+lat_dz[i])>z0){ dz=z0-lat_z[i]; }
    opt.addElement(dz,qf,qx,qy);
    i++; 
  }

  if ((i==lat_z.size())&&((lat_z[i-1]+lat_dz[i-1])<z0)){
      if (rank==0){ cout << "*** Warning: Matching position beyond end of lattice. Ignoring matching command" << endl; }
      matched=false;
      return; 
  }


  if (matched=opt.match(&mbetax,&malphax,&mbetay,&malphay,&phix,&phiy)){
    if (rank==0){
        cout << "Matching for periodic solution between z = 0 and z = "<<z0 << " :" << endl;
        cout << "   betax (m) : " << mbetax  << endl;
        cout << "   alphax    : " << malphax << endl;
        cout << "   phix (deg): " << phix    << endl; 
        cout << "   betay (m) : " << mbetay  << endl;
        cout << "   alphay    : " << malphay << endl;
        cout << "   phiy (deg): " << phiy    << endl; 
    }
    matched=true;
  } else {
    if (rank==0){ cout << "*** Warning: Cannot find a matching solution" << endl; }
    matched=false;
  }
  return;
}
  

void Lattice::unrollLattice(double delz)
{

  lat_aw.clear();
  lat_ku.clear();
  lat_kx.clear();
  lat_ky.clear();
  lat_gradx.clear();
  lat_grady.clear();
  lat_ax.clear();
  lat_ay.clear();
  lat_helical.clear();
  lat_dz.clear();
  lat_z.clear();



  map<double,int>::iterator it=layout.begin();
  map<double,int>::iterator iend=layout.end();
  iend--;
  

  // first fill lattice with undulator
  

  while(it !=iend){
      //default

      double z0=it->first;
      it++;
      double z1=it->first;

      double dz=z1-z0;

      int iele=this->findElement(z0,z1,"Undulator");
      if (iele==-1){                                    // outside of an undulator
        lat_dz.push_back(z1-z0);
        lat_aw.push_back(0);
        lat_ku.push_back(0);
        lat_kx.push_back(0);
        lat_ky.push_back(0);
        lat_gradx.push_back(0);
        lat_grady.push_back(0);
        lat_ax.push_back(0);
        lat_ay.push_back(0);
        lat_helical.push_back(0);

      } else {
        ID *und=(ID *)lat[iele];
        int nz=round(dz/delz);
        if (nz==0) {nz=1;}
        dz=dz/static_cast<double>(nz);
        for (int iz=0;iz<nz;iz++){
	  lat_dz.push_back(dz);
          lat_aw.push_back(und->aw);
	  double ku=4.*asin(1)/und->lambdau;
          lat_ku.push_back(ku);
          lat_kx.push_back(ku*ku*und->kx);
          lat_ky.push_back(ku*ku*und->ky);
          lat_gradx.push_back(ku*und->gradx);
          lat_grady.push_back(ku*und->grady);
          lat_ax.push_back(und->ax);
          lat_ay.push_back(und->ay);
          lat_helical.push_back(static_cast<int>(und->helical));
        }
      }
  }

  int nz=lat_aw.size();
  lat_z.resize(nz);
  lat_z[0]=0; 
  for (int i=1; i<nz;i++){
    lat_z[i]=(lat_z[i-1]+lat_dz[i-1]);
  }   

  // resize vector and assign default zero values 
  lat_qf.resize(nz);
  lat_qx.resize(nz);
  lat_qy.resize(nz);
  lat_delay.resize(nz);
  lat_lb.resize(nz);
  lat_ld.resize(nz);
  lat_lt.resize(nz);
  lat_cx.resize(nz);
  lat_cy.resize(nz);
  lat_mk.resize(nz);
  lat_ps.resize(nz);
  for (int i=0; i<nz;i++){
    lat_qf[i]=0;
    lat_qx[i]=0;
    lat_qy[i]=0;
    lat_delay[i]=0;
    lat_lb[i]=0;
    lat_ld[i]=0;
    lat_lt[i]=0;
    lat_cx[i]=0;
    lat_cy[i]=0;
    lat_mk[i]=0;
    lat_ps[i]=0;
  }

  int lastChicane=-1;


  for (int i=0;i<nz;i++){
    double z0=lat_z[i];
    double z1=z0+lat_dz[i];

    bool inUnd=(lat_aw[i]>0);
    bool inQuad=0;
    int iele=this->findElement(z0,z1,"Quadrupole");

    if (iele!=-1){     // found quadrupole
        Quadrupole *quad=(Quadrupole *)lat[iele];
        lat_qf[i]=quad->k1; 
        lat_qx[i]=quad->dx; 
        lat_qy[i]=quad->dy; 
	inQuad=true;
     }

    iele=this->findElement(z0,z1,"Chicane");
    if (iele!=-1){ 
      if (iele!=lastChicane){
	if ((!inUnd)&&(!inQuad)){
          Chicane *chicane=(Chicane *)lat[iele];
          lat_delay[i]=chicane->delay; 
          lat_lb[i]=chicane->lb; 
          lat_ld[i]=chicane->ld; 
          lat_lt[i]=chicane->l;    // save length of chicane in case that the chicane is split in parts
	  lastChicane=iele;
	} else {
	  cout << "*** Warning: Chicane inside Undulator or Quadrupole. Ignoring Chicane" << endl; 
	}
      }
    } 

    iele=this->findElement(z0,z1,"Corrector");
    if (iele!=-1){                                     // outside of a bend
        Corrector *cor=(Corrector *)lat[iele];
        lat_cx[i]=cor->cx; 
        lat_cy[i]=cor->cy; 
    } 

    iele=this->findElement(z0,z1,"Phaseshifter");
    if (iele!=-1){                                     // outside of a bend
        Phaseshifter *cor=(Phaseshifter *)lat[iele];
        lat_ps[i]=cor->phi; 
    } 

    // note that markers are only find at the beginning of the integration step. Thus a marker as the last element is ignored.
    iele=this->findMarker(z0,"Marker");
    if (iele!=-1){   
        Marker *mark=(Marker *)lat[iele];
        lat_mk[i]=mark->action; 
    } 
  }

  return;  

}



























