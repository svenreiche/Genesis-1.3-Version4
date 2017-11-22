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


bool Lattice::parse(string filename, string beamline, int rank)
{

  // release old lattice

  for (int i=0;i<lat.size();i++){
    delete lat[i];
  }
  lat.clear();

  LatticeParser parser;
  matched=false;

  if (rank == 0) { cout << "Parsing lattice file..." << endl; }
  bool err=parser.parse(filename,beamline,rank, lat);
  if (err==false) { return err; }
  
  layout.clear();
  
  for(int i=0; i<lat.size();i++){
    double z0=lat[i]->z;    
    if (lat[i]->type.compare("Drift")==0){
      continue;
    }
    layout[z0]=1;
    layout[z0+lat[i]->l]=1;
  }
  int last = lat.size()-1;
  if (lat[last]->type.compare("Drift")==0){
    double z0=lat[last]->z;
    layout[z0]=1;
    layout[z0+lat[last]->l]=1;
  }
  return true;
}



bool Lattice::writeLattice(hid_t fid, double delz, double lambda, double gamma, AlterLattice *alt)
{
  
  this->unrollLattice(delz);  
  this->calcSlippage(lambda,gamma);
  
  string group="Lattice";
  hid_t gid=H5Gcreate(fid,group.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  int ndata=lat_aw.size();
    
  this->writeDataDouble(gid,(char *)"z",&lat_z[0],ndata);
  this->writeDataDouble(gid,(char *)"dz",&lat_dz[0],ndata);
  this->writeDataDouble(gid,(char *)"aw",&lat_aw[0],ndata);
  this->writeDataDouble(gid,(char *)"ax",&lat_ax[0],ndata);
  this->writeDataDouble(gid,(char *)"ay",&lat_ay[0],ndata);
  this->writeDataDouble(gid,(char *)"ku",&lat_ku[0],ndata);
  this->writeDataDouble(gid,(char *)"kx",&lat_kx[0],ndata);
  this->writeDataDouble(gid,(char *)"ky",&lat_ky[0],ndata);
  this->writeDataDouble(gid,(char *)"gradx",&lat_gradx[0],ndata);
  this->writeDataDouble(gid,(char *)"grady",&lat_grady[0],ndata);
  this->writeDataDouble(gid,(char *)"qf",&lat_qf[0],ndata);
  this->writeDataDouble(gid,(char *)"qx",&lat_qx[0],ndata);
  this->writeDataDouble(gid,(char *)"qy",&lat_qy[0],ndata);
  this->writeDataDouble(gid,(char *)"delay",&lat_delay[0],ndata);
  this->writeDataDouble(gid,(char *)"lb",&lat_lb[0],ndata);
  this->writeDataDouble(gid,(char *)"ld",&lat_ld[0],ndata);
  this->writeDataDouble(gid,(char *)"lt",&lat_lt[0],ndata);
  this->writeDataDouble(gid,(char *)"cx",&lat_cx[0],ndata);
  this->writeDataDouble(gid,(char *)"cy",&lat_cy[0],ndata);
  this->writeDataInt(gid,(char *)"marker",&lat_mk[0],ndata);
  this->writeDataInt(gid,(char *)"helical",&lat_helical[0],ndata);
  this->writeDataDouble(gid,(char *)"slippage",&lat_slip[0],ndata);
  this->writeDataDouble(gid,(char *)"phaseshift",&lat_phase[0],ndata);
  

  H5Gclose(gid);

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

      if ((Lz>0) && (i>0)){ // apply last drift section
        tmp=Lz/2/gamma/gamma/lambda;  
        lat_slip[i-1]+=floor(tmp);   
        tmp-=floor(tmp);
        lat_phase[i-1]+=tmp*4*asin(1);   // auto phasing
	Lz=0; // clear last intra beam section 
      } 
      // reset drift parameters
    } else {
      Lz   +=lat_dz[i];
      tmp=lat_delay[i]/lambda;  // affect of chicane is always autophasing!!!
      lat_slip[i]=floor(tmp);
      lat_phase[i]=lat_ps[i];     // phase shifter goes here
    }
  }
  // correct for the end of the lattice that autophasing is applied in case for second, suceeding run
  tmp=Lz/2/gamma/gamma/lambda;  
  tmp-=floor(tmp);
  lat_phase[nz-1]+=tmp*4*asin(1);

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

    if (iele!=-1){                                     // outside of an undulator
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



























