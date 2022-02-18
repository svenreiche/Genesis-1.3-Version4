#include "FieldSolver.h"
#include "Field.h"
#include "Beam.h"

FieldSolver::FieldSolver()
{
  delz_save=0;
  difffilter_ = false;
  ngrid = 0;
}

FieldSolver::~FieldSolver(){
#ifdef FFTW
    if (ngrid > 0){
        delete[] in;
        delete[] out;
        fftw_destroy_plan(p);
        fftw_destroy_plan(pi);
    }
#endif
}


void FieldSolver::init(int ngrid_in){
    /**
     * Initialize the field solver to allocate memory etc.
     * @param ngrid_in - number of grid points in one dimension
     * @return none
     */

    if (ngrid !=ngrid_in) {
        ngrid = ngrid_in;
        c.resize(ngrid);
        r.resize(ngrid * ngrid);
        cbet.resize(ngrid);
        cwet.resize(ngrid);
        crsource.resize(ngrid * ngrid);
#ifdef FFTW
        if (ngrid > 0) {
            delete[] in;
            delete[] out;
            fftw_destroy_plan(p);
            fftw_destroy_plan(pi);
        }
        sigmoidx_.resize(ngrid);
        sigmoidy_.resize(ngrid);
        in = new complex<double>[ngrid * ngrid];
        out = new complex<double>[ngrid * ngrid];
        p = fftw_plan_dft_2d(ngrid, ngrid, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out),
                             FFTW_FORWARD, FFTW_MEASURE);
        pi = fftw_plan_dft_2d(ngrid, ngrid, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out),
                              FFTW_BACKWARD, FFTW_MEASURE);
#endif
    }
}

void FieldSolver::initSourceFilter(bool do_filter,double kx0, double ky0, double ksig) {
    /**
    * Initializes all parameter for filtering the source term
    * @param do_filter - filter to enable and disable filtering
    * @param kx0 - the relative spatial frequency in x, where sigmoid function is 1/2
    * @param ky0 - the relative spatial frequency in y, where sigmoid function is 1/2
    * @param ksig - the steepness of the sigmoid function
    * @returns None
    */
#ifdef FFTW
    difffilter_ = do_filter;
    if (ksig <= 0) {  // check for unphysical input
        difffilter_ = false;
    }
    if (difffilter_){  // pre-calculation of sigmoid function
        int nhalf = (ngrid-1)/2;
        double xc = static_cast<double>(nhalf)*kx0;
        double yc = static_cast<double>(nhalf)*ky0;
        double sig = static_cast<double>(nhalf)*ksig;
        for (int i=0; i<ngrid;i++){
            int j = ((i+nhalf) % ngrid) - ngrid; // sequence: 0 1 2 3 ... n -(n-1) -(n-2) .... -2 -1  -> FFT order
            double arg = (fabs(static_cast<double>(j))-xc)/sig;
            sigmoidx_[i] = 1./(1.+exp(arg));
            arg = (fabs(static_cast<double>(j))-yc)/sig;
            sigmoidy_[i] = 1./(1.+exp(arg));
        }
    }
#endif
}


void FieldSolver::advance(double delz, Field *field, Beam *beam, Undulator *und)
{
 
  for (int ii=0; ii<field->field.size();ii++){  // ii is index for the beam
    int i= (ii+field->first) % field->field.size();           // index for the field
  
    // clear source term
    for (int i=0; i< ngrid*ngrid; i++){
      crsource[i]=0;
    }

    int harm=field->getHarm();
    //   cout << "Harm: " << harm << " Coupling: " << und->fc(harm) << endl;
    // construc source term
    if (und->inUndulator()&& field->isEnabled()) { 
       double scl=und->fc(harm)*vacimp*beam->current[ii]*field->xks*delz;
       scl/=4*eev*beam->beam[ii].size()*field->dgrid*field->dgrid;
       complex<double> cpart;
       double part,weight,wx,wy;
       int idx;

       for (int ip=0;ip<beam->beam.at(ii).size();ip++){
	      double x    =beam->beam.at(ii).at(ip).x;
	      double y    =beam->beam.at(ii).at(ip).y;
	      double theta=static_cast<double>(harm)*beam->beam.at(ii).at(ip).theta;
	      double gamma=beam->beam.at(ii).at(ip).gamma;
          if (field->getLLGridpoint(x,y,&wx,&wy,&idx)){
            part=sqrt(und->faw2(x,y))*scl/gamma;
           // tmp  should be also normalized with beta parallel
            cpart=complex<double>(sin(theta),cos(theta))*part;
            weight=wx*wy;
 	        crsource[idx]+=weight*cpart;
            weight=(1-wx)*wy;
            idx++;
	        crsource[idx]+=weight*cpart;
            weight=wx*(1-wy);
            idx+=ngrid-1;
            crsource[idx]+=weight*cpart;
            weight=(1-wx)*(1-wy);
            idx++;
	        crsource[idx]+=weight*cpart;

          }
       } 
    }  // end of source term construction
    this->filterSourceTerm();
    this->ADI(field->field[i]);
  }
  return;
}

void FieldSolver::filterSourceTerm()
{
    /**
    * filter the source term to avoid emission under strong angle. The source term is transformed via
    * FFTW2D. Spectral components under a large angle are damped according to a sigmoid function.
    * @returns None
    */
    if (!difffilter_) { return; }

#ifdef FFTW
    for (int idx=0; idx <ngrid*ngrid;idx++){
        in[idx]=crsource[idx];   // field for the FFT
    }
    fftw_execute(p);
    // the sigmoid function can be pre-calculated
    int idx = 0;
    for (int ix=0; ix <ngrid; ix++){
        for (int iy; iy < ngrid; iy++){
            in[idx]=out[idx]*sigmoidx_[ix]*sigmoidy_[iy];
            idx++;
        }
    }
    fftw_execute(pi);
    double norm = static_cast<double>(ngrid*ngrid*ngrid*ngrid);
    for (int idx=0; idx <ngrid*ngrid;idx++){
        crsource[idx]=out[idx]*norm;
    }
#endif
}



void FieldSolver::ADI(vector<complex<double> > &crfield)
{
  int ix,idx;
  // implicit direction in x
  for (idx=0;idx<ngrid;idx++){
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+ngrid]-2.0*crfield[idx]);
  }
  for (idx=ngrid;idx<ngrid*(ngrid-1);idx++){
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+ngrid]-2.0*crfield[idx]+crfield[idx-ngrid]);
  }
  for (idx=ngrid*(ngrid-1);idx<ngrid*ngrid;idx++){
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx-ngrid]-2.0*crfield[idx]);
  }

  // solve tridiagonal system in x
  this->tridagx(crfield);

  // implicit direction in y
  for(ix=0;ix<ngrid*ngrid;ix+=ngrid){
    idx=ix;
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+1]-2.0*crfield[idx]);
    for(idx=ix+1;idx<ix+ngrid-1;idx++){
      r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx+1]-2.0*crfield[idx]+crfield[idx-1]);
    }
    idx=ix+ngrid-1;
    r[idx]=crsource[idx]+crfield[idx]+cstep*(crfield[idx-1]-2.0*crfield[idx]);
  }

  // solve tridiagonal system in y
  this->tridagy(crfield);


  return;
}


void FieldSolver::tridagx(vector<complex<double > > &u)
{ 
  for (int i=0;i<ngrid*(ngrid-1);i+=ngrid){
    u[i]=r[i]*cbet[0];
    for (int k=1; k<ngrid;k++){
      u[k+i]=(r[k+i]-c[k]*u[k+i-1])*cbet[k];
    }
    for (int k=ngrid-2;k>=0;k--){
      u[k+i]-=cwet[k+1]*u[k+i+1];
    }
  }
  return;
}

void FieldSolver::tridagy(vector<complex<double > > &u)
{
  for (int i=0; i<ngrid; i++){
    u[i]=r[i]*cbet[0];
  }
  for (int k=1;k<ngrid-1;k++){
    int n=k*ngrid;
    for (int i=0;i<ngrid;i++){
      u[n+i]=(r[n+i]-c[k]*u[n+i-ngrid])*cbet[k];
    }
  }
  for (int k=ngrid-2;k>=0;k--){
    int n=k*ngrid;
    for (int i=0; i<ngrid; i++){
      u[n+i]-=cwet[k+1]*u[n+i+ngrid];
    }
  }

  return;
}


void FieldSolver::getDiag(double delz,double dgrid, double xks, int ngrid_in)
{

  if (delz==delz_save){
    return;
  }
  delz_save=delz;
  ngrid=ngrid_in;


  double rtmp=0.25*delz/(xks*dgrid*dgrid); //factor dz/(4 ks dx^2)
  cstep = complex<double> ( 0, rtmp );

  double *mupp = new double[ngrid];
  double *mmid = new double[ngrid];
  double *mlow = new double[ngrid];
  complex<double> *cwrk1= new complex<double>[ngrid];
  complex<double> *cwrk2= new complex<double>[ngrid];

  mupp[0]=rtmp;
  mmid[0]=-2*rtmp;
  mlow[0]=0;
  for (int i=1;i<(ngrid-1);i++){
    mupp[i]=rtmp;
    mmid[i]=-2*rtmp;
    mlow[i]=rtmp;
  }
  mupp[ngrid-1]=0;
  mmid[ngrid-1]=-2*rtmp;
  mlow[ngrid-1]=rtmp;

  for (int i=0; i <ngrid; i++){
    cwrk1[i] =complex<double>(0,-mupp[i]);
    cwrk2[i] =complex<double>(1,-mmid[i]);
    c[i]     =complex<double>(0,-mlow[i]);
  }

 
  cbet[0]=1./cwrk2[0];
  cwet[0]=0.;
  for (int i=1; i<ngrid; i++){
    cwet[i]=cwrk1[i-1]*cbet[i-1];
    cbet[i]=1./(cwrk2[i]-c[i]*cwet[i]); 
    
  }

  delete[] mupp;
  delete[] mmid;
  delete[] mlow;
  delete[] cwrk1;
  delete[] cwrk2;
}
