#include "FieldSolver.h"
#include "Field.h"
#include "Beam.h"

// CL, added for crsource term dumping
#include <mpi.h>
#include <hdf5.h>
#include "HDF5_cwc.h"

FieldSolver::FieldSolver()
{
  delz_save=0;
  difffilter_ = false;
  ngrid = 0;
  hasPlan=false;

  cntr_=0;
}

FieldSolver::~FieldSolver(){
#ifdef FFTW
    if (hasPlan){
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
        if (hasPlan) {
            delete[] in;
            delete[] out;
            fftw_destroy_plan(p);
            fftw_destroy_plan(pi);
        }
        sigmoid_.resize(ngrid*ngrid);

        in = new complex<double>[ngrid * ngrid];
        out = new complex<double>[ngrid * ngrid];

        p = fftw_plan_dft_2d(ngrid, ngrid, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out),
                             FFTW_FORWARD, FFTW_MEASURE);
        pi = fftw_plan_dft_2d(ngrid, ngrid, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out),
                              FFTW_BACKWARD, FFTW_MEASURE);
        hasPlan=true;
		
#endif
    }
}

void FieldSolver::initSourceFilter(bool do_filter,double xc, double yc, double sig) {
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
    if (sig <= 0) {  // check for unphysical input
        difffilter_ = false;
    }
    if (difffilter_) {  // pre-calculation of sigmoid function
        int nhalf = (ngrid - 1) / 2;
        int idx = 0;
        for (int iy = 0; iy < ngrid; iy++) {
            double y = static_cast<double>(((iy + nhalf) % ngrid - nhalf)) / static_cast<double>(nhalf) / yc;
            for (int ix = 0; ix < ngrid; ix++) {
                double x = static_cast<double>(((ix + nhalf) % ngrid - nhalf)) / static_cast<double>(nhalf) / xc;
                double r = (sqrt(x * x + y * y) - 1) / sig;
                sigmoid_[idx] = 1. / (1 + exp(r));
                idx++;
            }
        }
    }
#endif
}


void FieldSolver::advance(double delz, Field *field, Beam *beam, Undulator *und)
{
  HDF5_CollWriteCore *pcwc = NULL;
  HDF5_CollWriteCore *pcwc_filt = NULL;
  int mpi_rank,mpi_size;
  hid_t fid;
  int do_dump=0;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int nstotal=mpi_size*field->field.size();
  int smin=mpi_rank*field->field.size();
  int smax=smin+field->field.size();

  cntr_++;
  do_dump = (cntr_%100)==0; // FIXME
  if(do_dump) {
    stringstream ss_fn;

    ss_fn << "x" << "." << cntr_ << ".crsource.h5";
    hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
    if (mpi_size>1){
      H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
    }
    fid=H5Fcreate(ss_fn.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
    H5Pclose(pid);

    const int datadim=3;
    vector<hsize_t> totalsize(datadim,0);
    totalsize[0] = nstotal;
    totalsize[1] = ngrid*ngrid;
    totalsize[2] = 2; /* re/im */

    pcwc = new HDF5_CollWriteCore;
    pcwc_filt = new HDF5_CollWriteCore;
    pcwc->create_and_prepare(fid, "crsource", " ", &totalsize, datadim);
    pcwc_filt->create_and_prepare(fid, "crsource_filtered", " ", &totalsize, datadim);
  }

  for (int ii=0; ii<field->field.size();ii++){  // ii is index for the beam
    int i= (ii+field->first) % field->field.size();           // index for the field
  
    // clear source term
    for (int k=0; k<ngrid*ngrid; k++){
      crsource[k]=0;
    }

    int harm=field->getHarm();
    //   cout << "Harm: " << harm << " Coupling: " << und->fc(harm) << endl;

    /*** construct source term ***/
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
    }
    /*** end of source term construction ***/

    struct dump_settings ds;
    ds.dump_en = do_dump;
    ds.pcwc = pcwc;
    ds.pcwc_filt = pcwc_filt;
    ds.curr_slice = smin+ii;
    ds.nstot = nstotal;

    this->filterSourceTerm(&ds);
    this->ADI(field->field[i]);
  }

  if(do_dump) {
    pcwc->close();
    pcwc_filt->close();
    H5Fclose(fid);
  }
  return;
}

void FieldSolver::dump_crsource(struct dump_settings *pds, HDF5_CollWriteCore *pcwc)
{
    const int datadim=3;
    vector<hsize_t> my_offset(datadim,0); /* is assigned in the loop (depends on current slice to be written) */
    vector<hsize_t> my_count(datadim,0);
    my_count[0] = 1;
    my_count[1] = ngrid*ngrid;
    my_count[2] = 2; /* re/im */


    vector<double> tmp(ngrid*ngrid*2);
    for(int i=0; i<ngrid*ngrid; i++) {
      tmp[2*i  ] = crsource[i].real();
      tmp[2*i+1] = crsource[i].imag();
    }
    tmp[0] = pds->curr_slice; // earmark this slice for debugging purposes

    /* initiate collective write (all processes on MPI communicator write their slice) */
    my_offset[0] = pds->curr_slice;
    pcwc->write(&tmp, &my_count, &my_offset);
}

void FieldSolver::filterSourceTerm(struct dump_settings *pds)
{
    /**
    * filter the source term to avoid emission under strong angle. The source term is transformed via
    * FFTW2D. Spectral components under a large angle are damped according to a sigmoid function.
    * @returns None
    */
    if (!difffilter_) { return; }

#ifdef FFTW
    if(pds->dump_en)
        dump_crsource(pds, pds->pcwc);

    for (int idx=0; idx <ngrid*ngrid;idx++){
        in[idx]=crsource[idx];   // field for the FFT
    }
    fftw_execute(p);
    // the sigmoid function can be pre-calculated
    for (int idx=0; idx <ngrid*ngrid; idx++){
        in[idx]=out[idx]*sigmoid_[idx];
    }
    fftw_execute(pi);
    double norm = 1./static_cast<double>(ngrid*ngrid);
    for (int idx=0; idx <ngrid*ngrid;idx++){
        crsource[idx]=out[idx]*norm;
    }

    if(pds->dump_en)
        dump_crsource(pds, pds->pcwc_filt);
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
