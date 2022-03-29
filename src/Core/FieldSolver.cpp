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

  crsource_dump_en_ = false;
  crsource_dump_every_ = 100;
  crsource_dump_rootname_ = "__";
  call_cntr_adv_ = 0; // variable counting the number of source term constructions (it is used for crsource dumping, if requested by user)
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

        unsigned int flags = FFTW_MEASURE;
        // DBG: if you want to have absolutely reproducible FFT results
        // flags = FFTW_ESTIMATE;
        p = fftw_plan_dft_2d(ngrid, ngrid,
               reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out),
               FFTW_FORWARD,  flags);
        pi = fftw_plan_dft_2d(ngrid, ngrid,
               reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out),
               FFTW_BACKWARD, flags);
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
void FieldSolver::initSourceFilter_DbgDumpSettings(bool dump_en_in, int step_in, string rootname_in) {
  // ignore illegal values for downscaler
  if(step_in<=0)
    return;

  crsource_dump_every_ = step_in;
  crsource_dump_rootname_ = rootname_in;
  crsource_dump_en_ = dump_en_in;
}


void FieldSolver::advance(double delz, Field *field, Beam *beam, Undulator *und)
{
  bool dump_at_this_step=false;
  struct dump_settings ds;

  // Write crsource dump file for this integration step?
  // Note that the result needs to be identical for all processes
  // on the MPI communicator (all nodes have to participate in HDF5
  // parallel I/O), otherwise the program will freeze
  call_cntr_adv_++;
  dump_at_this_step=false;
  if(crsource_dump_en_ && difffilter_)
    dump_at_this_step = ((call_cntr_adv_ % crsource_dump_every_) == 0);

  /* if crsource data is to be dumped, set up output file and collective I/O */
  ds.do_dump = dump_at_this_step;
  if(dump_at_this_step) {
    dump_file_open(&ds, field);
  }

  for (int ii=0; ii<field->field.size();ii++){  // ii is index for the beam
    int i= (ii+field->first) % field->field.size();           // index for the field
  
    // clear source term
    for (int k=0; k<ngrid*ngrid; k++){
      crsource[k]=0;
    }

    int harm=field->getHarm();
    //   cout << "Harm: " << harm << " Coupling: " << und->fc(harm) << endl;

    /*** construct source term for this slice ***/
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
    /*** end of source term construction for this slice ***/

    ds.curr_slice = ds.smin + ii;
    this->filterSourceTerm(&ds);
    this->ADI(field->field[i]);
  }

  if(dump_at_this_step) {
    dump_file_close(&ds);
  }
  return;
}



void FieldSolver::dump_file_open(struct dump_settings *pds, Field *field)
{
  HDF5_CollWriteCore *pcwc = NULL;
  HDF5_CollWriteCore *pcwc_filt = NULL;
  int mpi_rank,mpi_size;
  hid_t fid;
  stringstream ss_fn;
  string fn;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int nstotal=mpi_size*field->field.size();
  int smin=mpi_rank*field->field.size();
  int smax=smin+field->field.size();
  pds->smin=smin;
  pds->nstot=nstotal;

  // construct file name and open file for parallel write access
  ss_fn << crsource_dump_rootname_ << "." << call_cntr_adv_;
  if(field->getHarm()>1) {
    ss_fn << ".h" << field->getHarm();
  }
  ss_fn << ".crsource.h5";
  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  if (mpi_size>1){
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  fn = ss_fn.str();
  fid=H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
  H5Pclose(pid);
  if(mpi_rank==0) {
    cout << "Info: opened crsource dump file "<<fn <<endl;
  }

  const int datadim=3;
  vector<hsize_t> totalsize(datadim,0);
  totalsize[0] = nstotal;
  totalsize[1] = ngrid*ngrid;
  totalsize[2] = 2; /* re/im */

  pcwc = new HDF5_CollWriteCore;
  pcwc_filt = new HDF5_CollWriteCore;
  pcwc->create_and_prepare(fid, "crsource", " ", &totalsize, datadim);
  pds->pcwc = pcwc;
  pcwc_filt->create_and_prepare(fid, "crsource_filtered", " ", &totalsize, datadim);
  pds->pcwc_filt = pcwc_filt;

  dump_filter(pds, fid);

  pds->fid = fid;
}
void FieldSolver::dump_file_close(struct dump_settings *pds)
{
  int mpi_rank;

  pds->pcwc->close();
  pds->pcwc_filt->close();
  H5Fclose(pds->fid);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if(mpi_rank==0) {
    cout << "Info: closing crsource dump file" << endl;
  }
}

void FieldSolver::dump_crsource(struct dump_settings *pds, HDF5_CollWriteCore *pcwc)
{
    if(!pds->do_dump)
      return;


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
    // tmp[0] = pds->curr_slice; // label slice data for debugging purposes

    /* initiate collective write (all processes on MPI communicator write their slice) */
    my_offset[0] = pds->curr_slice;
    pcwc->write(&tmp, &my_count, &my_offset);
}

void FieldSolver::dump_filter(struct dump_settings *pds, hid_t pobj)
{
    if(!pds->do_dump)
      return;

    HDF5_CollWriteCore cwc;
    const int filt_datadim=1;
    vector<hsize_t> filt_totalsize(filt_datadim,0);
    vector<hsize_t> my_offset(filt_datadim,0);
    vector<hsize_t> my_count(filt_datadim,0);
    filt_totalsize[0] = ngrid*ngrid; // sigmoid is vector<double>

    cwc.create_and_prepare(pobj, "filter", " ", &filt_totalsize, filt_datadim);

    // Only process with rank=0 is writing (controlled by setting count to zero on all others).
    // Note that we still need to call write member of HDF5_CollWriteCore on all processes.
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(mpi_rank==0)
      my_count[0] = ngrid*ngrid;
    else
      my_count[0] = 0;

    // sigmoid_[0] = 1+mpi_rank; // test code: which node is writing the data???
    cwc.write(&sigmoid_, &my_count, &my_offset);
    cwc.close();
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
    if(pds->do_dump)
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

    if(pds->do_dump)
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
