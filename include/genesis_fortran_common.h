#ifndef __GENESIS_FORTRAN_COMMON__
#define __GENESIS_FORTRAN_COMMON__

#define npmax 1024
#define ncmax 513


extern "C"{

// common_blocks

  extern struct{
    double gamma[npmax],theta[npmax],xpart[npmax],ypart[npmax];
    double px[npmax],py[npmax],btpar[npmax],k2gg[npmax],k2pp[npmax];
    double k3gg[npmax],k3pp[npmax],ez[npmax],cpart1[2*npmax],xcuren;
    int npart,nbins;
  } beamcom_;

  extern struct{
    double aw0,xkx,xky,xku,awd,awdx,awdy,qfld,qdx,qdy,gammaref;
  } wigcom_;

  extern struct{
    double crfield[2*ncmax*ncmax],crsource[2*ncmax*ncmax];
    double crhm[2*ncmax*ncmax],cwet[2*ncmax],cbet[2*ncmax];
    double crmatc[2*ncmax];
    double xks;
    int ncar;
  } fieldcom_;
 

// fortran functions

  int track_(double *);
  int pushp_(double *,double *);
  int getdiag_(double *, double *, double *,int *);
  int field_(double*);
  double bessj_(int *,double *);


}

#endif
