#include "EFieldSolver.h"
#include "Beam.h"


EFieldSolver::EFieldSolver() {
    nz = 0;
    nphi = 0;
    ngrid = 0;
    rmax = 1;
    ks = 1;
    xcen = 0;
    ycen = 0;
    dr = 1;
    rbound =1;
    longrange = false;
    fcurrent.clear();
    fsize.clear();
}

EFieldSolver::~EFieldSolver()= default;

void EFieldSolver::init(double rmax_in, int ngrid_in, int nz_in, int nphi_in, double lambda, bool longr_in) {

    rmax = rmax_in;
    ngrid = ngrid_in;
    nz = nz_in;
    nphi = nphi_in;
    ks = 4 * asin(1) / lambda;
    longrange = longr_in;

    // adjust working arrays
    if (ngrid != csrc.size()) {
        csource.resize(ngrid*nz*(2*nphi+1));
        cfield.resize(ngrid*nz*(2*nphi+1));
        for (unsigned long i ; i < ngrid*nz*(2*nphi+1); i ++){
            cfield[i] = complex<double>(0,0);
        }
        vol.resize(ngrid);
        lmid.resize(ngrid);
        llow.resize(ngrid);
        lupp.resize(ngrid);
        celm.resize(ngrid);
        csrc.resize(ngrid);
        clow.resize(ngrid);
        cmid.resize(ngrid);
        cupp.resize(ngrid);
        gam.resize(ngrid);
        rlog.resize(ngrid);
    }

}

void EFieldSolver::longRange(Beam *beam, double gamma0, double aw) {
    auto nsize = beam->beam.size();
    // check if the beam size has been changed:
    if (nsize != beam->longESC.size()){
        beam->longESC.resize(nsize);
    }

    for (int i =0; i < nsize; i++){
        beam->longESC[i]=0;
    }
    if (!longrange) {
        return;
    }
    double gamma = gamma0/sqrt(1+aw*aw);

    int MPIsize=1;
    int MPIrank=0;
    if (!MPISingle){
        MPI_Comm_size(MPI_COMM_WORLD, &MPIsize); // assign rank to node
        MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank); // assign rank to node
    }

    // resize if needed also after harmonic conversion.
    if (fcurrent.size()!=(MPIsize*nsize)){
        fcurrent.resize(MPIsize*nsize);
        fsize.resize(MPIsize*nsize);
        work1.resize(nsize);
        work2.resize(nsize);
    }

    // fill local slices
    for (int i=0; i<nsize;i++){
        work1[i]=beam->current[i];
        work2[i]=beam->getSize(i);
    }

    // gathering information on full current and beam size profile
    MPI_Allgather(&work1.front(),nsize,MPI_DOUBLE,&fcurrent.front(),nsize,MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&work2.front(),nsize,MPI_DOUBLE, &fsize.front(),nsize,MPI_DOUBLE, MPI_COMM_WORLD);

    double scl = beam->slicelength/2./asin(1)/299792458.0/2/8.85e-12;  // convert to units of electron rest mass.

    for (int i=0; i < nsize; i++){
        double EFld = 0;
        auto isplit = nsize*MPIrank+i;
        for (int j = 0; j < isplit; j++) {
            double ds = static_cast<double>(j - isplit) * beam->slicelength*gamma;
            double coef = 1 +ds/sqrt(ds*ds+fsize[j]); // ds is negative
            EFld += coef*scl*fcurrent[j]/fsize[j];
        }
        for (auto j = isplit+1; j < MPIsize*nsize; j++) {
            double ds = static_cast<double>(j - isplit) * beam->slicelength*gamma;
            double coef = 1 -ds/sqrt(ds*ds+fsize[j]); // ds is negative
            EFld -= coef*scl*fcurrent[j]/fsize[j];
        }
        beam->longESC[i] = EFld;
    }
}

auto EFieldSolver::analyseBeam(vector<Particle> *beam){
/*
 *  calculates the center of the beam slice and its extension (3 times rms radius to construct the space charge grid
 */
    auto npart = beam->size();
    if (npart <1) {return npart;}

    // calculate center of beam slice
    xcen = 0;
    ycen = 0;
    for (int i = 0; i < npart; i++) {
        xcen += beam->at(i).x;
        ycen += beam->at(i).y;
    }
    xcen /= static_cast<double>(npart);
    ycen /= static_cast<double>(npart);

    // calculate radial extension of particle distribution
    rbound = 0;
    for (int i = 0; i < npart; i++) {
        double tx = beam->at(i).x - xcen;
        double ty = beam->at(i).y - ycen;
        rbound  += tx * tx + ty * ty;
    }
    rbound = sqrt(rbound/static_cast<double>(npart))*3; // calculating the rms radius and then assume 3 sigma for full beam
    return npart;
}

void EFieldSolver::constructLaplaceOperator(){
    // r_j are the center grid points. The boundaries are given +/-dr/2. Thus r_j=(j+1/2)*dr
    for (int j = 0; j < ngrid; j++) {
        vol[j] = dr * dr * (2 * j + 1);
        lupp[j] = 2 * (j + 1) / vol[j];  // Eq 3.30 of my thesis
        llow[j] = 2 * j / vol[j];      // Eq 3.29
        if (j == 0) {
            rlog[j] = 0;
        } else {
            rlog[j] = 2 * log((j + 1) / j) / vol[j];  // Eq. 3.31
        }
    }
    lupp[ngrid - 1] = 0;
}


double EFieldSolver::getEField(double x,double y,double theta) {

    double ez = 0;
    if (!this->hasShortRange()) {
        return ez;
    }
    auto tx = x - xcen;
    auto ty = y - ycen;
    auto r = sqrt(tx * tx + ty * ty);
    if (r > rmax * rbound) {
        return ez;
    }
    auto phi = atan2(ty, tx);
    auto ir = static_cast<unsigned long>(floor(r / dr));

    for (int l = 1 ; l < nz+1; l++) {
        for (int m = -nphi; m<=nphi; m++){
            auto phase = l * theta + m * phi;
            unsigned long idx = (l-1)*ngrid*(2*nphi+1)+(m+nphi)*ngrid;
            auto ctmp = cfield[idx+ir] * complex<double>(cos(phase), sin(phase));
            ez += 2 * ctmp.real();
        }
    }
    return ez;
}




void EFieldSolver::shortRange(vector<Particle> *beam, double current, double gz2) {

    // calculate center of beam slice and its extension
    auto npart = this->analyseBeam(beam);
    if (npart == 0) { return; }
    if (!this->hasShortRange()) { return; }

    // get radial grid spacing
    dr = rmax * rbound / (ngrid - 1); // calculate grid spacing
    this->constructLaplaceOperator(); // calculates the tri-diagonal laplace matrix of space charge field equation

    // construct source term
    double pi = 2. * asin(1);
    double xcuren = current;
    double coef = vacimp / eev * xcuren / static_cast<double>(npart) / pi;
    unsigned long nnphi = 2 * nphi + 1;
    // clear source term first
    for (unsigned long i; i < nz * nnphi * ngrid; i++) {
        csource[i] = complex<double>(0, 0);
    }

    for (auto part: *beam) {
        auto tx = part.x - xcen;
        auto ty = part.y - ycen;
        auto r = sqrt(tx * tx + ty * ty);
        if (r < rmax * rbound) {
            auto phi = atan2(ty, tx);
            auto ir = static_cast<unsigned long>(floor(r / dr));
            auto theta = part.theta;
            for (int l = 1; l < nz + 1; l++) {
                unsigned int idx0 = (l - 1) * nnphi * ngrid;
                for (int m = -nphi; m <= nphi; m++) {
                    auto idx = idx0 + (m + nphi) * ngrid +
                               ir;  // getting index of radial grid nested in loop over azimuthal angle phi and longitudinal harmonics
                    auto phase = l * theta + m * phi;
                    csource[idx] += complex<double>(cos(phase), -sin(phase));
                }
            }
        }
    }
    for (unsigned long j = 0; j < nz * nnphi; j++) {
        for (unsigned long i = 0; i < ngrid; i++) {
            csource[i + j * ngrid] *= complex<double>(0, coef / vol[i]);
        }
    }

    // solving field equation

    for (int l = 1; l < nz + 1; l++) {                 // loop over longitudinal modes
        for (int i = 1; i < ngrid; i++) {
            lmid[i] = -llow[i] - lupp[i] - l * l * ks * ks / gz2; // construct the diagonal terms for m=0
        }
        // solve tridiag(clow,cmid,cupp,csrc,celm,ngrid);
        for (int m = -nphi; m <= nphi; m++) {
            complex<double> bet = cmid[0];
            celm[0] = csrc[0] / bet;
            for (int i = 1; i < ngrid; i++) {
                gam[i] = cupp[i - 1] / bet;
                bet = cmid[i] - clow[i] * gam[i];
                celm[i] = (csrc[i] - clow[i] * celm[i - 1]) / bet;
            }
            for (int i = ngrid - 2; i > -1; i--) {
                celm[i] = celm[i] - gam[i + 1] * celm[i + 1];
            }

            auto idx = ((l - 1) * (2 * nphi + 1) + (m + nphi)) * ngrid;
            for (int i = 0; i < ngrid; i++) {
                cfield[idx + i] = celm[i];
            }
        }
    }
}

