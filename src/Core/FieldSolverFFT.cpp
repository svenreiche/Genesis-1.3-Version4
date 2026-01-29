#include "FieldSolverFFT.h"
#include "Field.h"
#include "Beam.h"



FieldSolverFFT::~FieldSolverFFT() = default;

void FieldSolverFFT::advance(double delz, Field *field, Beam *beam, Undulator *und) {

    for (unsigned long ii = 0; ii < field->field.size(); ii++) {  // ii is index for the beam

        // clear source term
        for (int ig = 0; ig < ngrid * ngrid; ig++) {
            crsource[ig] = 0;
        }

        // constructing source term
        int harm = field->getHarm();
        if (und->inUndulator() && field->isEnabled() && (harm % 2 == 1)) { // do not need to calculate for even harmonics
            double scl = und->fc(harm) * vacimp * beam->current[ii] * field->xks * delz;
            scl /= 4 * eev * static_cast<double>(beam->beam[ii].size()) * field->dgrid * field->dgrid;
            complex<double> cpart;
            double part, weight, wx, wy;
            int idx;

            for (auto & particle : beam->beam.at(ii)) {
                double x = particle.x;
                double y = particle.y;
                double theta = static_cast<double>(harm) * particle.theta;
                double gamma = particle.gamma;

                if (field->getLLGridpoint(x, y, &wx, &wy, &idx)) {

                    part = sqrt(und->faw2(x, y)) * scl / gamma;
                    // tmp  should be also normalized with beta parallel
                    cpart = complex<double>(sin(theta), cos(theta)) * part;

                    weight = wx * wy;
                    crsource[idx] += weight * cpart;
                    weight = (1 - wx) * wy;
                    idx++;
                    crsource[idx] += weight * cpart;
                    weight = wx * (1 - wy);
                    idx += ngrid - 1;
                    crsource[idx] += weight * cpart;
                    weight = (1 - wx) * (1 - wy);
                    idx++;
                    crsource[idx] += weight * cpart;
                }
            }
        }  // end of source term construction

        unsigned long i = (ii + field->first) % field->field.size();           // index for the field

        // get the FFT representation of the radiation field and the source term


        this->FFT(field->field[i]);
    }
}


void FieldSolverFFT::FFT(vector<complex<double> > &crfield)
{
    // Do the FFT of the field and source term
    for (unsigned long ii = 0; ii < crfield.size(); ii++) {
        in[ii] = crfield[ii];
    }
#ifdef FFTW
        fftw_execute(p);
#endif
    for (unsigned long ii = 0; ii < crfield.size(); ii++) {
        uf[ii] = out[ii];
        in[ii] = crsource[ii];
    }
#ifdef FFTW
        fftw_execute(p);
#endif

    for (unsigned long ii = 0; ii < crfield.size(); ii++) {
        sf[ii] = out[ii];
    }
    // filter source term
    if (doFilter_) {
        for (unsigned long ii = 0; ii < crfield.size(); ii++) {
            sf[ii]*=sigmoid_[ii];
        }
    }

    // do the actual propagation
    for (unsigned long ii = 0; ii < crfield.size(); ii++) {
        in[ii]=uf[ii]*exp(K2[ii]*delz_save) + 2.*sf[ii]; // - complex<double>(0,1.) *sf[ii];
    }
#ifdef FFTW
        fftw_execute(ip);
#endif
    double norm = 1./static_cast<double>(ngrid*ngrid);
    for (unsigned long ii = 0; ii < crfield.size(); ii++) {
        crfield[ii]=out[ii]*norm;
    }
}


void FieldSolverFFT::init(double delz,double dgrid, double xks, unsigned int ngrid_in) {

    delz_save = delz;
    if (!hasPlan) {
        ks = xks;
        ngrid = ngrid_in;
        dk = 4.*asin(1.)/(static_cast<double>(ngrid)*dgrid);
        in = new complex<double> [ngrid*ngrid];
        out= new complex<double> [ngrid*ngrid];
        uf.resize(ngrid*ngrid);
        sf.resize(ngrid*ngrid);
        K2.resize(ngrid*ngrid);
        sigmoid_.resize(ngrid*ngrid);

        double shift=-0.5*static_cast<double> (ngrid-1);
        for (int iy=0;iy<ngrid;iy++) {
            double dy=static_cast<double>(iy)+shift;
            double y = dy / static_cast<double>(ngrid) /yc ;
            for (int ix=0;ix<ngrid;ix++) {
                double dx=static_cast<double>(ix)+shift;
                double x = dx / static_cast<double>(ngrid) /xc ;
                int iiy=(iy+(ngrid+1)/2) % ngrid;
                int iix=(ix+(ngrid+1)/2) % ngrid;
                int ii=iiy*ngrid+iix;
                K2[ii] = complex<double>(0,-(dx*dx+dy*dy)*dk*dk/2./xks);
                double r = (sqrt(x * x + y * y) - 1) / sig;
                sigmoid_[ii] = 1. / (1 + exp(r));
            }
        }
        crsource.resize(ngrid* ngrid);


#ifdef FFTW
        p  = fftw_plan_dft_2d(ngrid,ngrid,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out),FFTW_FORWARD,FFTW_MEASURE);
        ip  = fftw_plan_dft_2d(ngrid,ngrid,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out),FFTW_BACKWARD,FFTW_MEASURE);
#endif
        hasPlan = true;
    }
}

void FieldSolverFFT::initSourceFilter(double xc_in, double yc_in, double sig_in,bool do_filter) {
    doFilter_ = do_filter;
    if (sig <= 0) {  // check for unphysical input
        doFilter_ = false;
    }
    xc=xc_in;
    yc=yc_in;
    sig=sig_in;
};