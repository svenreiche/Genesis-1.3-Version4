#include "FieldSolver.h"
#include "Field.h"
#include "Beam.h"

FieldSolver::FieldSolver()= default;

FieldSolver::~FieldSolver()= default;

void FieldSolver::advance(double delz, Field *field, Beam *beam, Undulator *und) {

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
        this->ADI(field->field[i]);
    }
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

}


void FieldSolver::tridagx(vector<complex<double > > &u) {
    for (int i = 0; i < ngrid * (ngrid - 1); i += ngrid) {
        u[i] = r[i] * cbet[0];
        for (int k = 1; k < ngrid; k++) {
            u[k + i] = (r[k + i] - c[k] * u[k + i - 1]) * cbet[k];
        }
        for (int k = ngrid - 2; k >= 0; k--) {
            u[k + i] -= cwet[k + 1] * u[k + i + 1];
        }
    }
}

void FieldSolver::tridagy(vector<complex<double > > &u) {
    for (int i = 0; i < ngrid; i++) {
        u[i] = r[i] * cbet[0];
    }
    for (int k = 1; k < ngrid - 1; k++) {
        int n = k * ngrid;
        for (int i = 0; i < ngrid; i++) {
            u[n + i] = (r[n + i] - c[k] * u[n + i - ngrid]) * cbet[k];
        }
    }
    for (int k = ngrid - 2; k >= 0; k--) {
        int n = k * ngrid;
        for (int i = 0; i < ngrid; i++) {
            u[n + i] -= cwet[k + 1] * u[n + i + ngrid];
        }
    }
}


void FieldSolver::getDiag(double delz,double dgrid, double xks, int ngrid_in) {

    if (delz == delz_save) {
        return;
    }
    delz_save = delz;
    ngrid = ngrid_in;


    double rtmp = 0.25 * delz / (xks * dgrid * dgrid); //factor dz/(4 ks dx^2)
    cstep = complex<double>(0, rtmp);

    auto *mupp = new double[ngrid];
    auto *mmid = new double[ngrid];
    auto *mlow = new double[ngrid];
    auto *cwrk1 = new complex<double>[ngrid];
    auto *cwrk2 = new complex<double>[ngrid];
    if (c.size() != ngrid) {
        c.resize(ngrid);
        r.resize(ngrid * ngrid);
        cbet.resize(ngrid);
        cwet.resize(ngrid);
        crsource.resize(ngrid * ngrid);
    }

    mupp[0] = rtmp;
    mmid[0] = -2 * rtmp;
    mlow[0] = 0;
    for (int i = 1; i < (ngrid - 1); i++) {
        mupp[i] = rtmp;
        mmid[i] = -2 * rtmp;
        mlow[i] = rtmp;
    }
    mupp[ngrid - 1] = 0;
    mmid[ngrid - 1] = -2 * rtmp;
    mlow[ngrid - 1] = rtmp;

    for (int i = 0; i < ngrid; i++) {
        cwrk1[i] = complex<double>(0, -mupp[i]);
        cwrk2[i] = complex<double>(1, -mmid[i]);
        c[i] = complex<double>(0, -mlow[i]);
    }


    cbet[0] = 1. / cwrk2[0];
    cwet[0] = 0.;
    for (int i = 1; i < ngrid; i++) {
        cwet[i] = cwrk1[i - 1] * cbet[i - 1];
        cbet[i] = 1. / (cwrk2[i] - c[i] * cwet[i]);

    }

    delete[] mupp;
    delete[] mmid;
    delete[] mlow;
    delete[] cwrk1;
    delete[] cwrk2;
}
