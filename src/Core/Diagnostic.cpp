//
// Created by reiche on 06.01.22.
//

#include <iostream>
#include <complex>

#include "Diagnostic.h"

// constructor/destructor
Diagnostic::Diagnostic(){}
Diagnostic::~Diagnostic(){}

void Diagnostic::init(int nz_in, int ns_in){

    nz = nz_in;
    ns = ns_in;
    val.clear();
    units.clear();
    for (auto group : dbeam){
        for (auto const &tag : group->getTags()){
            units[tag.first] = tag.second;
            val[tag.first].resize(ns*nz);
        }
    }
    iz=0;                              // set the counter
}

void Diagnostic::calc(Beam* beam,std::vector<Field*> *field,Undulator* und) {
    for (auto group: dbeam){
        group->getValues(beam,val,iz);
    }
    iz++;
}



//-----------------------------------------------------------
// actual implementation of diagnostic calculation

void DiagBeam::defineTags() {
    nharm = 2;
    npar = 8+nharm*2;
    tags.clear();
    tags["energy"]="mc^2";
    tags["energyspread"]="mc^2";
    tags["xposition"]="m";
    tags["xsize"]="m";
    tags["yposition"]="m";
    tags["ysize"]="m";
    tags["pxposition"]="rad";
    tags["pyposition"]="rad";
    tags["bunching"]=" ";
    tags["bunchingphase"]=" ";
    char buff[30];
    for (int iharm=1; iharm<nharm;iharm++){
        snprintf(buff,sizeof(buff),"bunching%d",iharm+1);
        tags[buff]=" ";
        snprintf(buff,sizeof(buff),"bunchingphase%d",iharm+1);
        tags[buff]="rad";
    }
}


void DiagBeam::getValues(Beam *beam,std::map<std::string,std::vector<double> >&val, int iz)
{

    int ns = beam->beam.size();
    int is=0;
    std::vector<complex<double> > b(nharm);

    for (auto const &slice :beam->beam){
        double g1=0;
        double g2=0;
        double x1=0;
        double x2=0;
        double y1=0;
        double y2=0;
        double px1=0;
        double py1=0;
        for (int iharm;iharm<nharm;iharm++){
            b[iharm] = 0;
        }
        for (auto const &par: slice){
           x1 += par.x;
           x2 += par.x*par.x;
           y1 += par.y;
           y2 += par.y*par.y;
           g1 += par.gamma;
           g2 += par.gamma*par.gamma;
           px1 += par.px;
           py1 += par.py;
           complex<double> phasor = complex<double> (cos(par.theta),sin(par.theta));
           b[0] = phasor;
           for(int iharm=1; iharm<nharm;iharm++){
               b[iharm]=b[iharm-1]*phasor;
           }
        }
        double norm=1.;
        if (slice.size()>0) {
            norm = 1. / static_cast<double>(slice.size());
        }
        int idx = iz*ns+is;
        val["energy"][idx]=g1*norm;  // make sure that this is the correct order as listed in the variable 'tags'
        val["energyspread"][idx]=sqrt(fabs(g2*norm-g1*g1*norm*norm));
        val["xsize"][idx]=sqrt( fabs(x2*norm-x1*x1*norm*norm));
        val["xposition"][idx]=x1*norm;
        val["ysize"][idx]=sqrt( fabs(y2*norm-y1*y1*norm*norm));
        val["yposition"][idx]=y1*norm;
        val["pxposition"][idx]=px1*norm;
        val["pyposition"][idx]=py1*norm;
        val["bunching"][idx] =abs(b[0]);
        val["bunchingphase"][idx] =atan2(b[0].imag(),b[0].real());
        char buff[100];
        for (int iharm = 1; iharm < nharm; iharm++){
            snprintf(buff,sizeof(buff),"bunching%d",iharm+1);
            b[iharm]*=norm;
            val[buff][idx]=abs(b[iharm]);
            snprintf(buff,sizeof(buff),"bunchingphase%d",iharm+1);
            val[buff][idx]=atan2(b[iharm].imag(),b[iharm].real());
        }
    }
}