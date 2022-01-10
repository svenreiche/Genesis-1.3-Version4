//
// Created by reiche on 06.01.22.
//

#include <iostream>
#include <complex>

#include "Diagnostic.h"
//----------------------------------------------------------
// main class which manage all the calculation on the beam, field or undulator.
// It allocates the memory and define the interface to the output class so that the
// diagnostics classes are only concern with the actual calculation.

// routine to register the output parameters and to allocate memory
void Diagnostic::init(int nz_in, int ns_in, int nfld){

    FilterDiagnostics filter;
    nz = nz_in;
    ns = ns_in;
    val.clear();
    units.clear();

    val.resize(1+nfld);   // at least the beam
    units.resize(1+nfld);

    zout.resize(nz);

    // loop through the different class for registration and allocate memory
    for (auto group : dbeam){
        for (auto const &tag : group->getTags(filter)){
            units.at(0)[tag.first] = tag.second.units;   // accumulate all units
            int size = ns * nz;
            if (tag.second.once) { size /= nz; }
            if (tag.second.global) {size /= ns; }
            val.at(0)[tag.first].resize(size);
        }
    }

    for (int ifld=0; ifld<nfld;ifld++){
        for (auto group : dfield){
            for (auto const &tag : group->getTags(filter)){
                units.at(1+ifld)[tag.first] = tag.second.units;   // accumulate all units
                int size = ns * nz;
                if (tag.second.once) { size /= nz; }
                if (tag.second.global) {size /= ns; }
                val.at(1+ifld)[tag.first].resize(size);
            }
        }
    }

    iz=0;                              // set the counter
}

// wrapper to do all the diagnostics calculation at a given integration step iz.
void Diagnostic::calc(Beam* beam,std::vector<Field*> *field,Undulator* und) {
    zout[iz]=und->getz();
    for (auto group: dbeam){
        group->getValues(beam,val[0],iz);
    }
    iz++;

    for (int ifld=0; ifld < field->size(); ifld++){
        for (auto group: dfield) {
            group->getValues(field->at(ifld),val[1+ifld],iz);
        }
    }
}


//-----------------------------------------------------------
// actual implementation of diagnostic calculation - beam

std::map<std::string,OutputInfo> DiagBeam::getTags(FilterDiagnostics & filter_in) {

    tags.clear();
    filter.clear();

    nharm = filter_in.beam.harm;
    global = filter_in.beam.global;

    if (filter_in.beam.energy) {
        filter["energy"] = true;
        tags["energy"] = {false, false, "mc^2"};
        tags["energyspread"] = {false, false, "mc^2"};
        if (global) {
            tags["global/energy"] = {true, false, "mc^2"};
            tags["global/energyspread"] = {true, false, "mc^2"};
        }
    }
    if (filter_in.beam.spatial) {
        filter["spatial"] = true;
        tags["xposition"] = {false, false, "m"};
        tags["xsize"] = {false, false, "m"};
        tags["yposition"] = {false, false, "m"};
        tags["ysize"] = {false, false, "m"};
        tags["pxposition"] = {false, false, "rad"};
        tags["pyposition"] = {false, false, "rad"};
        if (global) {
            tags["global/xposition"] = {true, false, "m"};
            tags["global/xsize"] = {true, false, "m"};
            tags["global/yposition"] = {true, false, "m"};
            tags["global/ysize"] = {true, false, "m"};
        }
    }
    // bunching as one of the fundamental parameter is always on
    tags["bunching"] = {false, false, " "};
    tags["bunchingphase"] = {false, false, "rad"};
    char buff[30];
    for (int iharm = 1; iharm < nharm; iharm++) {
        snprintf(buff, sizeof(buff), "bunching%d", iharm + 1);
        tags[buff] = {false, false, " "};
        snprintf(buff, sizeof(buff), "bunchingphase%d", iharm + 1);
        tags[buff] = {false, false, "rad"};
    }
    if (filter_in.beam.auxiliar) {
        filter["aux"] = true;
        tags["efield"] = {false, false, "eV/m"};
    }
    if (filter_in.beam.current) {
        tags["current"] = {false, false, "A"};
    } else{
        tags["current"]= {false,true,"A"};
    }
    tags["emitx"]={false,true,"m"};
    tags["emity"]={false,true,"m"};
    tags["betax"]={false,true,"m"};
    tags["betay"]={false,true,"m"};
    tags["alphax"]={false,true,"rad"};
    tags["alphay"]={false,true,"rad"};

    return tags;
}

void DiagBeam::getValues(Beam *beam,std::map<std::string,std::vector<double> >&val, int iz)
{

    int ns = beam->beam.size();
    int is=0;
    std::vector<complex<double> > b(nharm);
    double g_cur=0;
    double g_g1=0;
    double g_g2=0;
    double g_x1=0;
    double g_x2=0;
    double g_y1=0;
    double g_y2=0;

    for (auto const &slice :beam->beam){
        double g1=0;
        double g2=0;
        double x1=0;
        double x2=0;
        double y1=0;
        double y2=0;
        double px1=0;
        double py1=0;
        double px2=0;
        double py2=0;
        double xpx=0;
        double ypy=0;

        for (int iharm=0; iharm<nharm; iharm++){
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
           px2 += par.px*par.px;
           py2 += par.py*par.py;
           xpx += par.x*par.px;
           ypy += par.y*par.py;
           complex<double> phasor = complex<double> (cos(par.theta),sin(par.theta));
           b[0] = phasor;
           for(int iharm=1; iharm<nharm;iharm++){
               b[iharm]=b[iharm-1]*phasor;
           }
        }
        double norm=1.;
        if (!slice.empty()) {
            norm = 1. / static_cast<double>(slice.size());
        }
        // normalize the values
        x1 *= norm;
        x2 *= norm;
        y1 *= norm;
        y2 *= norm;
        px1 *= norm;
        px2 *= norm;
        py1 *= norm;
        py2 *= norm;
        g1 *= norm;
        g2 *= norm;
        xpx *=norm;
        ypy *= norm;
        for (int iharm=0; iharm<nharm;iharm++){
            b[iharm]*=norm;
        }

        //-------------------------------------------------------------------------------
        // save into the allocated memory space
        int idx = iz*ns+is;         // index for saving the data
        if (filter["energy"]){
            if (val.find("energy") != val.end()) {val["energy"][idx]=g1;}
            if (val.find("energyspread") != val.end()) {val["energyspread"][idx]=sqrt(fabs(g2-g1*g1));}
        }
        if (filter["spatial"]){
            if (val.find("xposition") != val.end()){val["xposition"][idx]=x1;}
            if (val.find("xsize") != val.end()){val["xsize"][idx]=sqrt( fabs(x2-x1*x1));}
            if (val.find("yposition") != val.end()){val["yposition"][idx]=y1;}
            if (val.find("xsize") != val.end()){val["ysize"][idx]=sqrt( fabs(y2-y1*y1));}
            if (val.find("pxposition") != val.end()){val["pxposition"][idx]=px1;}
            if (val.find("pyposition") != val.end()){ val["pyposition"][idx]=py1;}
        }

        if (val.find("bunching") != val.end()) {val["bunching"][idx] =abs(b[0]);}
        if (val.find("bunchingphase") != val.end()) {val["bunchingphase"][idx] =atan2(b[0].imag(),b[0].real());}
        char buff[100];
        for (int iharm = 1; iharm < nharm; iharm++) {
            snprintf(buff, sizeof(buff), "bunching%d", iharm + 1);
            if (val.find(buff) != val.end()) { val[buff][idx] = abs(b[iharm]); }
            snprintf(buff, sizeof(buff), "bunchingphase%d", iharm + 1);
            if (val.find(buff) != val.end()) { val[buff][idx] = atan2(b[iharm].imag(), b[iharm].real()); }
        }
        if (filter["auxiliar"]){
            if (val.find("efield") != val.end()) {val["efield"][idx]=beam->eloss[is];}
        }
        // here are all the values which are only evaluated once at the bieginning of the run with iz = 0
        if (tags["current"].once){
            if ((iz ==0) and (val.find("current") != val.end())) { val["current"][is] = beam->current[is]; }
        } else {
            if (val.find("current") != val.end()) { val["current"][idx] = beam->current[is]; }
        }
        if (iz == 0) {
            if (val.find("current") != val.end()) { val["current"][is] = beam->current[is]; }
            // because genesis works with momenta and not divergence, the emittance does not need energy
            double ex=sqrt(fabs((x2-x1*x1)*(px2-px1*px1)-(xpx-x1*px1)*(xpx-x1*px1)));
            double ey=sqrt(fabs((y2-y1*y1)*(py2-py1*py1)-(ypy-y1*py1)*(ypy-y1*py1)));
            if (val.find("emitx") != val.end()) { val["emitx"][is] = ex; }
            if (val.find("emity") != val.end()) { val["emity"][is] = ey; }
            if (val.find("betax") != val.end()) { val["betax"][is] = (x2-x1*x1)/ex*g1; }
            if (val.find("betay") != val.end()) { val["betay"][is] = (y2-y1*y1)/ey*g1; }
            if (val.find("alphax") != val.end()) { val["alphax"][is] = (xpx-x1*px1)/ex;}
            if (val.find("alphay") != val.end()) { val["alphay"][is] = (ypy-y1*py1)/ey; }
        }
        //-------------------------------------------------------------------------
        // gather moments for all slices with current as weighting factor
        if (global) {
            g_cur += beam->current[is];
            g_g1 += beam->current[is] * g1;
            g_g2 += beam->current[is] * g2;
            g_x1 += beam->current[is] * g1;
            g_x2 += beam->current[is] * g2;
            g_y1 += beam->current[is] * g1;
            g_y2 += beam->current[is] * g2;
        }
        // increment slice counter
        is++;
    }
    //-------------------------------------------------
    // complete gathering from all nodes if MPI size is larger than one
    if (global){
        int size = 1;
        if (!MPISingle) { MPI_Comm_size(MPI_COMM_WORLD, &size); }  // for future functionality to do scan functionalility of time-dependent runs!
                                                                    // this will be better controlled in the future with a dedicated MPI class (RAII)
        if (size>1){
            double temp;
            MPI_Allreduce(&g_cur, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_cur=temp;
            MPI_Allreduce(&g_g1, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_g1=temp;
            MPI_Allreduce(&g_g2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_g2=temp;
            MPI_Allreduce(&g_x1, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_x1=temp;
            MPI_Allreduce(&g_x2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_x2=temp;
            MPI_Allreduce(&g_y1, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_y1=temp;
            MPI_Allreduce(&g_y2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_y2=temp;
        }
        double norm =1;
        if (g_cur>0){norm=1./g_cur;} // normalize with the weighting function (for beam it is the current)
        g_g1*=norm;
        g_g2*=norm;
        g_x1*=norm;
        g_x2*=norm;
        g_y1*=norm;
        g_y2*=norm;

        if (filter["energy"]){
            if (val.find("global/energy") != val.end()) {val["global/energy"][iz]=g_g1;}
            if (val.find("global/energyspread") != val.end()) {val["global/energyspread"][iz]=sqrt(fabs(g_g2-g_g1*g_g1));}
        }
        if (filter["spatial"]){
            if (val.find("global/xposition") != val.end()){val["global/xposition"][iz]=g_x1;}
            if (val.find("global/xsize") != val.end()){val["global/xsize"][iz]=sqrt( fabs(g_x2-g_x1*g_x1));}
            if (val.find("global/yposition") != val.end()){val["global/yposition"][iz]=g_y1;}
            if (val.find("global/xsize") != val.end()){val["global/ysize"][iz]=sqrt( fabs(g_y2-g_y1*g_y1));}
        }
    }

}

//---------------------------------------------------------
// diagnostic calculation - field

std::map<std::string,OutputInfo> DiagField::getTags(FilterDiagnostics & filter_in){

    tags.clear();
    filter.clear();

    tags["power"]={false,false,"W"};
//    tags["dgrid"]={true,true,"m"};

    return tags;
}

void DiagField::getValues(Field *field,std::map<std::string,std::vector<double> >&val, int iz){

    int ns = field->field.size();
    int is = field->first;

    for (auto const &slice :field->field){
        double power=0;

        int idx = iz*ns+is;         // index for saving the data
        if (val.find("power") != val.end()) { val["power"][idx] = power; }
        is++;
        is = is % ns;
    }
    return;
}