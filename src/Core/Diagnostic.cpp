//
// Created by reiche on 06.01.22.
//

#include <iostream>
#include <complex>

#include "Diagnostic.h"
#include "Output.h"

//----------------------------------------------------------
// main class which manage all the calculation on the beam, field or undulator.
// It allocates the memory and define the interface to the output class so that the
// diagnostics classes are only concern with the actual calculation.

// routine to register the output parameters and to allocate memory
void Diagnostic::init(int rank, int size, int nz_in, int ns_in, int nfld,bool isTime, bool isScan, FilterDiagnostics &filter){

 //   FilterDiagnostics filter;
    nz = nz_in;
    ns = ns_in;
    time = isTime;
    scan = isScan;
    noff =rank*ns;
    ntotal = ns * size;
    val.clear();
    units.clear();
    single.clear();

    val.resize(3+nfld);   // global + lattice + beam + #fields
    units.resize(3+nfld);
    single.resize(3+nfld);
    zout.resize(nz);

    // loop through the different class for registration and allocate memory

    // undulator group (val.at(0) will be defined later in this->writeOutputFile

    // beam
    for (auto group : dbeam){
        for (auto const &tag : group->getTags(filter)){

            units.at(2)[tag.first] = tag.second.units;   // accumulate all units
            single.at(2)[tag.first] = tag.second.global;
            int size = ns * nz;
            if (tag.second.once) { size /= nz; }
            if (tag.second.global) {size /= ns; }
            val.at(2)[tag.first].resize(size);
        }
    }

    // field
    for (int ifld=0; ifld<nfld;ifld++){
        for (auto group : dfield){
            for (auto const &tag : group->getTags(filter)){
                units.at(3+ifld)[tag.first] = tag.second.units;   // accumulate all units
                single.at(3+ifld)[tag.first] = tag.second.global;   // accumulate all units

                int size = ns * nz;
                if (tag.second.once) { size /= nz; }
                if (tag.second.global) {size /= ns; }
                val.at(3+ifld)[tag.first].resize(size);
            }
        }
    }

    iz=0;                              // set the counter
}

void Diagnostic::addOutput(int groupID, std::string key, std::string unit, std::vector<double> &data){
    val.at(groupID)[key]=data;
    units.at(groupID)[key]=unit;
    single.at(groupID)[key]=true;
}

// adds some output and lfushes everything to a file
void Diagnostic::writeToOutputFile(std::string root, Beam *beam, vector<Field*> *field, Undulator *und)
{

    Output *out=new Output;
//    string file=root.append(".test");
    this->addOutput(0,"zplot","m", zout);
    this->addOutput(0,"z","m", und->z);
    this->addOutput(0,"dz","m", und->dz);
    this->addOutput(0,"aw"," ", und->aw);
    this->addOutput(0,"ax","m", und->ax);
    this->addOutput(0,"ay","m", und->ay);
    this->addOutput(0,"ku","m^{-1}", und->ku);
    this->addOutput(0,"kx"," ", und->kx);
    this->addOutput(0,"ky"," ", und->ky);
    this->addOutput(0,"qf","m^{-2}", und->qf);
    this->addOutput(0,"qx","m", und->qx);
    this->addOutput(0,"qy","m", und->qy);
    this->addOutput(0,"cx","rad", und->cx);
    this->addOutput(0,"cy","rad", und->cy);
    this->addOutput(0,"gradx","m^{-1}", und->gradx);
    this->addOutput(0,"grady","m^{-1}", und->grady);
    this->addOutput(0,"slippage"," ", und->slip);
    this->addOutput(0,"phaseshift","rad", und->phaseshift);
    this->addOutput(0,"chic_angle","degree", und->chic_angle);
    this->addOutput(0,"chic_lb","m", und->chic_lb);
    this->addOutput(0,"chic_ld","m", und->chic_ld);
    this->addOutput(0,"chic_lt","m", und->chic_lt);

    vector<double> global  (1);
    global[0]=und->getGammaRef();
    this->addOutput(1,"gamma0","mc^2",global);
    global[0]=beam->reflength;
    this->addOutput(1,"lambdaref","m",global);
    global[0]=beam->slicelength/beam->reflength;
    this->addOutput(1,"sample"," ",global);
    global[0]=beam->slicelength*ntotal;
    this->addOutput(1,"slen","m",global);
    global[0] = beam->one4one ? 1. : 0 ;
    this->addOutput(1,"one4one"," ",global);
    global[0] = time ? 1. : 0 ;
    this->addOutput(1,"time"," ",global);
    global[0] = scan ? 1. : 0 ;
    this->addOutput(1,"scan"," ",global);
    global.resize(ntotal);
    for (int i=0; i<ntotal; i++){
        global[i]=static_cast<double>(i)*beam->slicelength;
    }
    this->addOutput(1,"s","m",global);
    double e0=1023.842e-9/beam->reflength;
    double df=e0*beam->reflength/beam->slicelength/static_cast<double>(ntotal);
    if (ntotal == 1) {
        df=0;
    }
    e0 -=0.5*df*ntotal;
    for (int i=0; i<ntotal; i++){
        global[i]=e0+static_cast<double>(i)*df;
    }
    this->addOutput(1,"frequency","ev",global);

    out->open(root,noff,ns);
    out->writeMeta(und);
    out->writeGroup("Lattice",val[0], units[0],single[0]);
    out->writeGroup("Global",val[1], units[1],single[1]);
    out->writeGroup("Beam",val[2], units[2],single[2]);
    for (int i=3; i<val.size();i++){
        int h=field->at(i-3)->harm;
        if (h==1){
            out->writeGroup("Field",val[i], units[i],single[i]);
        } else {
            char buff[30];
            snprintf(buff, sizeof(buff), "Field%d", h);
            out->writeGroup(buff,val[i], units[i],single[i]);
        }
    }
    out->close();
    delete out;
    return;
}

// wrapper to do all the diagnostics calculation at a given integration step iz.
void Diagnostic::calc(Beam* beam,std::vector<Field*> *field,double z) {
    zout[iz]=z;
    for (auto group: dbeam){
        group->getValues(beam,val[2],iz);
    }

    for (int ifld=0; ifld < field->size(); ifld++){
        for (auto group: dfield) {
            group->getValues(field->at(ifld),val[3+ifld],iz);
        }
    }
    iz++;
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
            tags["Global/energy"] = {true, false, "mc^2"};
            tags["Global/energyspread"] = {true, false, "mc^2"};
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
            tags["Global/xposition"] = {true, false, "m"};
            tags["Global/xsize"] = {true, false, "m"};
            tags["Global/yposition"] = {true, false, "m"};
            tags["Global/ysize"] = {true, false, "m"};
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
        tags["efield"] = {false, false, "eV/m"}; // full field which changes particle energy
        tags["wakefield"] = {false, false, "eV/m"}; // effect from wakefields
        tags["LSCfield"] = {false, false, "eV/m"}; // effect from space charge field
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
//           complex<double> phasor = complex<double> (cos(par.theta),sin(par.theta));
//           complex<double> phasor_acc = phasor;
//           b[0] += phasor;
           for(int iharm=0; iharm<nharm;iharm++){
 //              phasor_acc *= phasor;
//               b[iharm]+=phasor_acc;
                b[iharm]+=complex<double> (cos((iharm+1)*par.theta),sin((iharm+1)*par.theta));
//                b[iharm]+=phasor*phasor;
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
            if (val.find("ysize") != val.end()){val["ysize"][idx]=sqrt( fabs(y2-y1*y1));}
            if (val.find("pxposition") != val.end()){val["pxposition"][idx]=px1;}
            if (val.find("pyposition") != val.end()){ val["pyposition"][idx]=py1;}
        }

        if (val.find("bunching") != val.end()) {val["bunching"][idx] =std::abs(b[0]);}
        if (val.find("bunchingphase") != val.end()) {val["bunchingphase"][idx] =atan2(b[0].imag(),b[0].real());}

        char buff[100];
        for (int iharm = 1; iharm < nharm; iharm++) {
            snprintf(buff, sizeof(buff), "bunching%d", iharm + 1);
            if (val.find(buff) != val.end()) { val[buff][idx] = std::abs(b[iharm]); }
            snprintf(buff, sizeof(buff), "bunchingphase%d", iharm + 1);
            if (val.find(buff) != val.end()) { val[buff][idx] = atan2(b[iharm].imag(), b[iharm].real()); }
        }
        if (filter["aux"]){
            if (val.find("efield") != val.end()) {val["efield"][idx]=beam->eloss[is]+beam->longESC[is];}
            if (val.find("wakefield") != val.end()) {val["wakefield"][idx]=beam->eloss[is];}
            if (val.find("LSCfield") != val.end()) {val["LSCfield"][idx]=beam->longESC[is];}
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
            if (val.find("alphax") != val.end()) { val["alphax"][is] = -(xpx-x1*px1)/ex;}
            if (val.find("alphay") != val.end()) { val["alphay"][is] = -(ypy-y1*py1)/ey; }
        }
        //-------------------------------------------------------------------------
        // gather moments for all slices with current as weighting factor
        if (global) {
            g_cur += beam->current[is];
            g_g1 += beam->current[is] * g1;
            g_g2 += beam->current[is] * g2;
            g_x1 += beam->current[is] * x1;
            g_x2 += beam->current[is] * x2;
            g_y1 += beam->current[is] * y1;
            g_y2 += beam->current[is] * y2;
        }
        // increment slice counter
        is++;
    }
    //-------------------------------------------------
    // complete gathering from all nodes if MPI size is larger than one
    if (global){
        int size = 1;
        if (!MPISingle) { MPI_Comm_size(MPI_COMM_WORLD, &size); }  // for future functionality to do scan functionality of time-dependent runs!
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
            if (val.find("Global/energy") != val.end()) {val["Global/energy"][iz]=g_g1;}
            if (val.find("Global/energyspread") != val.end()) {val["Global/energyspread"][iz]=sqrt(fabs(g_g2-g_g1*g_g1));}
        }
        if (filter["spatial"]){
            if (val.find("Global/xposition") != val.end()){val["Global/xposition"][iz]=g_x1;}
            if (val.find("Global/xsize") != val.end()){val["Global/xsize"][iz]=sqrt( fabs(g_x2-g_x1*g_x1));}
            if (val.find("Global/yposition") != val.end()){val["Global/yposition"][iz]=g_y1;}
            if (val.find("Global/xsize") != val.end()){val["Global/ysize"][iz]=sqrt( fabs(g_y2-g_y1*g_y1));}
        }
    }
}

//---------------------------------------------------------
// diagnostic calculation - field

FFTObj::FFTObj(int ngrid)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

	// this is severe error: print on *all* ranks
	if(ngrid<0) {
		cerr << "Error: enforcing ngrid>0" << endl;
		ngrid=1;
	}
	in_ = new complex<double>[ngrid * ngrid];
	out_ = new complex<double>[ngrid * ngrid];
	p_ = fftw_plan_dft_2d(ngrid, ngrid,
	    reinterpret_cast<fftw_complex *>(in_),
	    reinterpret_cast<fftw_complex *>(out_),
	    FFTW_FORWARD, FFTW_ESTIMATE /* FFTW_MEASURE */);
	ngrid_ = ngrid;
}

FFTObj::~FFTObj() {
	if(rank_==0) {
		cout << "~FFTObj" << endl;
	}

	fftw_destroy_plan(p_);
	delete [] in_;
	delete [] out_;
}

void DiagField::cleanup_FFT_resources(void)
{
	for(auto &[k,obj]: fftobj) {
		delete obj;
	}
	fftobj.clear();
}

int DiagField::obtain_FFT_resources(int ngrid, complex<double> **in, complex<double> **out, fftw_plan *pp)
{
	int rank=0;
	bool verbose=false;
	bool exists=false;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// see if the entry already exists
	auto end = fftobj.end();
	auto ele = fftobj.find(ngrid);
	exists = (ele!=end);

	if(!exists) {
		/* set up new entry */
		FFTObj *n = new FFTObj(ngrid);

		// FIXME: 'insert' also returns iterator to newly inserted element
		fftobj.insert({ngrid, n});
		ele = fftobj.find(ngrid);

		if(verbose && (rank==0)) {
			cout << "created FFT obj for ngrid="<<ngrid<<endl;
		}
	} else {
		if(verbose && (rank==0)) {
			cout << "getting existing FFT obj for ngrid="<<ngrid<<endl;
		}
	}

	// cout << "ngrid=" << ele->first << endl;
	*in = ele->second->in_;
	*out = ele->second->out_;
	*pp = ele->second->p_;

	return(0);
}


std::map<std::string,OutputInfo> DiagField::getTags(FilterDiagnostics & filter_in){

    tags.clear();
    filter.clear();
    global = filter_in.field.global;

    tags["power"]={false,false,"W"};
    if (global){
        tags["Global/energy"]={true,false,"J"};
    }

    if (filter_in.field.spatial) {
        filter["spatial"] = true;
        tags["xposition"] = {false, false, "m"};
        tags["xsize"] = {false, false, "m"};
        tags["yposition"] = {false, false, "m"};
        tags["ysize"] = {false, false, "m"};
        if (global) {
            tags["Global/xposition"] = {true, false, "m"};
            tags["Global/xsize"] = {true, false, "m"};
            tags["Global/yposition"] = {true, false, "m"};
            tags["Global/ysize"] = {true, false, "m"};
        }
    }
    if (filter_in.field.intensity) {
        filter["intensity"] = true;
        tags["intensity-nearfield"] = {false, false, "W/m^2"};
        tags["phase-nearfield"] = {false, false, "rad"};
        tags["intensity-farfield"] = {false, false, " "};
        tags["phase-farfield"] = {false, false, "rad"};
        if (global){
            tags["Global/intensity-nearfield"] = {true, false, "W/m^2"};
            tags["Global/intensity-farfield"] = {true, false, " "};
        }
    }
#ifdef FFTW
    if (filter_in.field.fft) {
        filter["fft"] = true;
        tags["xdivergence"] = {false, false, "rad"};
        tags["ydivergence"] = {false, false, "rad"};
        tags["xpointing"] = {false, false, "rad"};
        tags["ypointing"] = {false, false, "rad"};
        if (global) {
            tags["Global/xdivergence"] = {true, false, "rad"};
            tags["Global/ydivergence"] = {true, false, "rad"};
            tags["Global/xpointing"] = {true, false, "rad"};
            tags["Global/ypointing"] = {true, false, "rad"};
        }
    }
#endif
    // some basic parameter output
    tags["dgrid"]={true,true,"m"};
    tags["ngrid"]={true,true," "};
    return tags;
}

void DiagField::getValues(Field *field,std::map<std::string,std::vector<double> >&val, int iz){

    int ns = field->field.size();
    int is0 = 0;
    int ngrid = field->ngrid;

    double ks=4.*asin(1)/field->xlambda;
    double scl=field->dgrid*eev/ks;
    double scltheta=field->xlambda/ngrid/field->dgrid;
 //   std::cout  << "new: " << scltheta << " " << field->xlambda << " " << ngrid << std::endl;
    double shift=-0.5*static_cast<double> (ngrid-1);
//    double shift=-0.5*static_cast<double> (ngrid);

    // global variables
    double g_pow=0;
    double g_x1=0;
    double g_x2=0;
    double g_y1=0;
    double g_y2=0;
    double g_ff=0;
    double g_inten=0;

#ifdef FFTW
    double f_pow=0;
    double f_x1=0;
    double f_x2=0;
    double f_y1=0;
    double f_y2=0;

    complex<double> *in  = nullptr;
    complex<double> *out = nullptr;
    fftw_plan p;
    obtain_FFT_resources(ngrid, &in, &out, &p);
#endif


    for (auto const &slice :field->field) {
        int is = (ns + is0 - field->first) % ns;
        double power = 0;
        double x1 = 0;
        double x2 = 0;
        double y1 = 0;
        double y2 = 0;
        complex<double> loc;
        complex<double> ff = complex<double>(0, 0);
        for (int iy = 0; iy < ngrid; iy++) {
            double dy = static_cast<double>(iy) + shift;
            for (int ix = 0; ix < ngrid; ix++) {
                double dx = static_cast<double>(ix) + shift;
                int i = iy * ngrid + ix;
                loc = slice.at(i);
#ifdef FFTW
                in[i]=loc;   // field for the FFT
#endif
                double wei = loc.real() * loc.real() + loc.imag() * loc.imag();
                ff += loc;
                power += wei;
                x1 += dx * wei;
                x2 += dx * dx * wei;
                y1 += dy * wei;
                y2 += dy * dy * wei;
            }
        }

#ifdef FFTW
        double fpower=0;
        double fx1=0;
        double fx2=0;
        double fy1=0;
        double fy2=0;

        if (filter["fft"]){
            fftw_execute(p);
            for (int iy=0;iy<ngrid;iy++){
                double dy=static_cast<double>(iy)+shift;
                for (int ix=0;ix<ngrid;ix++){
                    double dx=static_cast<double>(ix)+shift;
                    int iiy=(iy+(ngrid+1)/2) % ngrid;
                    int iix=(ix+(ngrid+1)/2) % ngrid;
                    int ii=iiy*ngrid+iix;
                    loc=out[ii];
                    double wei=loc.real()*loc.real()+loc.imag()*loc.imag();
                    fpower+=wei;
                    fx1+=dx*wei;
                    fx2+=dx*dx*wei;
                    fy1+=dy*wei;
                    fy2+=dy*dy*wei;
                }
            }
        }
#endif

        if (global) {
            g_pow += power;
            g_x1 += x1;
            g_x2 += x2;
            g_y1 += y1;
            g_y2 += y2;
        }

        if (power > 0){
            x1/=power;
            x2/=power;
            y1/=power;
            y2/=power;
        }
        x2=sqrt(fabs(x2-x1*x1));
        y2=sqrt(fabs(y2-y1*y1));

#ifdef FFTW
        if (global){
            f_pow+=fpower;
            f_x1+=fx1;
            f_x2+=fx2;
            f_y1+=fy1;
            f_y2+=fy2;
        }

        if (fpower >0){
            fx1/=fpower;
            fx2/=fpower;
            fy1/=fpower;
            fy2/=fpower;
        }
 //       std::cout << "New: fpower: " << fpower << " fx1: " << fx1 <<  " scltheta: " << scltheta << std::endl;
	fx2=sqrt(fabs(fx2-fx1*fx1))*scltheta;
	fy2=sqrt(fabs(fy2-fy1*fy1))*scltheta;
        fx1*=scltheta;
        fy1*=scltheta;
#endif


        int i=(ngrid*ngrid-1)/2;
        loc=slice.at(i);
        double inten=loc.real()*loc.real()+loc.imag()*loc.imag();
        double intenphi=atan2(loc.imag(),loc.real());
        double farfield=ff.real()*ff.real()+ff.imag()*ff.imag();
        double farfieldphi =atan2(ff.imag(),ff.real());

        if (global){
            g_ff+=farfield;
            g_inten+=inten;
        }

        // scale to physical dimensions
        power*=scl*scl/vacimp; // scale to W
        x1*=field->dgrid;
        x2*=field->dgrid;
        y1*=field->dgrid;
        y2*=field->dgrid;
        inten*=eev*eev/ks/ks/vacimp;  // scale to W/m^2



        // save the data into the provided arrays
        int idx = iz*ns+is;         // index for saving the data
        if (val.find("power") != val.end()) { val["power"][idx] = power; }
        if (filter["spatial"]){
            if (val.find("xposition") != val.end()){val["xposition"][idx]=x1;}
            if (val.find("xsize") != val.end()){val["xsize"][idx]=x2;}
            if (val.find("yposition") != val.end()){val["yposition"][idx]=y1;}
            if (val.find("ysize") != val.end()){val["ysize"][idx]=y2;}
        }
        if (filter["intensity"]){
            if (val.find("intensity-nearfield") != val.end()){val["intensity-nearfield"][idx]=inten;}
            if (val.find("phase-nearfield") != val.end()){val["phase-nearfield"][idx]=intenphi;}
            if (val.find("intensity-farfield") != val.end()){val["intensity-farfield"][idx]=farfield;}
            if (val.find("phase-farfield") != val.end()){val["phase-farfield"][idx]=farfieldphi;}
        }
#ifdef FFTW
        if (filter["fft"]){
            if (val.find("xpointing") != val.end()){val["xpointing"][idx]=fx1;}
            if (val.find("xdivergence") != val.end()){val["xdivergence"][idx]=fx2;}
            if (val.find("ypointing") != val.end()){val["ypointing"][idx]=fy1;}
            if (val.find("ydivergence") != val.end()){val["ydivergence"][idx]=fy2;}
        }
#endif
        // increase the slice counter
        is0++;
    }
    // global variables
    //-------------------------------------------------
    // complete gathering from all nodes if MPI size is larger than one
    if (global) {
        int size = 1;
        if (!MPISingle) {
            MPI_Comm_size(MPI_COMM_WORLD, &size);
        }  // for future functionality to do scan functionalility of time-dependent runs!
        // this will be better controlled in the future with a dedicated MPI class (RAII)
        if (size > 1) {
            double temp;
            MPI_Allreduce(&g_pow, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_pow = temp;
            MPI_Allreduce(&g_x1, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_x1 = temp;
            MPI_Allreduce(&g_x2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_x2 = temp;
            MPI_Allreduce(&g_y1, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_y1 = temp;
            MPI_Allreduce(&g_y2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_y2 = temp;
            MPI_Allreduce(&g_ff, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_ff = temp;
            MPI_Allreduce(&g_inten, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            g_inten = temp;
#
#ifdef FFTW
            MPI_Allreduce(&f_pow, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            f_pow = temp;
            MPI_Allreduce(&f_x1, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            f_x1 = temp;
            MPI_Allreduce(&f_x2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            f_x2 = temp;
            MPI_Allreduce(&f_y1, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            f_y1 = temp;
            MPI_Allreduce(&f_y2, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            f_y2 = temp;
#endif
        }
        double norm = 1;
        if (g_pow > 0) { norm = 1. / g_pow; } // normalize with the weighting function (for beam it is the current)
        g_x1 *= norm;
        g_x2 *= norm;
        g_y1 *= norm;
        g_y2 *= norm;
        g_x2=sqrt(fabs(g_x2-g_x1*g_x1))*field->dgrid;
        g_y2=sqrt(fabs(g_y2-g_y1*g_y1))*field->dgrid;
        g_x1*=field->dgrid;
        g_y1*=field->dgrid;

#ifdef FFTW
        norm = 1;
        if (f_pow > 0) { norm = 1. / f_pow; } // normalize with the weighting function (for beam it is the current)
        f_x1 *= norm;
        f_x2 *= norm;
        f_y1 *= norm;
        f_y2 *= norm;
        f_x2=sqrt(fabs(f_x2-f_x1*f_x1))*scltheta;
        f_y2=sqrt(fabs(f_y2-f_y1*f_y1))*scltheta;
        f_x1*=scltheta;
        f_y1*=scltheta;
#endif
        g_pow *=scl*scl/vacimp*field->slicelength/3e8;
        norm = 1./static_cast<double>(size*field->field.size());
        g_inten *= norm*eev*eev/ks/ks/vacimp;
        g_ff*=norm;
        if (val.find("Global/energy") != val.end()) { val["Global/energy"][iz] = g_pow; }
        if (filter["spatial"]) {
            if (val.find("Global/xposition") != val.end()) { val["Global/xposition"][iz] = g_x1; }
            if (val.find("Global/xsize") != val.end()) { val["Global/xsize"][iz] = g_x2; }
            if (val.find("Global/yposition") != val.end()) { val["Global/yposition"][iz] = g_y1; }
            if (val.find("Global/xsize") != val.end()) { val["Global/ysize"][iz] = g_y2; }
        }

        if (filter["intensity"]) {
            if (val.find("Global/intensity-nearfield") !=
                val.end()) { val["Global/intensity-nearfield"][iz] = g_inten; }
            if (val.find("Global/intensity-farfield") != val.end()) { val["Global/intensity-farfield"][iz] = g_ff; }
        }
#ifdef FFTW
        if (filter["fft"]) {
            if (val.find("Global/xpointing") != val.end()) { val["Global/xpointing"][iz] = f_x1; }
            if (val.find("Global/xdivergence") != val.end()) { val["Global/xdivergence"][iz] = f_x2;}
            if (val.find("Global/ypointing") != val.end()) { val["Global/ypointing"][iz] = f_y1; }
            if (val.find("Global/ydivergence") != val.end()) {val["Global/ydivergence"][iz] = f_y2;}
        }
#endif
    }
    if ((iz == 0) && (val.find("dgrid") != val.end())) { val["dgrid"][0] = field->dgrid; }
    if ((iz == 0) && (val.find("ngrid") != val.end())) { val["ngrid"][0] = static_cast<double> (ngrid); }

    return;
}

