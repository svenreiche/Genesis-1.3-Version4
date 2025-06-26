//
// Created by reiche on 12/11/23.
//


#include "ImportTransformation.h"
#include "readMapHDF5.h"



ImportTransformation::ImportTransformation()= default;
ImportTransformation::~ImportTransformation()= default;

void ImportTransformation::usage(){
    std::cout << "List of keywords for IMPORTTRANSFORMATION" << std::endl;
    std::cout << "&importtransformation" << std::endl;
    std::cout << " string file = <empty>" << std::endl;
    std::cout << " string vector =<empty>" << std::endl;
    std::cout << " string matrix =<empty>" << std::endl;
    std::cout << " double slen = 0" << std::endl;
    std::cout << "&end" << std::endl << std::endl;
}


bool ImportTransformation::init(int rank_in, std::map<std::string,std::string> *arg,Beam *beam,Setup *setup) {

    rank = rank_in;
    if (beam->beam.empty()) {
        if (rank == 0) { cout << "*** Error: No beam defined yet for applying transformation" << endl; }
        return false;
    }

    lambda = setup->getReferenceLength();   // reference length for theta
    kr = 4.*asin(1)/lambda;
    gamma = setup->getReferenceEnergy();           // get default energy from setup input deck
    nslice = beam->beam.size();
    dslice = beam->slicelength;
    slen = 0;

    auto end = arg->end();
    char *ptr;
    std::string file, vector, matrix;


    if (arg->find("file") != end) {
        file = arg->at("file");
        arg->erase(arg->find("file"));
    }
    if (arg->find("vector") != end) {
        vector = arg->at("vector");
        arg->erase(arg->find("vector"));
    }
    if (arg->find("matrix") != end) {
        matrix = arg->at("matrix");
        arg->erase(arg->find("matrix"));
    }
    if (arg->find("slen") != end) {
        slen = strtod(arg->at("slen").c_str(), &ptr);
        arg->erase(arg->find("slen"));
    }

    if (!arg->empty()) {
        if (rank == 0) {
            std::cout << "*** Error: Unknown elements in &importtransformation" << std::endl;
            this->usage();
        }
        return false;
    }

    if (rank == 0) {
        cout << "Applying particle transformation from file: " << file << " ..." << endl;
    }


    // read the hdf5 file
    ReadMapHDF5 import;
    if (!import.open(rank, file)) {
        return false;
    }
    if (!import.readDataset(vector, rvec, true)) {
        return false;
    }
    if (!import.readDataset(matrix, rmat, false)) {
        return false;
    }
    import.close();
    // get count of supplied vectors or matrices, respectively
    nvec = rvec.size() / 6;
    nmat = rmat.size() / 36;

    // apply matrix transformation
    if (nmat > 0) {
        this->applyMatrix(beam->beam);
    }
    // apply centroid shift
    if (nvec > 0) {
        this->applyVector(beam->beam);
    }
    return true;
}

void ImportTransformation::applyMatrix(vector<vector<Particle>> &beam) {

    bool interpolate = false;
    if ((nmat > 1) && (slen>0)) {interpolate = true;}
    std::vector<double> r0 = {0,0,0,0,0,0};
    std::vector<double> r1 = {0,0,0,0,0,0};

    for (int islice = 0; islice<nslice; islice++){
        unsigned int ioff = 0;  // default first entry
        bool single = true;
        double wei = 1;

        if (interpolate){
            auto s = static_cast<double>(islice + rank*nslice)*dslice;  // current slice position
            if (s > static_cast<double>(nvec)*slen){
                ioff = nvec-1;
            } else {
                single = false;
                auto rat = s/slen;
                wei = 1-(rat-floor(rat));
                ioff = static_cast<unsigned long>(floor(rat));
            }
        }
        auto npart = beam.at(islice).size();
        for (int ipart = 0; ipart <npart; ipart++){
            r0.at(0)=beam.at(islice).at(ipart).x;
            r0.at(1)=beam.at(islice).at(ipart).px/gamma;
            r0.at(2)=beam.at(islice).at(ipart).y;
            r0.at(3)=beam.at(islice).at(ipart).py/gamma;
            r0.at(4)=beam.at(islice).at(ipart).theta/kr;
            r0.at(5)=beam.at(islice).at(ipart).gamma/gamma-1;

            for (int j=0;j<6;j++){
                r1.at(j) = 0;
                for (int l=0;l<6;l++){
                    r1.at(j)+=wei*rmat.at(ioff*36+j*6+l)*r0.at(l);
                }
            }

            if (!single) {
                for (int j = 0; j < 6; j++) {
                    for (int l = 0; l < 6; l++) {
                        r1.at(j) += (1 - wei) * rmat.at(ioff * 36 + 36 + j * 6 + l) * r0.at(l);
                    }
                }
            }

            beam.at(islice).at(ipart).x = r1.at(0);
            beam.at(islice).at(ipart).px = r1.at(1) * gamma;
            beam.at(islice).at(ipart).y = r1.at(2);
            beam.at(islice).at(ipart).py = r1.at(3) * gamma;
            beam.at(islice).at(ipart).theta = r1.at(4) * kr;
            beam.at(islice).at(ipart).gamma = 1+r1.at(5) * gamma;
        }
    }
}

void ImportTransformation::applyVector(vector<vector<Particle>> &beam) {

    bool interpolate = false;
    if ((nvec > 1) && (slen>0)) {interpolate = true;}

    for (int islice = 0; islice<nslice; islice++){
        unsigned int ioff = 0;  // default first entry
        bool single = true;
        double wei = 1;

        if (interpolate){
            auto s = static_cast<double>(islice + rank*nslice)*dslice;  // current slice position
            if (s > static_cast<double>(nvec)*slen){
                ioff = nvec-1;
            } else {
                single = false;
                auto rat = s/slen;
                wei = 1-(rat-floor(rat));
                ioff = static_cast<unsigned long>(floor(rat));
            }
        }
        auto npart = beam.at(islice).size();
        for (int ipart = 0; ipart <npart; ipart++){
            beam.at(islice).at(ipart).x  += wei*rvec.at(ioff*6);
            beam.at(islice).at(ipart).px += wei*rvec.at(ioff*6+1)*gamma;
            beam.at(islice).at(ipart).y  += wei*rvec.at(ioff*6+2);
            beam.at(islice).at(ipart).py += wei*rvec.at(ioff*6+3)*gamma;
            beam.at(islice).at(ipart).theta += wei*rvec.at(ioff*6+4)*kr;
            beam.at(islice).at(ipart).gamma += wei*rvec.at(ioff*6+5)*gamma;
        }
        if (!single) {
            for (int ipart = 0; ipart < npart; ipart++) {
                beam.at(islice).at(ipart).x += (1-wei) * rvec.at(ioff * 6+6);
                beam.at(islice).at(ipart).px += (1-wei) * rvec.at(ioff * 6 + 7) * gamma;
                beam.at(islice).at(ipart).y += (1-wei) * rvec.at(ioff * 6 + 8);
                beam.at(islice).at(ipart).py += (1-wei) * rvec.at(ioff * 6 + 9) * gamma;
                beam.at(islice).at(ipart).theta += (1-wei) * rvec.at(ioff * 6 + 10) * kr;
                beam.at(islice).at(ipart).gamma += (1-wei) * rvec.at(ioff * 6 + 11) * gamma;
            }
        }
    }
}


