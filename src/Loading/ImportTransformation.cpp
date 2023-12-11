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


bool ImportTransformation::init(int rank, int size, std::map<std::string,std::string> *arg,Beam *beam,Setup *setup)
{

    if (beam->beam.empty()){
        if (rank==0) {cout << "*** Error: No beam defined yet for applying transformation" << endl; }
        return false;
    }

    double lambda=setup->getReferenceLength();   // reference length for theta
    double gamma=setup->getReferenceEnergy();           // get default energy from setup input deck

    auto end=arg->end();
    char *ptr;
    std::string file,vector,matrix;
    double slen = 0;

    if (arg->find("file")!=end    )  {file=arg->at("file"); arg->erase(arg->find("file"));}
    if (arg->find("vector")!=end    ){vector=arg->at("vector"); arg->erase(arg->find("vector"));}
    if (arg->find("matrix")!=end    ){matrix=arg->at("matrix"); arg->erase(arg->find("matrix"));}
    if (arg->find("slen")!=end    )  {slen=strtod(arg->at("slen").c_str(),&ptr); arg->erase(arg->find("slen"));}

    if (!arg->empty()){
        if (rank==0){ std::cout << "*** Error: Unknown elements in &importtransformation" << std::endl; this->usage();}
        return false;
    }

    if (rank==0){
        cout << "Importing particle transformation from file: " << file << " ..." << endl;
    }

    std::cout << "slen: " << slen << std::endl;
    std::vector<double> rvec,rmat;

    // read the hdf5 file
    ReadMapHDF5 import;
    if (!import.open(rank,file)) {
        return false;
    }
    if (!import.readDataset(vector, rvec, true)) {
        return false;
    }
    if (!import.readDataset(matrix, rmat, false)) {
        return false;
    }
    import.close();

    // transform the electron distribution
    auto nvec = rvec.size()/6;
    auto nmat = rmat.size()/36;

    auto nslice = beam->beam.size();
    auto dslice = beam->slicelength;


    bool interpolate = false;
    if ((nvec > 1) && (slen>0)) {interpolate = true;}
    int skip = nvec > 0 ? 1 : 0;  // check if a transformation vector has beed defined at all
    double scl = 4.*asin(1)/lambda;
    for (int islice = 0; islice<nslice*skip; islice++){
        auto s = static_cast<double>(islice + rank*nslice)*dslice;  // current slice position
        unsigned int ioff = 0;
        if (interpolate){
            ioff = 0;
        }
        auto npart = beam->beam.at(islice).size();
        for (int ipart = 0; ipart <npart; ipart++){
            beam->beam.at(islice).at(ipart).x  += rvec.at(ioff*6);
            beam->beam.at(islice).at(ipart).px += rvec.at(ioff*6+1)*gamma;
            beam->beam.at(islice).at(ipart).y  += rvec.at(ioff*6+2);
            beam->beam.at(islice).at(ipart).py += rvec.at(ioff*6+3)*gamma;
            beam->beam.at(islice).at(ipart).theta += rvec.at(ioff*6+4)*scl;
            beam->beam.at(islice).at(ipart).gamma += rvec.at(ioff*+5)*gamma;
        }
    }

    return true;
}

