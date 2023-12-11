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
    std::string file,vector,matrix;
    double slen = 0;

    if (arg->find("file")!=end    )  {file=arg->at("file"); arg->erase(arg->find("file"));}
    if (arg->find("vector")!=end    ){vector=arg->at("vector"); arg->erase(arg->find("vector"));}
    if (arg->find("matrix")!=end    ){matrix=arg->at("matrix"); arg->erase(arg->find("matrix"));}
    if (arg->find("slen")!=end    )  {slen=atof(arg->at("slen").c_str()); arg->erase(arg->find("slen"));}

    if (!arg->empty()){
        if (rank==0){ std::cout << "*** Error: Unknown elements in &importtransformation" << std::endl; this->usage();}
        return false;
    }

    if (rank==0){
        cout << "Importing particle transformation from file: " << file << " ..." << endl;
    }

    ReadMapHDF5 import;
    auto status = import.read(rank,file,vector,matrix);
    if (!status){
        return status;
    }
    return true;
}

