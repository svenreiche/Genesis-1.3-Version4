//
// Created by reiche on 12/11/23.
//

#ifndef GENESIS_1_3_VERSION4_IMPORTTRANSFORMATION_H
#define GENESIS_1_3_VERSION4_IMPORTTRANSFORMATION_H

#include <iostream>
#include <string>
#include <map>
#include <cstdlib>

#include "Beam.h"
#include "Setup.h"

class ImportTransformation {
public:
    ImportTransformation();
    virtual ~ImportTransformation();
    bool init(int, std::map<std::string,std::string> *,Beam *,Setup *);

private:
    static void usage();
    void applyVector(vector<vector<Particle>> &beam);
    void applyMatrix(vector<vector<Particle>> &beam);

    double lambda{0},kr{0},gamma{0},dslice{0}, slen{0};
    unsigned long nslice{0};
    std::vector<double> rvec,rmat;
    unsigned long nvec{0},nmat{0},rank{0};

};


#endif //GENESIS_1_3_VERSION4_IMPORTTRANSFORMATION_H
