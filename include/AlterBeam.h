//
// Created by reiche on 4/17/25.
//

#ifndef GENESIS_1_3_VERSION4_ALTERBEAM_H
#define GENESIS_1_3_VERSION4_ALTERBEAM_H

#include <iostream>
#include <map>

#include "StringProcessing.h"
#include "Setup.h"
#include "GenTime.h"
#include "GenProfile.h"
#include "Beam.h"

class AlterBeam : public StringProcessing{
public:
    AlterBeam();
    virtual ~AlterBeam();

    bool init(int rank, int size, std::map<std::string,std::string>* arg, Beam *beam, Setup *setup, Time *time, Profile *prof);


private:
    static void usage();
};


#endif //GENESIS_1_3_VERSION4_ALTERBEAM_H
