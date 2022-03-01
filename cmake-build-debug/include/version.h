//
// Created by reiche on 05.01.22.
//

#ifndef GENESIS_1_3_VERSION4_VERSION_H
#define GENESIS_1_3_VERSION4_VERSION_H

#include <string>

class VersionInfo {
int vmajor = 4;
int vminor = 6;
int vrev = 1;
bool vbeta = true;


public:
    VersionInfo(){};
    ~VersionInfo(){};

    int Major(void) {return vmajor; };
    int Minor(void) {return vminor; };
    int Rev(void) {return vrev;};
    bool isBeta(void) {return vbeta;};
    const char *Build(void) {return "Compiled by reiche at 2022-03-01 10:33:53 [UTC] from Git Commit ID: a3c3f23dcf3a17f2031dc2c5a57d84915a56d824";};
};

#endif //GENESIS_1_3_VERSION4_VERSION_H
