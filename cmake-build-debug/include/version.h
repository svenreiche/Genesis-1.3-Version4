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
    const char *Build(void) {return "Compiled by reiche at 2022-03-04 15:14:46 [UTC] from Git Commit ID: 3dd7c2ae2ba4ff5e7ecf953c1d00919264fe0f38";};
};

#endif //GENESIS_1_3_VERSION4_VERSION_H
