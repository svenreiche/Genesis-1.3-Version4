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
    const char *Build(void) {return "Compiled by reiche at 2022-02-18 15:47:41 [UTC] from Git Commit ID: 3062d884629ff7fe960db6e391d5c13df52a56d3";};
};

#endif //GENESIS_1_3_VERSION4_VERSION_H
