//
// Created by reiche on 12/11/23.
//

#ifndef GENESIS_1_3_VERSION4_IMPORTTRANSFORMATION_H
#define GENESIS_1_3_VERSION4_IMPORTTRANSFORMATION_H

#include <iostream>
#include <string>
#include <map>

class ImportTransformation {
public:
    ImportTransformation();
    virtual ~ImportTransformation();
    bool init(int, int, std::map<std::string,std::string> *);

private:
    static void usage();
};


#endif //GENESIS_1_3_VERSION4_IMPORTTRANSFORMATION_H
