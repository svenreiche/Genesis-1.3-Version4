//
// Created by reiche on 12/8/23.
//

#ifndef GENESIS_1_3_VERSION4_READMAPHDF5_H
#define GENESIS_1_3_VERSION4_READMAPHDF5_H



#include <hdf5.h>
#include "HDF5base.h"
#include <string>

class ReadMapHDF5 : public HDF5Base {
public:
    ReadMapHDF5();
    virtual ~ReadMapHDF5();
    bool read(int, std::string, std::string, std::string);
    void close();

private:
    hid_t fid;
    int nwork;
    double *work;
    bool isOpen;
};


#endif //GENESIS_1_3_VERSION4_READMAPHDF5_H
