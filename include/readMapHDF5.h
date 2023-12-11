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
    ~ReadMapHDF5() override;
    bool read(int, const std::string&, const std::string&, const std::string&);

private:
    hid_t fid {-1};
    std::string file;
    int rank {0};
    int nwork {0};
    double *work {nullptr};

    // member functions
    bool reportShape(const std::string &dset, bool isVector);
    bool readVector(const string &dset, bool isVector);
};

inline bool ReadMapHDF5::reportShape(const std::string& dset, bool isVector){
    if(rank == 0) {
        std::cout << "*** Error: invalid shape of " << dset << " in " << file << std::endl;
        std::cout << "           Must be " ;
        if (isVector){
            std::cout << "(6) or (n,6)" << std::endl;
        } else {
            std::cout << "(6,6) or (n,6,6)" << std::endl;
        }
    }
    return false;
}
#endif //GENESIS_1_3_VERSION4_READMAPHDF5_H
