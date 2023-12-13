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
    bool readDataset(const std::string &name, std::vector<double> &data, bool isVector);
    bool open(unsigned long rank_in, const string &file_in);
    void close() const;

private:
    hid_t fid {-1};
    std::string file;
    unsigned long rank {0};
    int nwork {0};
    int nsize {0};
    bool isOpen {false};
    double *work {nullptr};

    // member functions
    bool reportShape(const std::string &dset, bool isVector);
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
