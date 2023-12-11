//
// Created by reiche on 12/8/23.
//

#include "readMapHDF5.h"

ReadMapHDF5::ReadMapHDF5()
= default;

ReadMapHDF5::~ReadMapHDF5()
{
    if (nwork>0){ delete [] work; }
}

void ReadMapHDF5::close() const{
    if (isOpen){ H5Fclose(fid); }
}

bool ReadMapHDF5::open(int rank_in, const std::string& file_in) {
    rank = rank_in;
    file = file_in;
    if ((fid = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) == H5I_INVALID_HID) {
        if (rank == 0) {
            std::cout << "*** Error: unable to open file " << file << endl;
        }
        return (false);
    }
    isOpen = true;
    return true;
}

bool ReadMapHDF5::readDataset(const std::string& dset, std::vector<double> &data, bool isVector) {

    if (dset.empty()){
        data.resize(0);
        return true;
    }

    std::vector<int> shape = {0, 0, 0};
    bool status = getFullDatasetSize(fid, dset.c_str(), shape);
    if (!status) {
        if (rank == 0) {
            std::cout << "*** Error: unable to open dataset " << dset << " in " << file << std::endl;
        }
        return status;
    }

    // check for matrix
    if (!isVector) {
        if (shape[2] == 0) {
            shape[2] = shape[1];
            shape[1] = shape[0];
            shape[0] = 1;
        }
    } else { // or vector
        shape[2] = 0;
        if (shape[1] == 0) {
            shape[1] = shape[0];
            shape[0] = 1;
        }
    }

    // check if all elements were defined
    if (shape[1] != 6) {
        return (this->reportShape(dset, isVector));
    }
    if ((!isVector) && (shape[2] != 6)) {
        return (this->reportShape(dset, isVector));
    }

    nsize = 6*shape[0];
    if (!isVector) { nsize *=6;}
    data.resize(nsize);
    readDataDouble(fid, const_cast<char *>(dset.c_str()), &data[0], nsize);
    return true;
}

