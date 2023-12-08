//
// Created by reiche on 12/8/23.
//

#include "readMapHDF5.h"

ReadMapHDF5::ReadMapHDF5()
{
    isOpen=false;
    nwork=-1;
}

ReadMapHDF5::~ReadMapHDF5()
{
    if (nwork>0){ delete [] work; }
}

void ReadMapHDF5::close(){
    if (isOpen){ H5Fclose(fid); }
}


bool ReadMapHDF5::read(int rank, std::string file, std::string dset_rvec, std::string dset_rmat){

    if ((fid=H5Fopen(file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT)) == H5I_INVALID_HID) {
        if(0==rank) {
            std::cout << "*** Error: unable to open file " << file << endl;
        }
        return(false);
    }

    auto shapevec=getFullDatasetSize(fid, dset_rvec.c_str());
    if (shapevec[0] <0) {
        if(0==rank) {
            std::cout << "*** Error: unable to open dataset " << dset_rvec << " in " << file << endl;
        }
        return(false);
    }





    auto shapemat=getFullDatasetSize(fid, dset_rmat.c_str());
    if (shapemat[0] <0){
        if(0==rank) {
            std::cout << "*** Error: unable to open dataset " << dset_rmat << " in " << file << endl;
        }
        return(false);
    }
    std::cout << dset_rvec << " - size 1st dim " << shapevec.at(0) << endl;
    std::cout << dset_rmat << " - size 1st dim " << shapemat.at(0) << endl;

    /*
 *
    if (nsize>nwork){ // allocate extra work array to hold field
        if (nwork>0) {delete [] work;}
        nwork=nsize;
        if (nwork<1){
            nwork=1;    // catch the error that the first slice can have no particles in one4one simulations
        }
        work=new double [nwork];
    }


    // allocate data in the beam record
    slice->resize(nsize);

    // get current
    sprintf(name,"slice%6.6d/current",islice);
    readDataDouble(fid,name,current,1);
    *current*=wei;

    sprintf(name,"slice%6.6d/gamma",islice);
    readDataDouble(fid,name,work,nsize);
    for (int i=0;i<nsize;i++){
        slice->at(i).gamma=work[i];
    }

    sprintf(name,"slice%6.6d/theta",islice);
    readDataDouble(fid,name,work,nsize);
    for (int i=0;i<nsize;i++){
        slice->at(i).theta=work[i];
    }

    sprintf(name,"slice%6.6d/x",islice);
    readDataDouble(fid,name,work,nsize);
    for (int i=0;i<nsize;i++){
        slice->at(i).x=work[i]*wei;
    }

    sprintf(name,"slice%6.6d/y",islice);
    readDataDouble(fid,name,work,nsize);
    for (int i=0;i<nsize;i++){
        slice->at(i).y=work[i]*wei;
    }

    sprintf(name,"slice%6.6d/px",islice);
    readDataDouble(fid,name,work,nsize);
    for (int i=0;i<nsize;i++){
        slice->at(i).px=work[i]*wei;
    }

    sprintf(name,"slice%6.6d/py",islice);
    readDataDouble(fid,name,work,nsize);
    for (int i=0;i<nsize;i++){
        slice->at(i).py=work[i]*wei;
    }
*/
    return true;
}

