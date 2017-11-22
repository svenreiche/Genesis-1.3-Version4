
#include "writeFieldHDF5.h"
#include "MPEProfiling.h"

// constructor destructor
WriteFieldHDF5::WriteFieldHDF5()
{
}

WriteFieldHDF5::~WriteFieldHDF5()
{
}

void WriteFieldHDF5::write(string fileroot, vector<Field *> *field){

  string file;
  MPI::Status status;

  size=MPI::COMM_WORLD.Get_size(); // get size of cluster
  rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node

  int ntotal=0;
  int nslice=field->at(0)->field.size();
  MPI::COMM_WORLD.Reduce(&nslice,&ntotal,1,MPI::INT,MPI::SUM,0);


  mpe.logIO(false,true,"Write Field Dump to File");  

  for (int i=0; i<field->size();i++){
    int harm=field->at(i)->harm;
    char charm[10];
    sprintf(charm,".h%d",harm);
    if (harm==1){
      file=fileroot;
    } else {
      file=fileroot+string(charm);
    }
    this->writeMain(file,field->at(i),ntotal);
  }
  mpe.logIO(true,true,"Write Field Dump to File");  
  return;
}





void WriteFieldHDF5::writeMain(string fileroot, Field *field,int ntotal){


  int slicecount=0;
  int tag=1;
  MPI::Status status;
  char filename[100];



  sprintf(filename,"%s.fld.h5",fileroot.c_str());

  if (rank==0){

    hid_t fid=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT); 
    this->writeGlobal(fid,field->xlambda,field->slicelength,field->s0,field->dgrid,ntotal);
    slicecount+=this->writeSlice(fid, field, slicecount);
    H5Fclose(fid);

    if (size > 1){
      	MPI::COMM_WORLD.Send(&slicecount, 1, MPI::INT, rank+1, tag);
    }

  } else {

    MPI::COMM_WORLD.Recv(&slicecount, 1, MPI::INT, rank-1, tag,status);

    hid_t fid=H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);   
    slicecount+=this->writeSlice(fid, field, slicecount);
    H5Fclose(fid);

    if (rank < (size-1)){
      	MPI::COMM_WORLD.Send(&slicecount, 1, MPI::INT, rank+1, tag);
    }
  }


  return;
}

void WriteFieldHDF5::writeGlobal(hid_t fid, double reflen, double slicelen, double s0, double dx, int count)
{

  this->writeDataDouble(fid, "wavelength", &reflen, 1);
  this->writeDataDouble(fid, "slicespacing", &slicelen, 1);
  this->writeDataDouble(fid, "refposition", &s0, 1);
  this->writeDataDouble(fid, "gridsize", &dx, 1);
  this->writeDataInt(fid, "slicecount", &count, 1);
  
  return;
}




int WriteFieldHDF5::writeSlice(hid_t fid, Field *field, int off)
{
  int nslice=field->field.size();
  
  int ngrid=field->ngrid;
  vector<double> work;
  work.resize(ngrid*ngrid);

  for (int islice=0;islice<nslice;islice++){
 
      char slicename[20];
      sprintf(slicename,"slice%6.6d",islice+off+1);
      hid_t gid=H5Gcreate(fid,slicename,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);


      for (int i=0; i<ngrid*ngrid;i++){ work[i]=field->field.at(islice).at(i).real();}  
      this->writeDataDouble(gid, "field-real", &work[0], ngrid*ngrid);

      for (int i=0; i<ngrid*ngrid;i++){ work[i]=field->field.at(islice).at(i).imag();}  
      this->writeDataDouble(gid, "field-imag", &work[0], ngrid*ngrid);

      H5Gclose(gid);
       

  }
  return nslice;
}


