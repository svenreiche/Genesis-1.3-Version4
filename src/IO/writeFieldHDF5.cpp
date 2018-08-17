
#include "writeFieldHDF5.h"

extern bool MPISingle;

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
  if (MPISingle){
    size=1;
    rank=0;
  }

  for (int i=0; i<field->size();i++){
    int harm=field->at(i)->harm;
    char charm[10];
    sprintf(charm,".h%d",harm);
    if (harm==1){
      file=fileroot;
    } else {
      file=fileroot+string(charm);
    }
    this->writeMain(file,field->at(i));
  }

  return;
}





void WriteFieldHDF5::writeMain(string fileroot, Field *field){


 

  char filename[100];
  sprintf(filename,"%s.fld.h5",fileroot.c_str()); 
  if (rank == 0) { cout << "Writing field distribution to file: " <<filename << " ..." << endl;} 

  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  if (size>1){
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  fid=H5Fcreate(filename,H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
  H5Pclose(pid);

  s0=rank;
  int ntotal=size*field->field.size();

  // write global data
  this->writeGlobal(field->xlambda,field->slicelength,field->s0,field->dgrid,field->ngrid,ntotal);

  // loop through slices
  
  int smin=rank*field->field.size();
  int smax=smin+field->field.size();

  int ngrid=field->ngrid;
  vector<double> work;
  work.resize(ngrid*ngrid);

  double ks=4.*asin(1)/field->xlambda;
  double scl=field->dgrid*eev/ks/sqrt(vacimp);


  for (int i=0; i<ntotal;i++){
    s0=-1;
    char name[16];
    sprintf(name,"slice%6.6d",i+1);
    hid_t gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    if ((i>=smin) && (i<smax)){
      s0=0;    // select the slice which is writing
    }

    int islice= (i+field->first) % field->field.size() ;   // include the rotation due to slippage

    if (s0==0){
      for (int j=0; j<ngrid*ngrid;j++){ 
	work[j]=scl*field->field.at(islice).at(j).real();
      }  
    }
    this->writeSingleNode(gid,"field-real",&work);     

    if (s0==0){
      for (int j=0; j<ngrid*ngrid;j++){ 
	work[j]=scl*field->field.at(islice).at(j).imag();
      }  
    }
    this->writeSingleNode(gid,"field-imag",&work);     


        
    H5Gclose(gid);
  }


  H5Fclose(fid);

 
  return;
}

void WriteFieldHDF5::writeGlobal(double reflen, double slicelen, double s0, double dx, int nx, int count)
{

  
  vector<double> tmp;
  tmp.resize(1);

  tmp[0]=reflen;
  this->writeSingleNode(fid,"wavelength",&tmp);
  tmp[0]=slicelen;
  this->writeSingleNode(fid,"slicespacing",&tmp);
  tmp[0]=s0;
  this->writeSingleNode(fid,"refposition",&tmp);
  tmp[0]=dx;
  this->writeSingleNode(fid,"gridsize",&tmp);
  
  vector<int> itmp;
  itmp.resize(1);

  itmp[0]=nx;
  this->writeSingleNodeInt(fid,"gridpoints",&itmp);

  itmp[0]=count;
  this->writeSingleNodeInt(fid,"slicecount",&itmp);
  

  
  return;
}


