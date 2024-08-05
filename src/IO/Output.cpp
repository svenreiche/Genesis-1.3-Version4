#include "Output.h"
#include <stdlib.h>

// Meta info
#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#include <fstream>
#include <streambuf>



extern bool MPISingle;

Output::Output(){
  fid=-1;
}

Output::~Output(){}

void Output::close(){

  // output of unclosed elements

   int norphans = H5Fget_obj_count(fid, H5F_OBJ_ALL);
   if (norphans > 1) { // expect 1 for the file we have not closed
       int i;
       H5O_info_t info;
       char name[64];
       hid_t *objects = (hid_t *) calloc(norphans, sizeof(hid_t));
       H5Fget_obj_ids(fid, H5F_OBJ_ALL, -1, objects);
       for (i=0; i<norphans; i++) {
	 //           H5Oget_info(objects[i], &info);
           H5Iget_name(objects[i], name, 64);
	   cout << "Slice " <<s0/ds << " : " << i+1 << " of " << norphans << " things still open: " << objects[i] << " with name " << name << endl;
	   // took out the H5O call since it changes with the new hdf5 release
       }
   }

 H5Fclose(fid);


}


bool Output::open(string file, int s0_in, int ds_in)
{
  // s0 = slice number of first slice for a given node
  // ds = number of slices, which are kept by a given node

  // set record range and allocate working memory
  s0=s0_in;
  ds=ds_in;

#if 0
  // create the file for parallel access
  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);

  if (!MPISingle){
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  fid=H5Fcreate(file.c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,pid);
  H5Pclose(pid);
#else
  return(create_outfile(&fid, file));
#endif
}


void Output::writeMeta(Undulator *und)
{
  hid_t gid=H5Gcreate(fid,"Meta",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  writeMetaWorker(und, gid);
  H5Gclose(gid);
}

void Output::writeMetaWorker(Undulator *und, hid_t gid)
{
  this->writeVersion(gid);
  /*
  VersionInfo vi;
  hid_t gidsub;
  vector<double> tmp(1,0);

  gidsub=H5Gcreate(gid,"Version",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  tmp[0]=vi.Major();
  this->writeSingleNode(gidsub,"Major"," ",&tmp);
  tmp[0]=vi.Minor();
  this->writeSingleNode(gidsub,"Minor"," ",&tmp);
  tmp[0]=vi.Rev();
  this->writeSingleNode(gidsub,"Revision"," ",&tmp);
  tmp[0]=0;
  if (vi.isBeta()) { tmp[0]=1;}
  this->writeSingleNode(gidsub,"Beta"," ",&tmp);
  string s_bi(vi.Build());
  this->writeSingleNodeString(gidsub,"Build_Info", &s_bi);
  H5Gclose(gidsub);
  */

  time_t timer;
  time(&timer);
  string tim (ctime(&timer));
  this->writeSingleNodeString(gid,"TimeStamp", &tim);

  /* store selected environment variables */
  // const char *env_vars[] = {"HOST", "PWD", NULL};
  const char *env_vars[] = {"HOST", NULL};
  for(const char **p_env = &env_vars[0]; *p_env!=NULL; p_env++) {
     char *env_val = getenv(*p_env);
     string env = "Undefined";
     if (env_val != NULL){
       env=env_val;
     }
     this->writeSingleNodeString(gid, *p_env, &env);
  }


  struct passwd *pws;
  string user = "username lookup failed";
  if (NULL != (pws=getpwuid(getuid()))) // 'getpwuid' system call returns nullptr in case lookup was unsuccessful...
    user = pws->pw_name;
  this->writeSingleNodeString(gid,"User", &user);


  const int cwd_buflen = 4096;
  char cwd_buf[cwd_buflen];
  string cwd="getcwd call failed";
  if (getcwd(cwd_buf, cwd_buflen)!=NULL) {
    cwd = cwd_buf;
  }
  this->writeSingleNodeString(gid,"cwd", &cwd);


  /*** copy input files into .out.h5 file ***/
  ifstream inFile (meta_inputfile.c_str());
  stringstream buffer;
  buffer << inFile.rdbuf();//read the file
  string str=buffer.str();
  this->writeSingleNodeString(gid,"InputFile", &str);
  inFile.close();

  ifstream inFile2 (meta_latfile.c_str());
  stringstream buffer2;
  buffer2 << inFile2.rdbuf();//read the file
  string str2=buffer2.str();
  this->writeSingleNodeString(gid,"LatticeFile", &str2);
  inFile2.close();

  reportDumps(gid, und);
  reportPlugins(gid, und);
  reportMPI(gid);
}

void Output::writeGlobal(Undulator *und, double gamma, double lambda, double sample, double slen, bool one4one, bool time, bool scan, int ntotal)
{
  this->writeMeta(und);
  vector<double> tmp;
  tmp.resize(1);
  hid_t gid;

  gid=H5Gcreate(fid,"Global",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  tmp[0]=gamma;
  this->writeSingleNode(gid,"gamma0"," ",&tmp);
  tmp[0]=lambda;
  this->writeSingleNode(gid,"lambdaref","m",&tmp);
  tmp[0]=sample;
  this->writeSingleNode(gid,"sample"," ",&tmp);
  tmp[0]=slen;
  this->writeSingleNode(gid,"slen","m",&tmp);
  tmp[0] = one4one ? 1. : 0 ;
  this->writeSingleNode(gid,"one4one"," ",&tmp);
  tmp[0] = time ? 1. : 0 ;
  this->writeSingleNode(gid,"time"," ",&tmp);
  tmp[0] = scan ? 1. : 0 ;
  this->writeSingleNode(gid,"scan"," ",&tmp);

  tmp.resize(ntotal);
  for (int i=0; i<ntotal; i++){
    tmp[i]=static_cast<double>(i)*sample*lambda;
  }
  this->writeSingleNode(gid,"s","m",&tmp);

  double e0=1023.842e-9/lambda;
  double df=e0/sample/static_cast<double>(ntotal);
  if (ntotal ==1) {
      df=0;
  }
  e0 = e0-0.5*df*ntotal;
  for (int i=0; i<ntotal; i++){
    tmp[i]=e0+static_cast<double>(i)*df;
  }
  this->writeSingleNode(gid,"frequency","eV",&tmp);
  H5Gclose(gid);
}


void Output::reportDumps(hid_t gid, Undulator *und)
{
  hid_t gid_dr;
  vector<int> tmp(1);
  size_t ndumps;
  size_t k;

  gid_dr=H5Gcreate(gid,"Fielddumps",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  ndumps = und->fielddumps_filename.size();
  tmp[0] = ndumps;
  writeSingleNodeInt(gid_dr, "ndumps", &tmp);
  if(ndumps>0) {
    writeSingleNodeInt(gid_dr, "intstep", &und->fielddumps_intstep);  // writeSingleNodeInt does not work when passing empty data vector
  }
  for(k=0; k<ndumps; k++)
  {
    stringstream objname;

    objname << "filename" << (k+1);
    writeSingleNodeString(gid_dr, objname.str(), &und->fielddumps_filename.at(k));
/*    objname.str("");
    objname << "intstep" << (k+1);
    tmp[0] = und->fielddumps_intstep.at(k);
    writeSingleNodeInt(gid_dr, objname.str(), &tmp); */
  }
  H5Gclose(gid_dr);

  gid_dr=H5Gcreate(gid,"Beamdumps",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  ndumps = und->beamdumps_filename.size();
  tmp[0] = ndumps;
  writeSingleNodeInt(gid_dr, "ndumps", &tmp);
  if(ndumps>0) {
    writeSingleNodeInt(gid_dr, "intstep", &und->beamdumps_intstep);  // writeSingleNodeInt does not work when passing empty data vector
  }
  for(k=0; k<ndumps; k++)
  {
    stringstream objname;

    objname << "filename" << (k+1);
    writeSingleNodeString(gid_dr, objname.str(), &und->beamdumps_filename.at(k));
/*    objname.str("intstep");
    objname << (k+1);
    tmp[0] = und->beamdumps_intstep.at(k);
    writeSingleNodeInt(gid_dr, objname.str(), &tmp); */
  }
  H5Gclose(gid_dr);
}

void Output::reportPlugins(hid_t gid, Undulator *und)
{
  size_t nplugins=und->plugin_info_txt.size();
  if(nplugins==0)
    return;

  hid_t gid_sub=H5Gcreate(gid,"Plugins",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  for(int kk=0; kk<nplugins; kk++) {
    stringstream objname;
    objname << "info_txt." << (kk+1);
    writeSingleNodeString(gid_sub, objname.str(), &und->plugin_info_txt.at(kk));
    objname.str(""); objname.clear();
    objname << "prefix." << (kk+1);
    writeSingleNodeString(gid_sub, objname.str(), &und->plugin_hdf5_prefix.at(kk));
  }
  H5Gclose(gid_sub);
}

void Output::reportMPI(hid_t gid)
{
  vector<double> tmp;
  tmp.resize(1);

  int mpisize;
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  tmp[0] = mpisize;
  this->writeSingleNode(gid,"mpisize", " ", &tmp);
}

void Output::writeLattice(Beam * beam,Undulator *und)
{
  hid_t gid;
  gid=H5Gcreate(fid,"Lattice",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  this->writeSingleNode(gid,"zplot","m",&beam->zpos);
  this->writeSingleNode(gid,"z","m",&und->z);
  this->writeSingleNode(gid,"dz","m",&und->dz);
  this->writeSingleNode(gid,"aw"," ",&und->aw);
  this->writeSingleNode(gid,"ax","m",&und->ax);
  this->writeSingleNode(gid,"ay","m",&und->ay);
  this->writeSingleNode(gid,"ku","1/m",&und->ku);
  this->writeSingleNode(gid,"kx"," ",&und->kx);
  this->writeSingleNode(gid,"ky"," ",&und->ky);
  this->writeSingleNode(gid,"qf","1/m^2",&und->qf);
  this->writeSingleNode(gid,"qx","m",&und->qx);
  this->writeSingleNode(gid,"qy","m",&und->qy);
  this->writeSingleNode(gid,"cx","rad",&und->cx);
  this->writeSingleNode(gid,"cy","rad",&und->cy);
  this->writeSingleNode(gid,"gradx","1/m",&und->gradx);
  this->writeSingleNode(gid,"grady","1/m",&und->grady);
  this->writeSingleNode(gid,"slippage"," ",&und->slip);
  this->writeSingleNode(gid,"phaseshift","rad",&und->phaseshift);
  this->writeSingleNode(gid,"chic_angle","degree",&und->chic_angle);
  this->writeSingleNode(gid,"chic_lb","m",&und->chic_lb);
  this->writeSingleNode(gid,"chic_ld","m",&und->chic_ld);
  this->writeSingleNode(gid,"chic_lt","m",&und->chic_lt);

  H5Gclose(gid);
}


void Output::writeGroup(std::string group, std::map<std::string,std::vector<double> >& data, std::map<std::string,std::string>& units,  std::map<std::string,bool>& single) {
    hid_t gid=H5Gcreate(fid,group.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    for (auto &[key, val]: data) {
        this->writeDataset(gid, key, val, units[key], single[key]);
    }
    H5Gclose(gid);
    return;
}

void Output::writeDataset(hid_t gid, string path, std::vector<double> &val, std::string unit, bool single)
{
    std::size_t pos = path.find("/");
    if (pos == std::string::npos) {
        if (single) {
            this->writeSingleNode(gid, path.c_str(),unit.c_str(), &val);
        } else {
            this->writeBuffer(gid, path.c_str(), unit.c_str(), &val);
        }
    }  else {
        /* recursive procedure: split path of HDF5 obj to write and (1) open/create first part, then (2) continue (and split again, if needed) */
        string group = path.substr(0,pos);
        hid_t gsub;
        if (this->groupExists(gid,group)){
            gsub = H5Gopen1(gid,group.c_str());
        } else {
            gsub=H5Gcreate(gid,group.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        }
        this->writeDataset(gsub, path.substr(pos+1),val,unit,single);
        H5Gclose(gsub);
    }
}


void Output::writeBeamBuffer(Beam *beam)
{
  hid_t gid, gidsub;

  // step 1 - create the group
  gid=H5Gcreate(fid,"Beam",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  // step 2 - write individual datasets
  if (beam->outputEnergy()){
    this->writeBuffer(gid, "energy"," ",&beam->gavg);
    this->writeBuffer(gid, "energyspread"," ", &beam->gsig);
  }
  if (beam->outputSpatial()){
    this->writeBuffer(gid, "xposition","m",&beam->xavg);
    this->writeBuffer(gid, "yposition","m",&beam->yavg);
    this->writeBuffer(gid, "pxposition","rad", &beam->pxavg);
    this->writeBuffer(gid, "pyposition","rad", &beam->pyavg);
    this->writeBuffer(gid, "xsize","m", &beam->xsig);
    this->writeBuffer(gid, "ysize","m", &beam->ysig);
  }
  this->writeBuffer(gid, "bunching"," ",&beam->bunch);
  this->writeBuffer(gid, "bunchingphase","rad", &beam->bphi);
  if (beam->outputAux()){
    this->writeBuffer(gid, "efield","eV/m", &beam->efld);
  }
  
  this->writeBuffer(gid, "betax","m",&beam->bx);
  this->writeBuffer(gid, "betay","m",&beam->by);
  this->writeBuffer(gid, "alphax","rad",&beam->ax);
  this->writeBuffer(gid, "alphay","rad",&beam->ay);
  this->writeBuffer(gid, "emitx","m",&beam->ex);
  this->writeBuffer(gid, "emity","m",&beam->ey);
  this->writeBuffer(gid, "current","A",&beam->cu);

  int bh=beam->getBunchingHarmonics();
  char bgroup[32];
  for (int i=1; i<bh;i++){
    sprintf(bgroup,"bunching%d",(i+1));
    this->writeBuffer(gid, bgroup, " ",  &beam->bh[i-1]);
    sprintf(bgroup,"bunchingphase%d",(i+1));
    this->writeBuffer(gid, bgroup,"rad",  &beam->ph[i-1]);
    }

  if(beam->get_global_stat()) {
    gidsub=H5Gcreate(gid,"Global",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if (beam->outputEnergy()){
      this->writeSingleNode(gidsub,"energy"," ", &beam->tgavg);
      this->writeSingleNode(gidsub,"energyspread"," ", &beam->tgsig);
    }
    if (beam->outputSpatial()){
      this->writeSingleNode(gidsub,"xposition","m", &beam->txavg);
      this->writeSingleNode(gidsub,"xsize","m", &beam->txsig);
      this->writeSingleNode(gidsub,"yposition","m", &beam->tyavg);
      this->writeSingleNode(gidsub,"ysize","m", &beam->tysig);
    }
    H5Gclose(gidsub);  
  }
  

  // step 3 - close group and done
  H5Gclose(gid);

  return;
}


void Output::writeFieldBuffer(Field *field)
{



  // step 1 - create the group
  hid_t gid, gidsub;
  char name[10];

  int harm=field->harm;
  if (harm==1){
     sprintf(name,"Field");
  } else {
     sprintf(name,"Field%d",harm);
  }
  gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  // step 2 - write individual datasets
  this->writeBuffer(gid, "power","W",&field->power);
  if (field->outputSpatial()){
    this->writeBuffer(gid, "xsize","m",&field->xsig);
    this->writeBuffer(gid, "ysize","m",&field->ysig);
    this->writeBuffer(gid, "xposition","m",&field->xavg);
    this->writeBuffer(gid, "yposition","m",&field->yavg);
    }
#ifdef FFTW
  if (field->outputFFT()){
    this->writeBuffer(gid, "xdivergence","rad",&field->txsig);
    this->writeBuffer(gid, "ydivergence","rad",&field->tysig);
    this->writeBuffer(gid, "xpointing","rad",&field->txavg);
    this->writeBuffer(gid, "ypointing","rad",&field->tyavg);
  }
#endif
  if (field->outputIntensity()){
    this->writeBuffer(gid, "intensity-nearfield","arb unit",&field->nf_intensity);
    this->writeBuffer(gid, "phase-nearfield","rad", &field->nf_phi);
    this->writeBuffer(gid, "intensity-farfield","arb unit",&field->ff_intensity);
    this->writeBuffer(gid, "phase-farfield","rad",&field->ff_phi);
  }
  
  if(field->get_global_stat()) {
    gidsub=H5Gcreate(gid,"Global",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    this->writeSingleNode(gidsub,"energy"," ", &field->energy);
    if (field->outputSpatial()){ 
      this->writeSingleNode(gidsub,"xsize","m", &field->gl_xsig);
      this->writeSingleNode(gidsub,"ysize","m", &field->gl_ysig);
      this->writeSingleNode(gidsub,"xposition","m", &field->gl_xavg);
      this->writeSingleNode(gidsub,"yposition","m", &field->gl_yavg);
    }
    if (field->outputIntensity()){
      this->writeSingleNode(gidsub,"intensity-nearfield","arb unit", &field->gl_nf_intensity);
      this->writeSingleNode(gidsub,"intensity-farfield","arb unit ", &field->gl_ff_intensity);
    }
#ifdef FFTW
    if (field->outputFFT()){
      this->writeSingleNode(gidsub,"xdivergence","rad", &field->gl_txsig);
      this->writeSingleNode(gidsub,"ydivergence","rad", &field->gl_tysig);
      this->writeSingleNode(gidsub,"xpointing","rad", &field->gl_txavg);
      this->writeSingleNode(gidsub,"ypointing","rad", &field->gl_tyavg);
    }
#endif
    H5Gclose(gidsub);  
  }
  

  vector<double> tmp;
  tmp.resize(1);
  tmp[0]=field->dgrid;
  this->writeSingleNode(gid,"dgrid","m",&tmp);
  
  vector<int> itmp;
  itmp.resize(1);
  itmp[0]=field->ngrid;
  this->writeSingleNodeInt(gid,"ngrid",&itmp);



  // step 3 - close group and done






  H5Gclose(gid);
  return;

 

}


