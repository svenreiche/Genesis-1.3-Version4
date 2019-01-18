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
   if (norphans > 1) { /* expect 1 for the file we have not closed */
       int i;
       H5O_info_t info;
       char name[64];
       hid_t *objects = (hid_t *) calloc(norphans, sizeof(hid_t));
       H5Fget_obj_ids(fid, H5F_OBJ_ALL, -1, objects);
       for (i=0; i<norphans; i++) {
           H5Oget_info(objects[i], &info);
           H5Iget_name(objects[i], name, 64);
	   cout << "Slice " <<s0/ds << " : " << i+1 << " of " << norphans << " things still open: " << objects[i] << " with name " << name << " of type " << info.type << endl; 
	   //           printf("%d of %zd things still open: %d with name %s of type %d", i, norphans, objects[i], name, info.type);
       }
   }



  H5Fclose(fid);

  
}


void Output::open(string file, int s0_in, int ds_in)
{
  


  // ns = total number of slices
  // nz = total number of integration steps
  // s0 = slice number of first slice for a given node
  // ds = number of slices, which are kept by a given node

  // set record range and allocate working memory
  s0=s0_in;
  ds=ds_in;
  
  // create the file for parallel access
  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  if (!MPISingle){
    H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  }
  fid=H5Fcreate(file.c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,pid); 
  H5Pclose(pid);

}


void Output::writeMeta()
{



  hid_t gid,gidsub;
  vector<double> tmp;
  tmp.resize(1);

  gid=H5Gcreate(fid,"Meta",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  gidsub=H5Gcreate(gid,"Version",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  tmp[0]=versionmajor;
  this->writeSingleNode(gidsub,"Major",&tmp);
  tmp[0]=versionminor;
  this->writeSingleNode(gidsub,"Minor",&tmp);
  tmp[0]=versionrevision;
  this->writeSingleNode(gidsub,"Revision",&tmp);
  tmp[0]=0;
  if (versionbeta) { tmp[0]=1;}
  this->writeSingleNode(gidsub,"Beta",&tmp);
  H5Gclose(gidsub);  
  
  time_t timer;
  time(&timer);
  string tim (ctime(&timer));
  this->writeSingleNodeString(gid,"TimeStamp", &tim);


  struct passwd *pws;
  pws = getpwuid(getuid());
  string user (pws->pw_name);
  this->writeSingleNodeString(gid,"User", &user);


  ifstream inFile (meta_inputfile->c_str());
  stringstream buffer;
  buffer << inFile.rdbuf();//read the file
  string str=buffer.str();
  this->writeSingleNodeString(gid,"InputFile", &str);
  inFile.close();

  ifstream inFile2 (meta_latfile->c_str());
  stringstream buffer2;
  buffer2 << inFile2.rdbuf();//read the file
  string str2=buffer2.str();
  this->writeSingleNodeString(gid,"LatticeFile", &str2);
  inFile2.close(); 



  H5Gclose(gid);


}

void Output::writeGlobal(double gamma, double lambda, double sample, double slen, bool one4one, bool time, bool scan)
{



  this->writeMeta();   
  vector<double> tmp;
  tmp.resize(1);
  hid_t gid;

  gid=H5Gcreate(fid,"Global",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  tmp[0]=gamma;
  this->writeSingleNode(gid,"gamma0",&tmp);
  tmp[0]=lambda;
  this->writeSingleNode(gid,"lambdaref",&tmp);
  tmp[0]=sample;
  this->writeSingleNode(gid,"sample",&tmp);
  tmp[0]=slen;
  this->writeSingleNode(gid,"slen",&tmp);
  tmp[0] = one4one ? 1. : 0 ;
  this->writeSingleNode(gid,"one4one",&tmp);
  tmp[0] = time ? 1. : 0 ;
  this->writeSingleNode(gid,"time",&tmp);
  tmp[0] = scan ? 1. : 0 ;
  this->writeSingleNode(gid,"scan",&tmp);

  H5Gclose(gid);

}



 
void Output::writeLattice(Beam * beam,Undulator *und)
{



  hid_t gid;
  gid=H5Gcreate(fid,"Lattice",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  this->writeSingleNode(gid,"zplot",&beam->zpos);
  this->writeSingleNode(gid,"z",&und->z);
  this->writeSingleNode(gid,"dz",&und->dz);
  this->writeSingleNode(gid,"aw",&und->aw);
  this->writeSingleNode(gid,"ax",&und->ax);
  this->writeSingleNode(gid,"ay",&und->ay);
  this->writeSingleNode(gid,"ku",&und->ku);
  this->writeSingleNode(gid,"kx",&und->kx);
  this->writeSingleNode(gid,"ky",&und->ky);
  this->writeSingleNode(gid,"qf",&und->qf);
  this->writeSingleNode(gid,"qx",&und->qx);
  this->writeSingleNode(gid,"qy",&und->qy);
  this->writeSingleNode(gid,"cx",&und->cx);
  this->writeSingleNode(gid,"cy",&und->cy);
  this->writeSingleNode(gid,"gradx",&und->gradx);
  this->writeSingleNode(gid,"grady",&und->grady);
  this->writeSingleNode(gid,"slippage",&und->slip);
  this->writeSingleNode(gid,"phaseshift",&und->phaseshift);
  this->writeSingleNode(gid,"chic_angle",&und->chic_angle);
  this->writeSingleNode(gid,"chic_lb",&und->chic_lb);
  this->writeSingleNode(gid,"chic_ld",&und->chic_ld);
  this->writeSingleNode(gid,"chic_lt",&und->chic_lt);

  H5Gclose(gid);


}

void Output::writeBeamBuffer(Beam *beam)
{



  // step 1 - create the group
  hid_t gid;
  gid=H5Gcreate(fid,"Beam",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  // step 2 - write individual datasets
  this->writeBuffer(gid, "energy",&beam->gavg);
  this->writeBuffer(gid, "energyspread",&beam->gsig);
  this->writeBuffer(gid, "xposition",&beam->xavg);
  this->writeBuffer(gid, "yposition",&beam->yavg);
  this->writeBuffer(gid, "pxposition",&beam->pxavg);
  this->writeBuffer(gid, "pyposition",&beam->pyavg);
  this->writeBuffer(gid, "xsize",&beam->xsig);
  this->writeBuffer(gid, "ysize",&beam->ysig);
  this->writeBuffer(gid, "bunching",&beam->bunch);
  this->writeBuffer(gid, "bunchingphase",&beam->bphi);
  this->writeBuffer(gid, "efield",&beam->efld);
  
  this->writeBuffer(gid, "betax",&beam->bx);
  this->writeBuffer(gid, "betay",&beam->by);
  this->writeBuffer(gid, "alphax",&beam->ax);
  this->writeBuffer(gid, "alphay",&beam->ay);
  this->writeBuffer(gid, "emitx",&beam->ex);
  this->writeBuffer(gid, "emity",&beam->ey);
  this->writeBuffer(gid, "current",&beam->cu);

  int bh=beam->getBunchingHarmonics();
  char bgroup[20];
  for (int i=1; i<bh;i++){
    sprintf(bgroup,"bunching%d",(i+1));
    this->writeBuffer(gid, bgroup,  &beam->bh[i-1]);
    sprintf(bgroup,"bunchingphase%d",(i+1));
    this->writeBuffer(gid, bgroup,  &beam->ph[i-1]);
    }

  // step 3 - close group and done

  H5Gclose(gid);

  return;
}


void Output::writeFieldBuffer(Field *field)
{



  // step 1 - create the group
  hid_t gid;
  char name[10];

  int harm=field->harm;
  if (harm==1){
     sprintf(name,"Field");
  } else {
     sprintf(name,"Field%d",harm);
  }
  gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  // step 2 - write individual datasets
  this->writeBuffer(gid, "power",&field->power);
  this->writeBuffer(gid, "xsize",&field->xsig);
  this->writeBuffer(gid, "ysize",&field->ysig);
  this->writeBuffer(gid, "xposition",&field->xavg);
  this->writeBuffer(gid, "yposition",&field->yavg);
  this->writeBuffer(gid, "intensity-nearfield",&field->nf_intensity);
  this->writeBuffer(gid, "phase-nearfield",&field->nf_phi);
  this->writeBuffer(gid, "intensity-farfield",&field->ff_intensity);
  this->writeBuffer(gid, "phase-farfield",&field->ff_phi);



  vector<double> tmp;
  tmp.resize(1);
  tmp[0]=field->dgrid;
  this->writeSingleNode(gid,"dgrid",&tmp);
  
  vector<int> itmp;
  itmp.resize(1);
  itmp[0]=field->ngrid;
  this->writeSingleNodeInt(gid,"ngrid",&itmp);



  // step 3 - close group and done






  H5Gclose(gid);
  return;

 

}


