
#include "writeBeamHDF5.h"

extern bool MPISingle;

#define SLICESELECT_DUMP    1
#define SLICESELECT_IGNORE  0

// constructor destructor
WriteBeamHDF5::WriteBeamHDF5()  {}
WriteBeamHDF5::~WriteBeamHDF5() {}


bool WriteBeamHDF5::write(string fileroot, Beam *beam, int stride)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node
  MPI_Comm_size(MPI_COMM_WORLD, &size); // assign rank to node
  if (MPISingle){
    size=1;
    rank=0;
  }

  string filename;
  filename = fileroot+".par.h5";
  if (rank == 0) { cout << "Writing particle distribution to file: " <<filename << " ..." << endl;} 
  if(!create_outfile(&fid, filename)) {
    return(false);
  }

  s0=rank;
  int ntotal=size*beam->beam.size();

  // write global data
  // this->writeGlobal(beam->nbins,beam->one4one,beam->reflength,beam->slicelength,beam->s0,ntotal);

  hid_t gid=H5Gcreate(fid,"Meta",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  this->writeVersion(gid);
  H5Gclose(gid);

  writeGlobal(beam, ntotal);

  /* write beam: loop through slices */
  vector<double> work,cur;
  cur.resize(1);
  int nwork=0;
  int npart=0;
  int local_nslice_written=0, global_nslice_written=0;
 
  int smin=rank*beam->beam.size();
  int smax=smin+beam->beam.size();

  for (int i=0; i<(ntotal);i++)
  {
    // if requested: selection of interesting slices
    if(beam->get_WriteFilter_active()) {
      int dodump = write_sliceselector(beam, i+1);
      if(dodump!=SLICESELECT_DUMP)
        continue;
    }

    local_nslice_written++;

    char name[16];
    sprintf(name,"slice%6.6d",i+1);
    gid=H5Gcreate(fid,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    int islice= i % beam->beam.size() ;   // count inside the given slice range

    // control which MPI rank writes its local particle data to the file
    s0=-1;
    if ((i>=smin) && (i<smax)){
      s0=0;    // select the slice which is writing
      npart = 0;
      for (int icount = 0; icount < beam->beam.at(islice).size(); icount += stride){
          npart++;
      }
    }

    int root = i /beam->beam.size();  // the current rank which sends the information of a slice
    if (size>1){
      MPI_Bcast(&npart,1,MPI_INT,root,MPI_COMM_WORLD);
    }

    if (npart != nwork){   // all cores do need to have the same length -> otherwise one4one crashes
	nwork=npart;
	work.resize(nwork);
    }

    if (s0==0){
      cur[0]=beam->current.at(islice);
    }
    this->writeSingleNode(gid,"current","A", &cur);



    //    if (nwork > 0 ){
      if (s0==0) {
	for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip*stride).gamma;}
      }
      this->writeSingleNode(gid,"gamma"  ," ", &work);

      if (s0==0) {
	for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip*stride).theta;}
      }
      this->writeSingleNode(gid,"theta"  ,"rad", &work);

      if (s0==0) {
	for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip*stride).x;}
      }
      this->writeSingleNode(gid,"x"  ,"m",&work);

      if (s0==0) {
	for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip*stride).y;}
      }
      this->writeSingleNode(gid,"y"  ,"m",&work);

      if (s0==0) {
	for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip*stride).px;}
      }
      this->writeSingleNode(gid,"px"  ,"rad",&work);

      if (s0==0) {
	for (int ip=0; ip<npart;ip++){work[ip]=beam->beam.at(islice).at(ip*stride).py;}
      }
      this->writeSingleNode(gid,"py"  ,"rad",&work);

      //  } 
    H5Gclose(gid);
  }

  H5Fclose(fid);

  // warn if no slice was written (wrong parameters for slice selector?)
  MPI_Reduce(&local_nslice_written, &global_nslice_written, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
  if(beam->get_WriteFilter_active()) {
    if ((global_nslice_written==0) && (rank==0)) {
      cout << "*** warning: no slice was written (beam dump slice selector parameters:"
           << " from=" << beam->get_WriteFilter_from()
           << " to=" << beam->get_WriteFilter_to()
           << " inc=" << beam->get_WriteFilter_inc()
           << ")" << endl;
    }
  }

  return(true);
}

// void WriteBeamHDF5::writeGlobal(int nbins,bool one4one, double reflen, double slicelen, double s0, int count)
void WriteBeamHDF5::writeGlobal(Beam *beam, int count)
{
  vector<double> tmp;
  tmp.resize(1);

  tmp[0]=beam->reflength;
  this->writeSingleNode(fid,"slicelength","m",&tmp);
  tmp[0]=beam->slicelength;
  this->writeSingleNode(fid,"slicespacing","m",&tmp);
  tmp[0]=beam->s0;
  this->writeSingleNode(fid,"refposition","m",&tmp);
  
  vector<int> itmp;
  itmp.resize(1);

  itmp[0]=beam->nbins;
  this->writeSingleNodeInt(fid,"beamletsize",&itmp);
  itmp[0]=count;
  this->writeSingleNodeInt(fid,"slicecount",&itmp);
  itmp[0]=static_cast<int>(beam->one4one);
  this->writeSingleNodeInt(fid,"one4one",&itmp);

  if(beam->get_WriteFilter_active())
  {
    int slice_min=-1;
    int slice_max=-1;
    int slice_inc = beam->get_WriteFilter_inc();

    write_sliceselector_geteffective_minmax(beam, &slice_min, &slice_max);

    itmp[0]=slice_min;
    writeSingleNodeInt(fid,"slicerange_min",&itmp);
    itmp[0]=slice_max;
    writeSingleNodeInt(fid,"slicerange_max",&itmp);
    itmp[0]=slice_inc;
    writeSingleNodeInt(fid,"slicerange_inc",&itmp);
  }
  
  return;
}

// filter function for selective dumping of slices
int WriteBeamHDF5::write_sliceselector(Beam *beam, int idslice)
{
  // note that slice count in this function begins with 1 (= as in HDF5 file containing beam dump)
  int slice_min=-1;
  int slice_max=-1;
  int inc = beam->get_WriteFilter_inc();

  write_sliceselector_geteffective_minmax(beam, &slice_min, &slice_max);

  if((slice_min<=idslice) && (idslice<=slice_max)) {
    int t = idslice-slice_min;
    if((t%inc)==0)
      return(SLICESELECT_DUMP);
  }
  return(SLICESELECT_IGNORE);
}

void WriteBeamHDF5::write_sliceselector_geteffective_minmax(Beam *beam, int *min, int *max)
{
  // note that slice count in this function begins with 1 (= as in HDF5 file containing beam dump)
  int slice_total = size*beam->beam.size(); // total slice number (all processes)
  int slice_min = 1;
  int slice_max = slice_total;

  // If the from/to range parameter is negative (=default),
  // then the global first/last slice determines the value.
  // Using the 'inc' parameter, one can the write a subset
  // of all slices, for instance every 10th slice, to the file. 
  if(beam->get_WriteFilter_from() >= 0)
    slice_min = beam->get_WriteFilter_from();
  if(beam->get_WriteFilter_to() >= 0)
    slice_max = beam->get_WriteFilter_to();

  *min = slice_min;
  *max = slice_max;
}
