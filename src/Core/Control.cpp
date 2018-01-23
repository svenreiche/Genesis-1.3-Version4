#include <sstream>
#include "Control.h"
#include "writeFieldHDF5.h"
#include "writeBeamHDF5.h"


#ifdef VTRACE
#include "vt_user.h"
#endif



Control::Control()
{
  nwork=0;
}


Control::~Control()
{
}


bool Control::applyMarker(Beam *beam, vector<Field*>*field, Undulator *und)
{

  bool sort=false;

  int marker=und->getMarker();
  // possible file names
  int istepz=und->getStep();
  stringstream sroot;
  sroot << "." << istepz;

  if ((marker & 1) != 0){
    WriteFieldHDF5 dump;
    dump.write(root+sroot.str(),field);
  }
  
  if ((marker & 2) != 0){
    WriteBeamHDF5 dump;
    dump.write(root+sroot.str(),beam);
  }
  
  if ((marker & 4) != 0){
    sort=true;   // sorting is deferred after the particles have been pushed by Runge-Kutta
  }

  // bit value 8 is checked in und->advance()
  

  return sort;
}




bool Control::init(int inrank, int insize, const char *file, Beam *beam, vector<Field*> *field, Undulator *und, Output *out)
{

  rank=inrank;
  size=insize;
  accushift=0;
  stringstream sroot(file);
  root=sroot.str();
  root.resize(root.size()-7);  // remove the extension ".h5"

  hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(pid,MPI_COMM_WORLD,MPI_INFO_NULL);
  hid_t fid=H5Fopen(file,H5F_ACC_RDONLY,pid);
  H5Pclose(pid);


  int nzfld,nzpar,nsort;
  double zstop;

  readDataDouble(fid,(char *)"/Global/lambdaref",&reflen,1);
  readDataDouble(fid,(char *)"/Global/sample",&sample,1);
  readDataDouble(fid,(char *)"/Global/slen",&slen,1);
  readDataDouble(fid,(char *)"/Global/sref",&sref,1);
  readDataInt(fid, (char *)"/Global/nzout",&nzout,1);
  readDataInt(fid, (char *)"/Global/nzfld",&nzfld,1);
  readDataInt(fid, (char *)"/Global/nzpar",&nzpar,1);
  readDataInt(fid, (char *)"/Global/nsort",&nzpar,1);
  readDataDouble(fid,(char *)"/Global/zstop",&zstop,1);
 
  int tmp;
  readDataInt(fid, (char *)"/Global/one4one",&tmp,1);
  one4one=!(tmp==0);
  readDataInt(fid, (char *)"/Global/time",&tmp,1);
  timerun=!(tmp==0);
  readDataInt(fid, (char *)"/Global/scan",&tmp,1);
  scanrun=!(tmp==0);

 

  und->init(fid);  // read lattice from file
  und->updateMarker(nzfld,nzpar,nsort,zstop);   // include output from the &track namelist  
  
  H5Fclose(fid);


 


  MPI::COMM_WORLD.Barrier(); // synchronize all nodes till root has finish writing the output file

  nslice=beam->beam.size();
  MPI::COMM_WORLD.Bcast(&nslice,1,MPI::INT,0);
  noffset=rank*nslice;
  nslice=beam->beam.size();

  MPI::COMM_WORLD.Reduce(&nslice,&ntotal,1,MPI::INT,MPI::SUM,0);
  MPI::COMM_WORLD.Bcast(&ntotal,1,MPI::INT,0);


  //  cout <<"Rank: " << rank << " Slices: " << nslice << " Offset: " << noffset << endl;  

  if (rank==0){
    if(scanrun) { 
       cout << "Scan run with " << ntotal << " slices" << endl; 
    } else {
       if(timerun) { 
         cout << "Time-dependent run with " << ntotal << " slices" << " for a time window of " << slen*1e6 << " microns" << endl; 
       } else { 
         cout << "Steady-state run" << endl;
       }
    }
  }
  

  if (rank==0) { cout << "Opening Output File..."  << endl; }

  out->open(file,noffset,nslice);
  beam->initDiagnostics(und->outlength());
  beam->diagnostics(true,0);
  beam->diagnosticsStart();
  for (int i=0; i<field->size();i++){
      field->at(i)->initDiagnostics(und->outlength());
      field->at(i)->diagnostics(true);  // initial values
  }	

  return true;  
}



void Control::applySlippage(double slippage, Field *field)
{

#ifdef VTRACE
  VT_TRACER("Slippage");
#endif  

  if (timerun==false) { return; }

 
  // update accumulated slippage
  accushift+=slippage;

  // allocate working space

  if(nwork<field->ngrid*field->ngrid*2){
    nwork=field->ngrid*field->ngrid*2;
    work=new double [nwork];
  } 
  
  MPI::Status status;

  // following routine is applied if the required slippage is alrger than 80% of the sampling size

  int direction=1;

  while(abs(accushift)>(sample*0.8)){
      // check for anormal direction of slippage (backwards slippage)
      if (accushift<0) {direction=-1;} 

      accushift-=sample*direction;

      // get adjacent node before and after in chain
      int rank_next=rank+1;
      int rank_prev=rank-1;
      if (rank_next >= size ) { rank_next=0; }
      if (rank_prev < 0 ) { rank_prev = size-1; }	

      // for inverse direction swap targets
      if (direction<0) {
	int tmp=rank_next;
        rank_next=rank_prev;
        rank_prev=tmp; 
      }

      int tag=1;
   
      // get slice which is transmitted
      int last=(field->first+field->field.size()-1)  %  field->field.size();
      // get first slice for inverse direction
      if (direction<0){
	last=(last+1) % field->field.size();  //  this actually first because it is sent backwards
      }

      if (size>1){
        if ( (rank % 2)==0 ){                   // even nodes are sending first and then receiving field
           for (int i=0; i<nwork/2; i++){
	     work[2*i]  =field->field[last].at(i).real();
	     work[2*i+1]=field->field[last].at(i).imag();
	   }
	   MPI::COMM_WORLD.Send(work, nwork, MPI::DOUBLE, rank_next, tag);
	   MPI::COMM_WORLD.Recv(work, nwork, MPI::DOUBLE, rank_prev, tag, status);
	   for (int i=0; i<nwork/2; i++){
	     complex <double> ctemp=complex<double> (work[2*i],work[2*i+1]);
	     field->field[last].at(i)=ctemp;
	   }
	} else {                               // odd nodes are receiving first and then sending

	  MPI::COMM_WORLD.Recv(work, nwork, MPI::DOUBLE, rank_prev, tag, status);
	  for (int i=0; i<nwork/2; i++){
	    complex <double> ctemp=complex<double> (work[2*i],work[2*i+1]);
	    work[2*i]  =field->field[last].at(i).real();
	    work[2*i+1]=field->field[last].at(i).imag();
	    field->field[last].at(i)=ctemp;
	  }
	  MPI::COMM_WORLD.Send(work, nwork, MPI::DOUBLE, rank_next, tag);
	}
      }

      // first node has emptz field slipped into the time window
      if ((rank==0) && (direction >0)){
        for (int i=0; i<nwork/2;i++){
	  field->field[last].at(i)=complex<double> (0,0);
        }
      }

      if ((rank==(size-1)) && (direction <0)){
        for (int i=0; i<nwork/2;i++){
	  field->field[last].at(i)=complex<double> (0,0);
        }
      }

      // last was the last slice to be transmitted to the succeding node and then filled with the 
      // the field from the preceeding node, making it now the start of the field record.
      field->first=last;
      if (direction<0){
	field->first=(last+1) % field->field.size();
      }
  }

}
