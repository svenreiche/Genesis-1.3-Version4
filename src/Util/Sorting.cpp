#include "Sorting.h"

extern bool MPISingle;
 
#include <algorithm>

Sorting::Sorting()
{
  dosort=false;
  doshift=false;
}

Sorting::~Sorting(){
}



void Sorting::init(int rank_in,int size_in, bool doshift_in, bool dosort_in)
{
  rank=rank_in;  // rank of the node (0-size-1)
  size=size_in;  // size of the mpi run
  dosort=dosort_in;        // flag to do the sorting at all
  doshift=doshift_in;      // flag to apply the global shifting, meaning that the radiationfield is rather moved than the beam
  return;

}


void Sorting::configure(double s0_in, double slicelen_in, double sendmin_in, double sendmax_in, double keepmin_in, double keepmax_in, bool frame_in)
{
  s0=s0_in;      // offset in particle distribution of the given node
  slen=slicelen_in; // distance between two adjacent slice in particle distribution vector
  sendmin=sendmin_in;  // s < sendmin -> particle is send to previous node
  sendmax=sendmax_in;  // same in the other direction
  keepmin=keepmin_in;  // s < keepmin -> particle is deleted from internal slice
  keepmax=keepmax_in;  // same in the other direction
  globalframe=frame_in; // flag whether particle position is global or not
  return;

}





int Sorting::sort(vector <vector <Particle> > * recdat){

  //   cout << "Rank: " << rank << " s0: " << s0 << " slen: " << slen << " sendmin: " << sendmin << " sendmax: " << sendmax << " keepmin: " << keepmin << " keepmax: " << keepmax << " globaL : " << globalframe << endl;

  if (!dosort) {
    if (rank==0) {cout << "*** Warning: Sorting only enabled for one-2-one simulations" << endl;}
    return 0;
  }
  if (rank==0) {cout << "Sorting..." << endl; }
  int shift =0;
  // step one - calculate global shift and see whether it can be adjusted by shifting radiation instead
  //////  int shift = this->centerShift(recdat);


  // step two - push outside particles to other nodes
  this->globalSort(recdat);
  // step three - sort within the given node
  this->localSort(recdat);


  return shift;

}


void Sorting::localSort(vector <vector <Particle> > * recdat)  // most arguments are now part of the class
{
 
  Particle p;  

  double invslen=1./slen;
  cout << "Testfloor: " << floor(-0.2) << " and " << static_cast<int>(floor(-0.2)) << endl;
  // note that global sorting comes first. Therefore all particles are staying in the same domain 
  
  vector<int> count,count2;
  count.resize(recdat->size());
  count2.resize(recdat->size());
  for  (int i=0; i< count.size(); i++){ 
     count[i]=0;
     count2[i]=0;
  }


  for (int a=0;a<recdat->size();a++){  //Run over the slices 
    int b=0;
    while ( b < recdat->at(a).size()){
      double theta=recdat->at(a).at(b).theta;
      int atar=static_cast<int>(floor(theta*invslen));   // relative target slice. atar = 0 -> stays in same slice

      if (atar!=0) {     // particle not in the same slice

	  p.theta=recdat->at(a).at(b).theta-slen*(atar);
	  p.gamma=recdat->at(a).at(b).gamma;
	  p.x    =recdat->at(a).at(b).x;
	  p.y    =recdat->at(a).at(b).y;
	  p.px   =recdat->at(a).at(b).px;
	  p.py   =recdat->at(a).at(b).py;
	  recdat->at(a+atar).push_back(p);  // pushing particle in correct slice

          count[a+atar]+=1;
          count2[a]-=1;

	  // eliminating the current particle
 	  int ilast=recdat->at(a).size()-1;
      	  recdat->at(a).at(b).theta=recdat->at(a).at(ilast).theta;
	  recdat->at(a).at(b).gamma=recdat->at(a).at(ilast).gamma;
	  recdat->at(a).at(b).x    =recdat->at(a).at(ilast).x;
	  recdat->at(a).at(b).y    =recdat->at(a).at(ilast).y;
	  recdat->at(a).at(b).px   =recdat->at(a).at(ilast).px;
	  recdat->at(a).at(b).py   =recdat->at(a).at(ilast).py;
	  recdat->at(a).pop_back();
      } else {
	b++;
      }
    }
  }	  
	  //	  count++;
      
    



  for (int i=0; i<count.size(); i++){
    cout<<"Rank: " << rank << " Slice: " << i << " received: " << count[i] << " send: "<<count2[i]<<endl;
  } 
   
  return;
}



// routine which moves all particles, which are misplaced in the given domain of the node to other nodes.
// the methods is an iterative bubble sort, pushing excess particles to next node. There the fitting particles are collected the rest moved further.

void Sorting::globalSort(vector <vector <Particle> > *rec)
{

  this->fillPushVectors(rec);   // here is the actual sorting to fill the vectore pushforward and pushbackward
  if (rank==(size-1)) { pushforward.clear(); }        
  if (rank==0) { pushbackward.clear(); }
  if (size==1) { return; } // no need to transfer if only one node is used.
  
  cout << "Rank: " << rank << " - Forward: " << pushforward.size()/6 << " - Backward: " << pushbackward.size()/6 << endl;

  int maxiter=size-1;  
  int nforward=pushforward.size();
  int nbackward=pushbackward.size();
  int ntotal=nforward+nbackward;
  int nreduce=0;
  MPI::COMM_WORLD.Allreduce(&ntotal,&nreduce,1,MPI::INT,MPI::SUM);  // get the info on total number of particles shifted among nodes

  if (nreduce == 0){ return; }	

  while(maxiter>0){
    if (rank==0) {cout << "Sorting: Transferring " << nreduce/6 << " particles to other nodes at iteration " << size-maxiter << endl;}



  // step one - pairing ranks (0,1) (2,3) etc, the last rank, if uneven last element should clear its pushforward.
     bool transfer = true;
     if (((rank % 2) == 0) && (rank == (size -1))) { transfer = false; }  // last rank for an uneven count of cores -> no forward transmission

     if ((rank % 2)==0) {   // even ranks sending
       if (transfer) {this->send(rank+1,&pushforward);}        // sends its forward particles to higher node
       pushforward.clear();                                    // no need for the record, data has been sent
       if (transfer) {this->recv(rank+1,rec,&pushbackward);}   // catch particles which are sent back. In self->recv either get the particles or put them in backword array
     }  else {   // odd ranks receiving - there will be always a smaller even rank therefore always receiving
	   this->recv(rank-1,rec,&pushforward);
           this->send(rank-1,&pushbackward);
	   pushbackward.clear();
     }
      
  // step two - pairing ranks (1,2) (3,4) etc

     transfer == true;
     if (((rank % 2) == 1) && (rank == (size -1))) { transfer = false; }  // last core for an even number
     if (rank==0) { transfer = false; }                                  // as well as first core.
     
     if ((rank % 2)==1) {  // odd ranks sending
       if (transfer) {this->send(rank+1,&pushforward);}
	 pushforward.clear();
	 if (transfer) {this->recv(rank+1,rec,&pushbackward);}
     } else {  // even one receiving - here check for the very first node
       if (transfer){
	   this->recv(rank-1,rec,&pushforward); 
           this->send(rank-1,&pushbackward);
           pushbackward.clear();
       }
     }
 
     maxiter--;
     nforward=pushforward.size();
     nbackward=pushbackward.size();
     ntotal=nforward+nbackward;
     nreduce=0;
     MPI::COMM_WORLD.Allreduce(&ntotal,&nreduce,1,MPI::INT,MPI::SUM);
     if (nreduce == 0){  return; }
  }
  pushforward.clear();
  pushbackward.clear();
  return;



  

}


void Sorting::send(int target, vector<double> *data)
{
  int ndata=data->size();
  MPI::COMM_WORLD.Send( &ndata,1,MPI::INT,target,0);
  if (ndata == 0) {
    return;
  }
  MPI::COMM_WORLD.Send(&data->front(),ndata,MPI::DOUBLE,target,0);
}


 void Sorting::recv(int source, vector <vector <Particle> > *rec ,vector<double> *olddata)
{

  double shift=slen;
  if(globalframe){shift=0;}

   MPI::Status status;
   int ndata=0;


   MPI::COMM_WORLD.Recv(&ndata,1,MPI::INT,source,0,status);
   if (ndata==0) {  // no data received.
     return;
   }

   vector<double> data (ndata);
   MPI::COMM_WORLD.Recv(&data.front(),ndata,MPI::DOUBLE,source,0,status); // geting the particles from adjacent node.
 

   //Determines whether the data needs to be pushed forward or backwards or stored in the correct slices
   int np=rec->size();                      // number of slices
   if (source>rank) {                       // data are coming from higher node -> particles are pushed backward
     for (int a=0;a<ndata;a+=6) {
       double s = s0+slen*(np-1)+data[a];  // get the actual positionassume that backward the particles are placed in the last slice !
       if (s<sendmin){
         olddata->push_back(data[a]+shift*np);  // if the particle isstillpushed through than it phase is adjusted by slen*slicenumber
	 for (int b=1;b<6;b++){
           olddata->push_back(data[a+b]);
         }
       } 
       if (s>=keepmin){
         Particle par;
         par.theta=data[a];
         par.gamma=data[a+1];
         par.x    =data[a+2];
         par.y    =data[a+3];
         par.px   =data[a+4];
         par.py   =data[a+5];
         rec->at(np-1).push_back(par); // add at last slice;
       }
     }
   } else {    // particles are coming from a lower node -> particles are pushed forward
     for (int a=0;a<ndata;a+=6) {
       double s = s0+data[a];  // get the actual position assume that forward the particles are placed in the first slice !
       if (s>sendmax){
	 olddata->push_back(data[a]-shift*np);
	 for (int b=1;b<6;b++){
           olddata->push_back(data[a+b]);
         }
       } 
       if (s<=keepmax){
         Particle par;
         par.theta=data[a];
         par.gamma=data[a+1];
         par.x    =data[a+2];
         par.y    =data[a+3];
         par.px   =data[a+4];
         par.py   =data[a+5];
         rec->at(0).push_back(par); // add to first slice
       }
     }
   }

   return;
}




void Sorting::fillPushVectors(vector< vector <Particle> >*rec)
{
  cout << "Rank: " << rank << " sendmin: " << sendmin << " sendmax: " << sendmax << endl;
  cout << "Rank: " << rank << " keepmin: " << keepmin << " keepmax: " << keepmax << endl;
  //step one - fill the push vectors
  pushforward.clear();
  pushbackward.clear();
  
  int nsize=rec->size();
  double shift=slen;  // flag to indicate correction in position because each slice has its own position ( 3pi in slice 5 is pi in slice 6}
  if (globalframe) {shift = 0;} // don't change position if it is a global frame (e.g. when importing elegant distibution)
  
  int count = 0;

  for (int i = 0; i < nsize; i++){  // loop over slices
    for (int j = 0; j < rec->at(i).size(); j++){ // loop over particles in slice
      double s = s0+slen*i+rec->at(i).at(j).theta;  // get the actual position
      if (s<sendmin){
	pushbackward.push_back(rec->at(i).at(j).theta+(i+1)*shift); 
	// example: in first slice phase is -1. Particle needs to be sent to previous node to the last slice there. Therefore
        // the offset is just the slice  length (normally 2 pi). In the second slice it needs to be a value of -7 < - 2 pi.
	pushbackward.push_back(rec->at(i).at(j).gamma);
	pushbackward.push_back(rec->at(i).at(j).x);
	pushbackward.push_back(rec->at(i).at(j).y);
	pushbackward.push_back(rec->at(i).at(j).px);
	pushbackward.push_back(rec->at(i).at(j).py);
      }
      if (s>sendmax){
	pushforward.push_back(rec->at(i).at(j).theta-(nsize-i)*shift); 
	// example: last slice has value of 7 > 2 pi to be pushed to next node in the first slice there. The offset correction would be - 2 pi, which is given by
	// -(nlen-i)*slicelength  with i = nlen-1 for the last slice
	pushforward.push_back(rec->at(i).at(j).gamma);
	pushforward.push_back(rec->at(i).at(j).x);
	pushforward.push_back(rec->at(i).at(j).y);
	pushforward.push_back(rec->at(i).at(j).px);
	pushforward.push_back(rec->at(i).at(j).py);

      }
    }

    int j=0;

    while (j<rec->at(i).size()){
        double s = s0+slen*i+rec->at(i).at(j).theta;  // get the actual position
        if ((s<keepmin)||(s>keepmax)){
          count++;
 	  int ilast=rec->at(i).size()-1;
	  rec->at(i).at(j).theta=rec->at(i).at(ilast).theta;
	  rec->at(i).at(j).gamma=rec->at(i).at(ilast).gamma;
	  rec->at(i).at(j).x    =rec->at(i).at(ilast).x;
	  rec->at(i).at(j).y    =rec->at(i).at(ilast).y;
	  rec->at(i).at(j).px   =rec->at(i).at(ilast).px;
	  rec->at(i).at(j).py   =rec->at(i).at(ilast).py;
	  rec->at(i).pop_back();
	} else {
	  j++;
	}
    }

  }
  cout << "Rank: " <<rank << " Deleted: " << count << " Forward: " << pushforward.size()/6 << " Backward: " <<pushbackward.size()/6 << endl;
}


int Sorting::centerShift(vector <vector <Particle> > * recdat)
{
  if (!doshift){ return 0;}

  double invslen=1./slen;
  double shift = 0;
  double part  = 0;

  // calculate on the node the amount of total transferred particles

  for (int a=0;a<recdat->size();a++){  //Run over the slices 
    for (int b=0;b<recdat->at(a).size();b++) {  //Loop over the partiles in the slice
      part+=1.0;
      shift+=floor(recdat->at(a).at(b).theta*invslen);   // relative target slice. = 0 -> stays in same slice
    }
  }

  double tshift,tpart;
  if (MPISingle){
    tshift=shift;
    tpart=part;

  } else {
     MPI::COMM_WORLD.Allreduce(&shift,&tshift,1,MPI::DOUBLE,MPI::SUM);  // get the info on total number of particles shifted among nodes
     MPI::COMM_WORLD.Allreduce(&part,&tpart,1,MPI::DOUBLE,MPI::SUM);  // get the info on total number of particles shifted among nodes
  }

  int nshift=-static_cast<int>(round(tshift/tpart));  // instead of moving more than half the particles forward, it is easier to move the field backwards. Also theta needs to be corrected

  for (int a=0;a<recdat->size();a++){  //Run over the slices 
    for (int b=0;b<recdat->at(a).size();b++) {  //Loop over the partiles in the slice
      recdat->at(a).at(b).theta+=static_cast<double>(nshift)*slen;   // adjust theta position because the field is moved instead of particles
    }
  }

  return nshift;
}
