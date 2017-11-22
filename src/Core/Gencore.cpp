#include "Gencore.h"
#include "MPEProfiling.h"
int Gencore::run(const char *file, Beam *beam, vector<Field*> *field)
{


        //-------------------------------------------------------
        // init MPI and get size etc.
        //


        int size=1;
        int rank=0;
       
        size=MPI::COMM_WORLD.Get_size(); // get size of cluster
        rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node


	if (rank==0) {
          cout << endl << "Running Core Simulation..." << endl;
        }

        //-----------------------------------------
	// init beam, field and undulator class

        Control   *control=new Control;;
        Undulator *und=new Undulator;
        Output    *out=new Output;

	control->init(rank,size,file,beam,field,und,out); 


        //------------------------------------------
        // main loop
	       	
	while(und->advance(rank)){
	  double delz=und->steplength();

	  // ----------------------------------------
	  // step 1 - apply most marker action  (always at beginning of a step)

	  //	  bool sort=control->applyMarker(beam, field, und);


	  // ---------------------------------------
	  // step 2 - Advance electron beam

	  mpe.logCalc(false,true,"Core Calculation");
	  beam->track(delz,field,und);
	  mpe.logCalc(true,true,"Core Calculation");

	  // -----------------------------------------
	  // step 3 - Beam post processing, e.g. sorting

	  bool sort=false;

	  if (sort){
	    int shift=beam->sort();

	    if (shift!=0){
	      for (int i=0;i<field->size();i++){
		control->applySlippage(shift, field->at(i));  
	      }
	    }
	  }
  
	  // ---------------------------------------
	  // step 4 - Advance radiation field

	  mpe.logCalc(false,false,"Core Calculation");
	  for (int i=0; i<field->size();i++){
	    field->at(i)->track(delz,beam,und);
          }
          mpe.logCalc(true,false,"Core Calculation");


	  //-----------------------------------------
	  // step 5 - Apply slippage

	  for (int i=0;i<field->size();i++){
	    control->applySlippage(und->slippage(), field->at(i));  
	  }

	  //-------------------------------
	  // step 6 - Calculate beam parameter stored into a buffer for output

	  mpe.logCalc(false,false,"Diagnostic Calculation");
	  beam->diagnostics(und->outstep(),und->getz());
	  for (int i=0;i<field->size();i++){
	    field->at(i)->diagnostics(und->outstep());
	  }
          mpe.logCalc(true,false,"Diagnostic Calculation");


          
        }
       
     
        //---------------------------
        // end and clean-up 

        mpe.logIO(false,true,"Write Diagnostic Output to File");  
	out->writeBeamBuffer(beam);
	for (int i=0; i<field->size();i++){
	  out->writeFieldBuffer(field->at(i));
        }
        mpe.logIO(true,true,"Write Diagnostic Output to File");


      	out->close();
	//	if (rank==0){
	  

        delete und;
        delete out;  
      
        if (rank==0){
	  cout << endl << "Core Simulation done." << endl;
        }


        return 0;

}
