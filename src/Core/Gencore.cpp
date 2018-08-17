#include "Gencore.h"

extern bool MPISingle;

#ifdef VTRACE
#include "vt_user.h"
#endif


int Gencore::run(const char *file, Beam *beam, vector<Field*> *field, Undulator *und,bool isTime, bool isScan, bool supressOutput)
{


        //-------------------------------------------------------
        // init MPI and get size etc.
        //
#ifdef VTRACE
  VT_TRACER("Core");
#endif  
        int size=1;
        int rank=0;

	if (!MPISingle){
           size=MPI::COMM_WORLD.Get_size(); // get size of cluster
           rank=MPI::COMM_WORLD.Get_rank(); // assign rank to node
        }

	if (rank==0) {
          cout << endl << "Running Core Simulation..." << endl;
        }

        //-----------------------------------------
	// init beam, field and undulator class

        Control   *control=new Control;

	control->init(rank,size,file,beam,field,und,isTime,isScan); 


        //------------------------------------------
        // main loop
	       	
	while(und->advance(rank)){
	  double delz=und->steplength();

	  // ----------------------------------------
	  // step 1 - apply most marker action  (always at beginning of a step)

	  bool sort=control->applyMarker(beam, field, und);


	  // ---------------------------------------
	  // step 2 - Advance electron beam

	  beam->track(delz,field,und);

	  // -----------------------------------------
	  // step 3 - Beam post processing, e.g. sorting


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

	  for (int i=0; i<field->size();i++){
	    field->at(i)->track(delz,beam,und);
          }


	  //-----------------------------------------
	  // step 5 - Apply slippage

	  for (int i=0;i<field->size();i++){
	    control->applySlippage(und->slippage(), field->at(i));  
	  }

	  //-------------------------------
	  // step 6 - Calculate beam parameter stored into a buffer for output

	  beam->diagnostics(und->outstep(),und->getz());
	  for (int i=0;i<field->size();i++){
	    field->at(i)->diagnostics(und->outstep());
	  }

          
        }
     
        //---------------------------
        // end and clean-up 

	// perform last marker action

	bool sort=control->applyMarker(beam, field, und);
	if (sort){
	    int shift=beam->sort();

	    if (shift!=0){
	      for (int i=0;i<field->size();i++){
		control->applySlippage(shift, field->at(i));  
	      }
	    }
	}


	// write out diagnostic arrays

	if (rank==0){
	  cout << "Writing output file..." << endl;
	}

        if (!supressOutput) {
   	   control->output(beam,field,und);
        }


	delete control;
      
        if (rank==0){
	  cout << endl << "Core Simulation done." << endl;
        }


        return 0;

}
