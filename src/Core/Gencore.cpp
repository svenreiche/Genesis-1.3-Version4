#include "Gencore.h"

extern bool MPISingle;



int Gencore::run(const char *file, Beam *beam, vector<Field*> *field, Undulator *und,bool isTime, bool isScan, FilterDiagnostics &filter)
{


    //-------------------------------------------------------
    // init MPI and get size etc.
    //
    int size=1;
    int rank=0;
	if (!MPISingle){
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node
	    MPI_Comm_size(MPI_COMM_WORLD, &size); // assign rank to node
    }

	if (rank==0) {
        cout << endl << "Running Core Simulation..." << endl;
    }

    //-----------------------------------------
	// init beam, field and undulator class

    Control   *control=new Control;
    control->init(rank,size,file,beam,field,und,isTime,isScan);

    Diagnostic diag;
    diag.init(rank, size, und->outlength(), beam->beam.size(),field->size(),isTime,isScan,filter);
    diag.calc(beam, field, und->getz());  // initial calculation

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

      if (und->outstep()) {
          diag.calc(beam, field, und->getz());
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

	diag.writeToOutputFile(file, beam,field,und);
    control->output(beam,field,und,diag);

	delete control;
      
    if (rank==0){
	  cout << endl << "Core Simulation done." << endl;
    }


    return 0;

}
