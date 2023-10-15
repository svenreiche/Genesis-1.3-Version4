#include "Gencore.h"

#ifdef USE_DPI
  #include "DiagnosticHookS.h"
#endif

extern bool MPISingle;

bool Gencore::run(const char *file, Beam *beam, vector<Field*> *field, Setup *setup, Undulator *und,bool isTime, bool isScan, FilterDiagnostics &filter)
{
    // function returns 'true' if everything is ok


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
#ifdef USE_DPI
    for(int kk=0; kk<setup->diagpluginfield_.size(); kk++) {
	if(rank==0) {
            cout << "Setting up DiagFieldHook for libfile=\"" << setup->diagpluginfield_.at(kk).libfile << "\", obj_prefix=\"" << setup->diagpluginfield_.at(kk).obj_prefix << "\"" << endl;
        }
        DiagFieldHook *pdfh = new DiagFieldHook(); /* !do not delete this instance, it will be destroyed when DiagFieldHook instance is deleted! */
        bool diaghook_ok = pdfh->init(&setup->diagpluginfield_.at(kk));
	pdfh->set_runid(setup->getCount()); // propagate run id so that it can be used in the plugins, for instance for filename generation
        if(diaghook_ok) {
            diag.add_field_diag(pdfh);
            if(rank==0) {
                cout << "DONE: Registered DiagFieldHook" << endl;
            }
        } else {
            delete pdfh;
            if(rank==0) {
                cout << "failed to set up DiagFieldHook, not registering" << endl;
            }
        }
    }

	for(int kk=0; kk<setup->diagpluginbeam_.size(); kk++) {
		if(rank==0) {
			cout << "Setting up DiagBeamHook for libfile=\"" << setup->diagpluginbeam_.at(kk).libfile << "\", obj_prefix=\"" << setup->diagpluginbeam_.at(kk).obj_prefix << "\"" << endl;
		}
		DiagBeamHook *pdbh = new DiagBeamHook(); /* !do not delete this instance, it will be destroyed when DiagBeamHook instance is deleted! */
		bool diaghook_ok = pdbh->init(&setup->diagpluginbeam_.at(kk));
		pdbh->set_runid(setup->getCount()); // propagate run id so that it can be used in the plugins, for instance for filename generation
		if(diaghook_ok) {
			diag.add_beam_diag(pdbh);
			if(rank==0) {
				cout << "DONE: Registered DiagBeamHook" << endl;
			}
		} else {
			delete pdbh;
			if(rank==0) {
				cout << "failed to set up DiagBeamHook, not registering" << endl;
			}
		}
	}
#endif
    diag.init(rank, size, und->outlength(), beam->beam.size(),field->size(),isTime,isScan,filter);
    diag.calc(beam, field, und->getz());  // initial calculation

	/*************/
	/* MAIN LOOP */
	/*************/
	while(und->advance(rank))
	{
	  double delz=und->steplength();

	  // ----------------------------------------
	  // step 1 - apply most marker action  (always at beginning of a step)
	  bool error_IO=false;
	  bool sort=control->applyMarker(beam, field, und, error_IO);
	  if(error_IO) {
	    return(false);
	  }


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

	  //beam->diagnostics(und->outstep(),und->getz());
	  //for (int i=0;i<field->size();i++){
	  //  field->at(i)->diagnostics(und->outstep());
	  //}

	  if (und->outstep()) {
	    diag.calc(beam, field, und->getz());
	  }
	}
     
        //---------------------------
        // end and clean-up 

	// perform last marker action
        bool error_IO=false;
	bool sort=control->applyMarker(beam, field, und, error_IO);
	if(error_IO) {
	  return(false);
	}
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
	// control->output(beam,field,und,diag);
	if(!diag.writeToOutputFile(file, beam,field,und)) {
	  delete control;
	  return(false);
	}

	delete control;
      
    if (rank==0){
	  cout << endl << "Core Simulation done." << endl;
    }
    return(true);
}
