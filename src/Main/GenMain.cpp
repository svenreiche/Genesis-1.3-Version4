#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <ctime>


#include <fenv.h>
#include <signal.h>

#include <mpi.h>



// genesis headerfiles & classes
//#include "CodeTracing.h"

#include "Beam.h"
#include "Field.h"
#include "EField.h"

#include "Parser.h"
#include "GenProfile.h"
#include "Setup.h"
#include "AlterSetup.h"
#include "Lattice.h"
#include "GenTime.h"
#include "Gencore.h"
#include "LoadField.h"
#include "LoadBeam.h"
#include "AlterLattice.h"
#include "Track.h"
#include "SDDSBeam.h"
#include "SponRad.h"
#include "dump.h"
#include "ImportBeam.h"
#include "ImportField.h"
#include "writeBeamHDF5.h"
#include "writeFieldHDF5.h"
#include "Collective.h"
#include "Wake.h"
#include "Diagnostic.h"
#include "SemaFile.h"
#include "FieldManipulator.h"

#include <sstream>

using namespace std;

const double vacimp = 376.73;
const double eev    = 510999.06; 
const double ce     = 4.8032045e-11;

// info for the meta group in the hdf5 output file
string meta_inputfile;
string meta_latfile;

bool MPISingle;  // global variable to do mpic or not

//vector<string> event;
//vector<double> evtime;
//double evt0;

int genmain (string mainstring, map<string,string> &comarg, bool split) {

    meta_inputfile=mainstring;
    int ret=0;
    MPISingle=split;
	int rank,size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node
    MPI_Comm_size(MPI_COMM_WORLD, &size); // assign rank to node
	if (MPISingle){
	  rank=0;
	  size=1;
    }


    time_t timer;
	clock_t clocknow;
	clock_t clockstart = clock();
	//	evt0 = double(clockstart);
	//	event.push_back("start");
	//	evtime.push_back(0);

    if (rank==0){         // the output of version and build has been moved to the wrapper mainwrap.cpp
        time(&timer);
        cout << "Starting Time: " << ctime(&timer)<< endl;
        cout << "MPI-Comm Size: " << size << " nodes" << endl << endl;
    }

    //---------------------------------------------------------
    // Instance of beam and field class to carry the distribution

    vector<Field *> field;   // vector of various field components (harmonics, vertical/horizonthal components)
    Beam  *beam =new Beam;


    //----------------------------------------------------------
    // main loop extracting one element with arguments at a time
    //
    // It is assumed that (except for a few well-defined reasons)
    // the loop is only left if there is an error while processing
    // the input file.
    bool successful_run=false; // successful simulation run?
    Parser parser;
    string element;
    map<string,string> argument;

    // some dummy argument used earlier
    string latstring="";
    string outstring="";
    int in_seed=0;


    //-------------------------------------------
    // instances of main classes
    Setup *setup=new Setup;
	AlterLattice *alt=new AlterLattice;
    Lattice *lattice=new Lattice;
    Profile *profile=new Profile;
	Series  *seq    =new Series;
    Time *timewindow=new Time;
    FilterDiagnostics filter;

    //-----------------------------------------------------------
    // main loop for parsing
    parser.open(mainstring,rank);

        while(true){
          bool parser_result = parser.parse(&element,&argument);
          if(parser_result==false) {
            // parser returned 'false', analyse reason
            if(parser.fail()) {
              successful_run=false;
            } else {
              successful_run=true; // probably reached eof (nothing was parsed, but also no parsing error detected) => FIXME: implement more analysis
            }
            break;
          }

	  //----------------------------------------------
	  // log event
	  //	  clocknow=clock();
	  //	  event.push_back(element);
	  //	  evtime.push_back(double(clocknow-clockstart));

      //----------------------------------------------
	  // setup & parsing the lattice file

      if (element.compare("&setup")==0){
          // overwriting from the commandline input
          for (const auto &[key,val] :comarg){
              // consistency check: display info message if new element is added to argument map (it may fail in Setup::init)
              if(argument.find(key)==argument.end()) {
                  if(rank==0) {
                      cout << "info message: adding parameter '"<<key<<"' to &setup" << endl;
                  }
              }
              argument[key] = val;
          }
          if (!setup->init(rank,&argument,lattice, filter)){ break;}
          meta_latfile=setup->getLattice();

          /* successfully processed "&setup" block => generate start semaphore file if requested to do so by user */
          if (setup->getSemaEnStart()) {
            string semafile_fn_started;
            SemaFile sema;
            if (rank==0) {
              if(setup->getSemaFN(&semafile_fn_started)) {
                semafile_fn_started+="start"; // "started"-type semaphore file have additional filename suffix
                sema.put(semafile_fn_started);
                cout << endl << "generating 'start' semaphore file " << semafile_fn_started << endl << endl;
              } else {
                cout << endl << "error: not writing 'start' semaphore file, filename not defined" << endl << endl;
              }
            }
          }

          continue;
      }

      //----------------------------------------------
	  // modifying run

      if (element.compare("&alter_setup")==0){
	    AlterSetup *altersetup= new AlterSetup;
        if (!altersetup->init(rank,&argument,setup,lattice,timewindow,beam,&field)){ break;}
        delete altersetup;
        continue;
      }

      //----------------------------------------------
	  // modifying the lattice file

          if (element.compare("&lattice")==0){
            if (!alt->init(rank,size,&argument,lattice,setup,seq)){ break;}
            continue;  
          }  


          //---------------------------------------------------
          // adding sequence elements

          if ((element.compare("&sequence_const")==0)||
              (element.compare("&sequence_power")==0)||
              (element.compare("&sequence_random")==0)){
            if (!seq->init(rank,&argument,element)){ break; }
            continue;
	  }




          //---------------------------------------------------
          // adding profile elements

          if ((element.compare("&profile_const")==0)||
              (element.compare("&profile_gauss")==0)||
              (element.compare("&profile_file")==0)||
              (element.compare("&profile_file_multi")==0) ||
              (element.compare("&profile_polynom")==0)||
              (element.compare("&profile_step")==0)){
            if (!profile->init(rank,&argument,element)){ break; }
            continue;
	  }

          //----------------------------------------------------
          // defining the time window of simulation

	  if (element.compare("&time")==0){
            if (!timewindow->init(rank,size,&argument,setup)){ break;}
            continue;  
          }  

          //----------------------------------------------------
          // internal generation of the field

	  if (element.compare("&field")==0){
	    LoadField *loadfield=new LoadField;
            if (!loadfield->init(rank,size,&argument,&field,setup,timewindow,profile)){ break;}
	    delete loadfield;
            continue;  
          }  

          //----------------------------------------------------
          // field manipulation

	  if (element.compare("&field_manipulator")==0){
	    FieldManipulator *q = new FieldManipulator;
            if (!q->init(rank,size,&argument,&field,setup,timewindow,profile)){ break;}
	    delete q;
            continue;  
          }  
	 
          //----------------------------------------------------
          // setup of space charge field

	  if (element.compare("&efield")==0){
   	    EField *efield=new EField;
            if (!efield->init(rank,size,&argument,beam,setup,timewindow)){ break;}
	    delete efield;
            continue;  
          }  

          //----------------------------------------------------
          // setup of spontaneous radiation

	  if (element.compare("&sponrad")==0){
            SponRad *sponrad=new SponRad;
            if (!sponrad->init(rank,size,&argument,beam)){ break;}
	    delete sponrad;
            continue;  
          }  

          //----------------------------------------------------
          // setup wakefield

	  if (element.compare("&wake")==0){
	    Wake *wake = new Wake;
            if (!wake->init(rank,size,&argument,timewindow, setup, beam,profile)){ break;}
	    delete wake;
	    continue;  
          }  

          //----------------------------------------------------
          // internal generation of beam

	  if (element.compare("&beam")==0){
            LoadBeam *loadbeam=new LoadBeam;
            if (!loadbeam->init(rank,size,&argument,beam,setup,timewindow,profile,lattice)){ break;}
	    delete loadbeam;
            continue;  
          }  

          //----------------------------------------------------
          // external generation of beam with an sdds file

	  if (element.compare("&importdistribution")==0){
            SDDSBeam *sddsbeam=new SDDSBeam;
            if (!sddsbeam->init(rank,size,&argument,beam,setup,timewindow,lattice)){ break;}
	    delete sddsbeam;
            continue;  
          }  

          //----------------------------------------------------
          // tracking - the very core part of Genesis

	  if (element.compare("&track")==0){
            Track *track=new Track;
	    if (!track->init(rank,size,&argument,beam,&field,setup,lattice,alt,timewindow,filter)){ break;}
            delete track;
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&sort")==0){
	    beam->sort();
            continue;  
          }  


          //----------------------------------------------------
          // write beam, field or undulator to file

	  if (element.compare("&write")==0){
            Dump *dump=new Dump;
	    if (!dump->init(rank,size,&argument,setup,beam,&field)){ break;}
            delete dump;
            continue;  
          }  


          //----------------------------------------------------
          // import beam from a particle dump

	  if (element.compare("&importbeam")==0){
            ImportBeam *import=new ImportBeam;
	    if (!import->init(rank,size,&argument,beam,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  


          //----------------------------------------------------
          // import field from a field dump

	  if (element.compare("&importfield")==0){
            ImportField *import=new ImportField;
	    if (!import->init(rank,size,&argument,&field,setup,timewindow)){ break;}
            delete import;
            continue;  
          }  


          //----------------------------------------------------
          // stop execution of input file here (useful for debugging)

          if (element.compare("&stop")==0) {
            if (rank==0) {
              cout << endl << "*** &stop element: User requested end of simulation ***" << endl;
            }
            successful_run=true; // we reached "stop" command without error
            break;
          }



          //-----------------------------------------------------
          // error because the element typ is not defined

          if (rank==0){
            cout << "*** Error: Unknown element in input file: " << element << endl; 
	  }
          break;
        }

        /*
         * Semaphore file: Decide now what to do in the end (and prepare the needed information)
         * -> if requested by user (filename is provided by user *or* derived from rootname)
         * -> if run was successful
         *
         * Actual file generation is delayed until the very end as the
         * simulation could still crash when free-ing objects.
         */
        bool semafile_en=false;
        string semafile_fn;
        bool semafile_fn_valid=false;
        if ((semafile_en=setup->getSemaEnDone())) {
          if(successful_run) {
            if(setup->getSemaFN(&semafile_fn)) {
              semafile_fn_valid=true;
            }
          }
        }


        /*** clean up ***/
        delete timewindow;
        delete seq;
        delete profile;
        delete lattice;
        delete alt;
        delete setup;
        delete beam;

        // release memory allocated for fields
        for(int i=0; i<field.size(); i++)
          delete field[i];


        /* NOW, generate the semaphore file */
        if (semafile_en) {
          if(successful_run) {
            SemaFile sema;
            if (rank==0) {
              if(semafile_fn_valid) {
                sema.put(semafile_fn);
                cout << endl << "generating semaphore file " << semafile_fn << endl;
              } else {
                cout << endl << "error: not writing semaphore file, filename not defined" << endl;
              }
            }
          } else {
            if (rank==0) {
              cout << endl << "not writing semaphore file, as the run likely was not successful" << endl;
            }
          }
        }


	//	event.push_back("end");
	//	evtime.push_back(double(clocknow-clockstart));
	clocknow=clock();

 	if (rank==0) {

	  double elapsed_Sec=double(clocknow-clockstart)/CLOCKS_PER_SEC;

          time(&timer);
          cout << endl<< "Program is terminating..." << endl;
	  cout << "Ending Time: " << ctime(&timer);
	  cout << "Total Wall Clock Time: " << elapsed_Sec << " seconds" << endl;
          cout << "-------------------------------------" << endl;


	  /* tracing report
	  cout << "Tracing Summary" << endl;
	  cout << "==========================" << endl;
	  cout << setw(10) << "Event" << setw(8) << "dT (s)" << endl;
	  cout << "--------------------------" << endl;
	  for (int i=0; i<evtime.size(); i++){
	    cout << setw(10) << event[i] << setw(8) << evtime[i]/CLOCKS_PER_SEC << " " << endl;
	  }
	  cout << "--------------------------" << endl;
	  */
        }


        ret = (successful_run) ? 0 : 1;
        return(ret);
}
