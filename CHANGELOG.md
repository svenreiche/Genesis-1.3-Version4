# Change Log

## Released:
### [4.5.1] - 2022/01/11 
Version 4.5.1 has been released. Part of it is this Change Log file, which will be updated in the future releases.

## Unreleased:
### [4.6.1-beta] - 2022/02/16
- Moved the calculation for the output file from the Beam and Field class to the class Diagnostic.cpp. This class interacts with the class  in Output.cpp so that the output class does not need to know the name of the field to be written. This informationis provided by the new class in Diagnostic.cpp. That way it gets easier to add additional output, where only on location needs to be changed (namely the derived classes in DiagnostUser.cpp)
- Command line arguments are now overwriting the correpsonding input the main input file. Supported are the output file root name, the lattice file, the beam line and the seed for the shot noise power.
- The compilation is changed to CMAKE, which allows for some configuration before compilation. It automatically searches for the required libraries and no manual configuration of the Makefile is needed any longer. 
- 20220218: New command line parser
- 20220218: Added support for semaphore file. If requested, this file is written at the end of a successful simulation run.
- 20220218: Exit code of GENESIS binary now depends on status (0 for success, 1 in case of error)
- 20220329: Added additional "start" semaphore file that is written after &setup block was processed


