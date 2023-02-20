# Change Log

## Released:
### [4.5.1] - 2022/01/11 
Version 4.5.1 has been released. Part of it is this Change Log file, which will be updated in the future releases.

## Unreleased:

### [4.6.2-beta] - 2022/12/30
- added parameter "outputdir" to "&setup": allows to define directory (relative or absolute) where to store the simulation output

### [4.6.2-beta] - 2022/12/06
- add support for self-consistent calculation of the space charge field on the scale longer than the radiation wavelength. The calculation is controlled
with the boolean flag `longrange` in the efield namelist. It adds also the individual output dataset `LSCfield` in the `beam` group of the output file.
A second boolean flag `reducedLF` controls if the Lorentz Factor is reduced by the undulator field with the scaling factor of sqrt(1+aw^2).
At the moment I cannot decide if the correction is needed since the theory predicts it but experiments at SwissFEL with ESASE shows no dependence on aw when scanning the K-value between 0.8 and 3.7
So I leave in the option to use this correction or not.

### [4.6.1-beta] - 2022/12/06
- add "field_manipulator" feature, currently it can be used to scale to power of the light field

### [4.6.1-beta] - 2022/12/05
- Refactor the name "Time.h" and "Time.cpp" to "GenTime.h" and "GenTime.cpp" to prevent the conflict with the std-library Time.h
- Fixed in Diagnostic.cpp the evaluation of a complext value element with fabs. Replaced it with std::abs

### [4.6.1-beta] - 2022/11/23
- Fixed the bug that importing harmonics with the importfield namelist cause always the error that wavelength and reference length weren't matching.  

### [4.6.1-beta] - 2022/02/16
- Moved the calculation for the output file from the Beam and Field class to the class Diagnostic.cpp. This class interacts with the class  in Output.cpp so that the output class does not need to know the name of the field to be written. This informationis provided by the new class in Diagnostic.cpp. That way it gets easier to add additional output, where only on location needs to be changed (namely the derived classes in DiagnostUser.cpp)
- Command line arguments are now overwriting the correpsonding input the main input file. Supported are the output file root name, the lattice file, the beam line and the seed for the shot noise power.
- The compilation is changed to CMAKE, which allows for some configuration before compilation. It automatically searches for the required libraries and no manual configuration of the Makefile is needed any longer. 
- 20220218: New command line parser
- 20220218: Added support for semaphore file. If requested, this file is written at the end of a successful simulation run.
- 20220218: Exit code of GENESIS binary now depends on status (0 for success, 1 in case of error)
- 20220329: Added additional "start" semaphore file that is written after &setup block was processed


