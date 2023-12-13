# Change Log

### [4.6.5] - 2023/12/06
- added the ability to reduce the amount of particles when writing to a beam dump
- added the baility to apply bema shifts or matrix transformation to the electron beam by supplying the corresponding vectors and matrices by an HDF5 file 

### [4.6.4] - 2023/10/02
- bug fixed which didn't allow to compile Genesis without the FFTW3 library
- bug fixed for default behaviour of the undulator roll-off parameter kx and ky when only the type helical is defined in lattice file
- bug fix for accesing the vector for the z position when parsing the markers in the lattice file.

### [4.6.3] - 2023/09/12
- added output for the minimum and maximum value of the electron beam for each slice and step. This includes the parameters x,y, px, py and energy. The output is suprressed if the parameter exclude_aux_output is set to true
- Examples have been added for a step-by-step guide on setting up input decks.

### [4.6.2] - 2023/02/28
- improved class `Beam` to reduce memory consumption in simulations with harmonic upconversion

### [4.6.2] - 2022/12/30
- added parameter "outputdir" to "&setup": allows to define directory (relative or absolute) where to store the simulation output

### [4.6.2] - 2022/12/06
- add support for self-consistent calculation of the space charge field on the scale longer than the radiation wavelength. The calculation is controlled
with the boolean flag `longrange` in the efield namelist. It adds also the individual output dataset `LSCfield` in the `beam` group of the output file.
A second boolean flag `reducedLF` controls if the Lorentz Factor is reduced by the undulator field with the scaling factor of sqrt(1+aw^2).
At the moment I cannot decide if the correction is needed since the theory predicts it but experiments at SwissFEL with ESASE shows no dependence on aw when scanning the K-value between 0.8 and 3.7
So I leave in the option to use this correction or not.

### [4.6.1] - 2022/12/06
- add "field_manipulator" feature, currently it can be used to scale to power of the light field

### [4.6.1] - 2022/12/05
- Refactor the name "Time.h" and "Time.cpp" to "GenTime.h" and "GenTime.cpp" to prevent the conflict with the std-library Time.h
- Fixed in Diagnostic.cpp the evaluation of a complext value element with fabs. Replaced it with std::abs

### [4.6.1] - 2022/11/23
- Fixed the bug that importing harmonics with the importfield namelist cause always the error that wavelength and reference length weren't matching.  

### [4.6.1] - 2022/02/16
- Moved the calculation for the output file from the Beam and Field class to the class Diagnostic.cpp. This class interacts with the class  in Output.cpp so that the output class does not need to know the name of the field to be written. This informationis provided by the new class in Diagnostic.cpp. That way it gets easier to add additional output, where only on location needs to be changed (namely the derived classes in DiagnostUser.cpp)
- Command line arguments are now overwriting the correpsonding input the main input file. Supported are the output file root name, the lattice file, the beam line and the seed for the shot noise power.
- The compilation is changed to CMAKE, which allows for some configuration before compilation. It automatically searches for the required libraries and no manual configuration of the Makefile is needed any longer. 
- 20220218: New command line parser
- 20220218: Added support for semaphore file. If requested, this file is written at the end of a successful simulation run.
- 20220218: Exit code of GENESIS binary now depends on status (0 for success, 1 in case of error)
- 20220329: Added additional "start" semaphore file that is written after &setup block was processed


