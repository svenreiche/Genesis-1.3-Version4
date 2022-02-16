# GENESIS 1.3 - Current Development and Upcoming Features

Currently working on:

- Filtering of the source term in the field equation. Under certain circumstanses, some radiaiton modes are generated which are rather defined by the underlying grid, when these field components are reflected at the boundary. These modes, which diffracts quickly will have little interaction in the FEL process and can be damp with additional parameters
- The lattice will  have the option to use labels for numeric values, e.g. in a wiggler with ```aw=@aw```. The explicit value of ```@aw``` can be given within the input file, e.g. ```aw: VARIABLE=1.234;``` or in the main input deck. This will also the path to include undulator errors in the simulation
- A new command line argument will generate a semaphore file during execution which are deleted at the end of the run. This helps to check if runs have been sucessful or have been crashing.
 
