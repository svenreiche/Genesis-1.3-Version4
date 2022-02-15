# GENESIS 1.3 - Compilation and Running Genesis 1.3


- [Compilation](#Compilation)
- [Running Genesis 1.3](#Running)


<a name="Compilation">**Compilation**</a>

In the root directory there is the file `Makefile` which is the makefile for the compilation on
Unix and Linux machines. This file can be used to help setting up the installation. In the file there is
the `OBJECTS` list, which contains all required source code files. There are located in several
directories in the subfolder `src`. To find them, the compiler search path has to be extended
as with the `VPATH` directive under unix. The compilation instruction should convert each
individual source file into an objects file and then link them with the HDF5 and OpenMPI
libraries. Note that at some systems the required definition of the include and libraries files
are wrapped in the compiler executive. E.g. in the given Makefile the compiler `h5pcc`
supports directly both required libraries. If this wrapper command does not exist on your system and both
libraries (HDF5 and OpenMPI) are not supported by the standard search path, they have
to be set explicitly with the `-I`, `-L`, and `-l` directive for the include search path, library search
path and the library name, respectively. Genesis 1.3 supports the library `FFTW`, enabling it
by setting the directive `-DFFTW` in the make file. It is needed if the user wants to calculate the instantaneous
divergence of the radiation field during run time.

If the make file is set-up correctly, type `make` under Unix systems to compile, which should
produce the executable `gencore` in the local directory. If source files have changed
the code can be compiled again. Only those files, which have been edited will be recompiled
and then linked with the unchanged object files. However, if some include files have changed
it is recommended to recompile from scratch. The command `make clean` will clean up the
directory by deleting all object files and the current executable. To install into your local
library, which is the `bin` folder in your home directory, just type `make install`. You can
edit the makefile to fit your system.

On platforms, which do not support GNU compiler directly, you have to build up your com-
pilation instruction according to the Makefile.

To complete the installation it is recommended to install the scripts in the directory `sdds2hdf`
in a directory (e.g. /bin), which is included in your search path. The files will convert El-
egant output distribution into HDF5 format. More info can be found later in the manual.
Also, the directory `xgenesis` includes some functions, which allows Matlab to parse the
output file of Genesis and to plot the results. They are not required for running Genesis
itself. To use them add the given folder to your Matlab Path or place the file in a directory,
which is already in the search path.

<a name="Running">**Running Genesis 1.3**</a>

Genesis 1.3 is a command line executable and requires at least a single input argument, which is the
filename of the input deck:

```
genesis4 [-o output-rootname] [-l lattice-filename] [-s seed] input-filename
```

The arguments for the lattice filename, the rootname for all output filenames and the seed for the random number generator are optional.
If defined, they act as the default values in the `setup` namelist and thus can be omitted in the main input file.
Any further information must be provided within the input deck. Starting Genesis in this
way will effectively invoke a parallel job with only one single core. To start the code with
just more than one core the launching command `mpirun` needs to be called, including the
number of cores. In its simplest call it would be:

```
mpirun -np xxx genesis4 [-o output-rootname] [-l lattice-filename] [-s seed] input-filename
```

where xxx is the number of cores to be requested.

Once Genesis started, the master core, which is the core with the MPI rank 0, will provide
the output to `stdout`. The output can vary for different input decks. Normally it is reported, once
the code generates or alters particle distributions or the radiation field, tracks the beam
through the lattice or when a time window is specified.

The execution produces output files whenever a tracking command is issued or a field or
particle dumped is requested. All in common is the root filename which is defined in the input
deck. Normal output files have the extension `.out.h5` for the first tracking or `.Runxxx.out.h`
for the xxxth run (except for the first, where the addition of ’Run’ is omitted). The file format
is HDF5. The explicit format is described in a later chapter of this manual. Field or particle
dumps are using the root filename as well as the default but under certain circumstances this
behavior can be overruled. For more information see the corresponding namelist elements
in the next chapter.

For time-dependent simulations the code requires a large amount of memories to keep the
entire radiation field and particle distribution in memory. This is the significant difference
to previous versions where the code was able to progress through the bunch sequentially.
However this excluded features such long range space charge fields, exchanging particles between slices, and
large delaying chicanes in the older versions. For an estimate of the memory demand the
radiation field requires

```
Nb=ns∗ng*ng∗16 Bytes,
```

where ns is the number of slices and ng the number of grid points. For the electron distribution
it is

```
Nb=ns∗ne∗48 Bytes,
```

with ne the number of macro particles per slice.
Both numbers should be added up and for safety multiplied by a factor of 2 for some
overhead (aligning arrays in memory, temporary working arrays, arrays to buffer the output
information etc.). The memory request will be distributed  over the nodes. Non time-dependent
simulation should be small enough to run on any computer.

In time-dependent runs the number of slices are equally distributed over all cores. The
number of slices are even increased to guarantee symmetry, extending the time-window of
the simulations. More information about that can be found in the next section. Often it
is useful to calculate the number of generated slices in advance, in particular if advanced
option such as subharmonic conversion are used. In the case of subharmonic conversion
the slices per core should be an integer of the subharmonic number. Otherwise artefacts of
numerically enhanced bunching causes wrong results.



