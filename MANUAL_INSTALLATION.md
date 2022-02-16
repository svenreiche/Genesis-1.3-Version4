# GENESIS 1.3 - Compilation and Running Genesis 1.3


- [Compilation](#Compilation)
- [Running Genesis 1.3](#Running)


<a name="Compilation">**Compilation**</a>

Genesis supports now the automatic configuration with CMAKE. Following commands will build Genesis from source on Linux platforms. The minimal command to compile the code from the source code root directory is:
```
mkdir build
chdir build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
Depending on your compiler and system configuration, additional information may be given to CMAKE. Most important is the used compiler. This can be done with the additional definition ```-DCMAKE_CXX_COMPILER=####``` in the cmake command overwriting the default compiler, where ```###``` is the compiler, e.g. `mpicxx` or `CC`. If debugging is needed the CMAKE_BUILD_TARGET should be changed to `Debug`.
The executable is found in the build directory. For a successful build the libraries openmpi/mpich and parallel HDF5 are needed. The FFTW3 is recommended but might become mandatory in upcoming releases.


<a name="Running">**Running Genesis 1.3**</a>

Genesis 1.3 is a command line executable and requires at least a single input argument, which is the
filename of the input deck:

```
genesis4 [-o output-rootname] [-l lattice-filename] [-s seed] [-b beamline] input-filename
```

The arguments for the lattice filename, the beamline, the rootname for all output filenames and the seed for the random number generator are optional.
when defined they will overwrite any values in the ```setup``` namelist from the input file.
Starting Genesis in this
way will effectively invoke a parallel job with only one single core. To start the code with
just more than one core the launching command `mpirun` needs to be called, including the
number of cores. In its simplest call it would be:

```
mpirun -np xxx genesis4 [-o output-rootname] [-l lattice-filename] [-s seed] [-b beamline] input-filename
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



