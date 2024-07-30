# GENESIS 1.3 - Compilation and Running Genesis 1.3


- [Conda](#Conda)
- [Compilation](#compilation)
- [Running](#running)


## Conda

### Using Miniforge

For users new to conda, we recommend downloading and installing the latest [Miniforge](https://conda-forge.org/miniforge/) for your platform.

This is a customized version of miniconda that includes the fast mamba solver by default and automatically sources packages from `conda-forge`.

Open a terminal and enable `conda`-related commands by doing the following:
```bash
conda init
```
To use conda commands, you now need to open a new session of your terminal program.

To install the OpenMPI version of Genesis 4, perform the following:

```bash
conda create -n genesis4 genesis4=*=mpi_openmpi*
```

To install the MPICH version of Genesis 4, perform the following:

```bash
conda create -n genesis4 genesis4=*=mpi_mpich*
```

To use Genesis4 from that environment, first activate the environment:

```bash
conda activate genesis4
genesis4 --help
```
### Existing conda installation

Users with Anaconda or another version of conda already installed may simply specify the channel ``conda-forge`` during the environment creation.

To install the OpenMPI version of Genesis 4, perform the following:

```bash
conda create -n genesis4 -c conda-forge genesis4=*=mpi_openmpi*
```

To install the MPICH version of Genesis 4, perform the following:

```bash
conda create -n genesis4 -c conda-forge genesis4=*=mpi_mpich*
```

To use Genesis4 from that environment, first activate the environment:

```bash
conda activate genesis4
genesis4 --help
```
## Compilation 

### macOS 

Installation on macOS requires a suitable compiler and dependencies, which can be provided by [MacPorts](https://www.macports.org). With a working MacPorts, install these:
```bash
sudo port install gcc12
sudo port select gcc mp-gcc12
sudo port install hdf5 +openmpi 
sudo port install fftw-3
```

Then build Genesis:
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Unix
(For detailed instructions for Ubuntu 22.04LTS including a list of required packages, see [here](build_notes/ubuntu_22.04/).)

Genesis supports now the automatic configuration with CMAKE. Following commands will build Genesis from source on Linux platforms. The minimal command to compile the code from the source code root directory is:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Depending on your compiler and system configuration, additional information may be given to CMAKE. Most important is the used compiler. This can be done with the additional definition ```-DCMAKE_CXX_COMPILER=####``` in the cmake command overwriting the default compiler, where ```###``` is the compiler, e.g. `mpicxx` or `CC`. If debugging is needed the CMAKE_BUILD_TARGET should be changed to `Debug`.
The executable is found in the build directory. For a successful build the libraries openmpi/mpich and parallel HDF5 are needed. The FFTW3 is recommended but might become mandatory in upcoming releases.

To complete the installation it is recommended to install the scripts in the directory `sdds2hdf` in a directory (e.g. /bin), which is included in your search path. The files will convert Elegant output distribution into HDF5 format. More info can be found later in the manual.
Also, the directory `xgenesis` includes some functions, which allows Matlab to parse the output file of Genesis and to plot the results. They are not required for running Genesis itself. To use them add the given folder to your Matlab Path or place the file in a directory, which is already in the search path.





## Running

Genesis 1.3 is a command line executable and requires at least a single input argument, which is the filename of the input deck:

```
genesis4 [-o output-rootname] [-l lattice-filename] [-s seed] input-filename
```

The arguments for the lattice filename, the rootname for all output filenames and the seed for the random number generator are optional.
If defined, they act as the default values in the `setup` namelist and thus can be omitted in the main input file.
Any further information must be provided within the input deck. Starting Genesis in this way will effectively invoke a parallel job with only one single core. To start the code with just more than one core the launching command `mpirun` needs to be called, including the number of cores. In its simplest call it would be:

```
mpirun -np xxx genesis4 [-o output-rootname] [-l lattice-filename] [-s seed] input-filename
```

where xxx is the number of cores to be requested.

Once Genesis started, the master core, which is the core with the MPI rank 0, will provide the output to `stdout`. The output can vary for different input decks. Normally it is reported, once the code generates or alters particle distributions or the radiation field, tracks the beam through the lattice or when a time window is specified.

The execution produces output files whenever a tracking command is issued or a field or particle dumped is requested. All in common is the root filename which is defined in the input deck. Normal output files have the extension `.out.h5` for the first tracking or `.Runxxx.out.h` for the xxxth run (except for the first, where the addition of ’Run’ is omitted). The file format is HDF5. The explicit format is described in a later chapter of this manual. Field or particle dumps are using the root filename as well as the default but under certain circumstances this behavior can be overruled. For more information see the corresponding namelist elements in the dedicated chapter.

For time-dependent simulations the code requires a large amount of memories to keep the entire radiation field and particle distribution in memory. This is the significant difference to previous versions where the code was able to progress through the bunch sequentially.
However this excluded features such long range space charge fields, exchanging particles between slices, and large delaying chicanes in the older versions. For an estimate of the memory demand the radiation field requires

```
Nb=ns*ng*ng*16 Bytes,
```

where ns is the number of slices and ng the number of grid points. 

For the electron distribution it is

```
Nb=ns*ne*48 Bytes,
```

with ne the number of macro particles per slice.
Both numbers should be added up and for safety multiplied by a factor of 2 for some overhead (aligning arrays in memory, temporary working arrays, arrays to buffer the output information etc.). The memory request will be distributed  over the nodes. Non time-dependent simulation should be small enough to run on any computer.

In time-dependent runs the number of slices are equally distributed over all cores. The number of slices are even increased to guarantee symmetry, extending the time-window of the simulations. More information about that can be found in the next section. Often it is useful to calculate the number of generated slices in advance, in particular if advanced option such as subharmonic conversion are used. In the case of subharmonic conversion the slices per core should be an integer of the subharmonic number. Otherwise artefacts of numerically enhanced bunching causes wrong results.


<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>
