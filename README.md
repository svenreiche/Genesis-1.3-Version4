# GENESIS 1.3

A scientific code to simulate the interaction of electrons within the magnetic field of an undulator to produce coherent radiation, based on the Free-electron LAser process.

News:
Version 4.5.1 has been released (11.1.2022)

Change Log:
[here](CHANGELOG.md)

Manual:
[here](MANUAL.md)

Current Development:
[here](DEVELOPMENT.md)

# Compilation

Genesis supports now the automatic configuration with CMAKE. Following commands will build Genesis from source on Linux platforms
```
mkdir build
chdir build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
If a differnet compiler is needed the additional definition -DCMAKE_CXX_COMPILER=#### in the cmake command overwrites the default compiler. If debugging is needed the CMAKE_BUILD_TARGET should be Debug.

The executable is found in the build directory. Note that for a successful build the libraries openmpi and parallel HDF5 are needed. The FFTW3 is recommended since it might become mandatory in upcoming releases.


