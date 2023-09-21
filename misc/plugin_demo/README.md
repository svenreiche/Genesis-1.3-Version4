# README
Christoph Lechner, European XFEL
Version: 21 Sept, 2023

Note that this document describes the current developer version. The interface is likely to change in the future.

## Building it

We begin in the directory containing the source code of the demo plugin.
```
mkdir build
cmake ..
make
```

If the build process was successful, your build directory should now contain the shared library `libdemo.so`.

## Including it into simulation
The demo plugin computes the power of the electromagnetic field. The resulting object in the output file should be identical to `/Field/power`. There is a Python script to verify this.

To enable the demo plugin, add the following lines to your `.in` file (before the first `&track` in which the plugin should be used):
```
&add_plugin_fielddiag
    libfile = ./libdemo.so
    obj_prefix = plugindemo
#   interface_verbose = 1
&end
```
Note that currently the full path to the `.so` file containing the compiled plugin code has to be provided. One possible solution is to place a symlink pointing to the location of the `.so` file into the working directory of your simulation.

While setting up the plugin as the tracking begins, there is some informative output:
```
[..]
Running Core Simulation...
Time-dependent run with 2048 slices for a time window of 2.5392 microns
Initial analysis of electron beam and radiation field...
Setting up DiagFieldHook for libfile="./libdemo.so", obj_prefix="plugindemo"
DiagFieldHook::init
Rank 0: Loaded the library, handle is 0x2338c20
Got the library
Rank 0: Factory location is 0x2b4316dd0646
Rank 0: Destroyer location is 0x2b4316dd0694
Rank 0: Calling factory
Rank 0: Calling get_infos
Rank 0: Got class instance
   info_txt="Demo for plugin diagnostics"
   do_multi=0
   provides object named my_power
   provides object named abc
DONE: Registered DiagFieldHook
  Calculation: 0% done
[..]
```

Once the simulation run is finished, the output file contains new objects (output from `h5ls -r name-of-your-file.out.h5`):
```
[..]
/Field/plugindemo        Group
/Field/plugindemo/abc    Dataset {26, 2048}
/Field/plugindemo/my_power Dataset {26, 2048}
/Field/power             Dataset {26, 2048}
[..]
```
