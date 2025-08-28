# README
Christoph Lechner, European XFEL
Version: April 2024

## Overview
In this directory you find the source code of a test environment for field analysis plugins (currently, plugins for analysis of the electron beam are not supported). This test environment allows to perform tests without having to run a complete "GENESIS 1.3" simulation run.

## Building it
To build follow the usual procedure to compile "GENESIS 1.3", but with the additional parameter `-DUSE_DPI=ON`:
```
% cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_DPI=ON ..
% make
```
This compiles "GENESIS 1.3", v4 and the plugin testbench which re-uses some of the code of "GENESIS 1.3", v4 to perform the tests.

## Basic usage
### Single process
The parameters are specified in a configuration file. This directory contains one example configuration file, `demo_params.txt`. The testbench expects a single argument, the name of the configuration file to use.
```
./g4_tb demo_params.txt
```
At the end of the execution, a file named `plugin_data.txt` is generated (at the moment the name of this file is hardcoded). The formatting of the file is optimized for import into Python/numpy.

### In MPI environment
FEL simulations with "GENESIS 1.3", v4 typically run in an MPI environment. The testbench also supports operation in MPI enviroments.

Depending on the configuration of your MPI environment this can as simple as 
```
mpirun ./g4_tb demo_params.txt
```

Note that in this case the parameter `nslice` in the configuration file specifies the number of slices per process on the MPI communicator. The diagnostics data written to the generated output file contains the data generated in all processes, but only the data generated on MPI rank 0 is written to stdout.

## Usage with gdb
Now, we illustrate how to use the testbench to debug a plugin with `gdb`. This example was done with git commit ID 3295805.

While the testbench supports MPI, we don't use it here for sake of simplicity. One standard approach to debug the field analysis code is to place a breakpoint at the beginning of the field processing code in the plugin:
```
% gdb ./g4_tb
GNU gdb (GDB) 9.2
Copyright (C) 2020 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
Type "show copying" and "show warranty" for details.
This GDB was configured as "x86_64-pc-linux-gnu".
Type "show configuration" for configuration details.
For bug reporting instructions, please see:
<http://www.gnu.org/software/gdb/bugs/>.
Find the GDB manual and other documentation resources online at:
    <http://www.gnu.org/software/gdb/documentation/>.

For help, type "help".
Type "apropos word" to search for commands related to "word"...
Reading symbols from ./g4_tb...
(gdb) b DiagFieldHookedDemo::doit
Function "DiagFieldHookedDemo::doit" not defined.
Make breakpoint pending on future shared library load? (y or [n]) y
Breakpoint 1 (DiagFieldHookedDemo::doit) pending.
(gdb) r demo_params.txt
Starting program: .../build_20240328/misc/plugin_testbench/g4_tb demo_params.txt
[Thread debugging using libthread_db enabled]
Using host libthread_db library "/lib64/libthread_db.so.1".
[Detaching after fork from child process 75512]
[New Thread 0x2aaab0691700 (LWP 75525)]
[New Thread 0x2aaab0a55700 (LWP 75533)]
opened config file 'demo_params.txt.
processed config file.

Setting up DiagFieldHook for libfile="./libdemo.so", obj_prefix="plugin"
DiagFieldHook::init
Rank 0: Loaded the library, handle is 0x6abf50
Got the library
Rank 0: FieldDiag Factory location is 0x2aaab31af686
Rank 0: FieldDiag Destroyer location is 0x2aaab31af6d4
Rank 0: BeamDiag Factory location is 0
Rank 0: BeamDiag Destroyer location is 0
Rank 0: Calling factory
Rank 0: Calling get_infos
Rank 0: Got class instance
   info_txt="Demo for plugin diagnostics"
   do_multi=0
   provides object named my_power
   provides object named abc
DONE: Registered DiagFieldHook

Thread 1 "g4_tb" hit Breakpoint 1, DiagFieldHookedDemo::doit (this=0x6ac620, pd=0x7fffffffccf0) at .../misc/plugin_demo/DiagFieldPowerDemo.cc:46
46              verify_datastructure(pd);
(gdb)
```
The execution stops inside of the plugin being studied. From here on you can use standard debugging techniques.
