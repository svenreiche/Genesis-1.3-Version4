# GENESIS 1.3 - Output Files

Beside the input and lattice file, Genesis uses exclusively the HDF5 format for the output files, which are the main output file, and the particle and field dump. HDF5 is a well supported format by many programs or programming languages which makes it easy to read in the data. As an example, Matlab requires only a single line to read an entire dataset such as the 2d array for the power, with sample points along the undulator and beam frame, already converted to the right data type. This should make it easy to write your own postprocesser. Nevertheless the Genesis distribution comes with some Matlab routines to ease the parsing and display of the output files.

## Main Output File

The main output file is written when ever the&tracknamelist is called. The name is given by the specified rootname and the extension’.out.h5’for the first call. If there are more than one tracking commands the root is extended by the run number, such as’.Run3.out.h5’ for the third time of tracking.

In the root level of the HDF5 files there are multiple groups, to organize the output data. The group ’Global’ list the basic configuration of the simulation such as reference energy and wavelength or whether it is a time-dependent or scan run. 
The group `Lattice` contains all lattice information. The undulator strength, quadrupole field and other are resolved with the resolution of the requested integration step size, which is given in the dataset `dz`. For the z-position there exist two dataset. The regular onezhas the same length and describes the lattice quantities from the `positionz[i]toz[i]+dz[i]` of the i-th integration step. The dataset `zplot` is used for plotting the beam or field parameters along the undulator. Note that those are evalutated before the integration started, so that there can be one more entry than the lattice datasets. Also if the output is reduced by the output step option in the tracking command, the length of `zplot` is shorter because it has to match the length of the beam and field parameters. 
Another group is `Meta` which contains some additional information, which are not stricktly related to the input case. These are the version number, creation date and user, which has evoked the simulations.

The main data is contained in the `Beam` and `Field` record. In the case that additional harmonics have been selected, there are additional `Field` groups, where the harmonic number has been added. As an example `Field3` would refer to the third harmonic. 
The `Beam` group contains parameters which are evaluated at each integration steps. They are: energy, energy spread, bunching, bunching phase, centroid position in x, y, px and py and the size in x and y. 
Other parameters are only evaluated at the beginning, which are the optical function, emittance and current profile, but in the future they might be becoming also larger, e.g. when sorting is enabled and the current profile can change. In theFieldgroup one can find the power, position and size in x and y as well as the on-axis intensity and phase in the near field (central grid point) or far field (summation over all grid points)

In general the units of the datasets should be given as attributes to the given dataset. This is not yet implemented.

## Particle and Field Dump.

Particle dumps are following the data structure that each slice has its own 6D distribution with the datasets `gamma`, `theta`, `x`, `y`, `px`, and `py`. In additions the local current value is written as well. The slice number is encoded in the group name, which is the composition of `slice` and the 6 digit representation of the slice numbers. Preceding spaces are changed to zeros. E.g. the group name of the 7th slice is `slice000007`. The total number of slices is given by the dataset `slicecount`. In addition some extra information is given on the root level of the reference frame are defined by the starting position `refposition` and the sample length per slice `slicespacing`. The reference length `slicelength` can be used to convert the ponderomotive phase into a longitudinal position in units of meter. Whether the distribution has been generated with beamlets or resolving all electrons are indicated by the parameters `beamletsize` and `one4one`.

The field dump has very similar structure than the beam dump, except that instead of `slicelength` the field is `wavelength` to distinguish between the fundamental or harmonic radiation, while the beam is always measured against the fundamental wavelength for the ponderomotive phase. The field `gridsize` defines the extension of the 2D grid from the origin to the edge along of the major axis. The `slice` groups contains the real and imaginary part of the wavefront, however for programming reason the data record is one dimensional with $ngrid^2$ elements. Any processing or display should convert it to a symmetric 2D array.

## Particle Distribution

Genesis can import particle distribution which are generated by Elegant. However it does not support the parsing of SDDS files directly. Instead, it provides a simple shell script to convert the typical SDDS distribution into a HDF5 file. This script expects that the SDDS toolkit is installed because it extracts the columns `t`, `p`, `x`, `y`, `xp`, and `yp` from the SDDS files and converts them into datasets with the same name in the HDF5 file. The short program will add the extension’.h5’ to the newly created file. This format is the one Genesis expects to be imported. The corresponding charge is supplied in the namelist and therefore does not need to be part of the HDF5 file.


## The Postprocessor XGENESIS

Although the output format of Genesis is HDF5, where entire datasets can be read with a single line in Matlab, the source code distribution also contains some Matlab functions, which help with the postprocessing. They can be found in the directory `xgenesis` from the root directory of the source code. It is recommended that the search path of Matlab extends to that directory or that the files are copied to a locations where the search path already points to.

The following is a brief description of the existing functions. More will come in the future.

## xgeninit.

The routine requires only the name of the output file, which then is parsed for the most essential data, which is the size inzandsand the dataset tree structure. This function needs to be called once. All other functions will make use of the data, which has been parsed here. If another file needs to be read then `xgeninit` needs to be called again.

## xgenplot

The main routine to display data from the output file and returns the data. It can generate multiple plots at the same time and the return data type is a cell array, with each element referring to a dataset in the plot. The dataset itself is a cell array as well, where the first element is the x-coordinate and the second the y.

The function requires up to three input arguments. In special cases, e.g. results from steady-state simulations some options are restricted such as plotting along the time-window coordinate is not meaningful and thus not supported.

The first element mostly refers to the corresponding dataset in the input deck. The function uses regular expression to match the dataset names with the requested plot data. As an example plotting `size` will plot the sizes in x and y of the electron beam and each supported radiation field, so at least 4 plots and even more if harmonics were included. To restrict the plots the requested plot can be more explicit, according to the tree structure in the output file. As an example the input argument `/Beam/.size` would plot the x and y sizes of only the electron beam while `Field.*/.size` does it for the radiation of the fundamental and all harmonics. Besides the datasets in the output file there is one derived dataset, which is the spectrum. The postprocessor can calculate it either from the near field intensity and phase on axis or the equivalent pair in the far field. If only `spectrum` is requested than both are plotted. To select one of them, just specify as it would be done for the intensity, e.g. `spectrum-nearfield`.

If no further arguments are supplied then the first slice of the specified data is plotted along the undulator axis. This is in particular enforced for any lattice function, such as the undulator field. Also if the output data is from a steady-state file then the output is always along undulator axis, despite that a profile might be requested.

The second argument specify the plotting direction: either along undulator direction or along bunch frame or frequency. If the data is a 2D array, e.g. such as the power after a time-dependent run, than plotting a line need as the third argument the specification of the orthogonal position. Here the post processor chooses the line, which lies the closest to the specified point. The default direction and plotting mode is `normal` plotting along the undulator. To plot the power as an example at $3\,\mu m$ of the bunch frame the plot command is

```
xgeninit('power','normal',3e-6)
```.

The direction along the bunch frame is indicated with `profile` such as

```
xgeninit('power','profile',11.5)
```

to plot around the position 11.5 m within the undulator lattice.

For 2D datasets from time-dependent simulations or scansxgenplotallows some condensed plotting along the undulator: `max` for the maximum value of the given profile, `mean` the mean value, and `rms` the rms value. The evolution of the radiation bandwidth is calculated best, using the `rms` option. In fact, RMS values are only meaningful for primary datasets, which are power profile, spectrum and current profile. For all other parameters, such as beam size it is recommended to use the mode `weighted` which is defined as

$$\bar{P}(z)=\frac{\int P(s,z)\cdot W(s,z) ds}{\int W(s,z)ds}$$

where $P$ is the parameter to be display and $W$ the weighting function. For electron beam parameters it is the current while for radiation parameter it is the power profile.

For a two-dimensional output the modes `2d` and `2dnorm` can be used. The latter differs from the former that for each step along the undulator the resulting profile is normalized by its mean value. This allows to exclude the dominant exponential change in some parameters such as power, bunching or energy.

## xgenwigner

This function calculates the 2D Wigner distribution at a position within the undulator, which is supplied as the input argument. The Wigner distribution is a very compact display of time and spectral properties, namely the projection onto the time axis is the power profile and the projection on the frequency axis is the power spectrum.

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>


## Where to Find Further Information

https://gitlab.ethz.ch/Genesis/Genesis4

S. Reiche, ”Numerical Studies for a Single Pass High Gain Free-Electron Laser”, DESY
print, DESY-THESIS-2000-012 (2000)

### FEL Basics
- E.L. Saldin, E.A. Schneidmiller and M.V.Yurkov, ”The Physics of Free Electron Lasers”, Springer, (2000).
- J.B. Murphy and C. Pellegrini, ”Introduction to the Physics of the FEL”, Proceedings of the South Padre Island Conference, Springer (1986) 163. 

### Numeric Basics
- W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery, ”Numerical Recipes in FORTRAN”, Cambridge University Press, (1988)

- W.F. Ames, ”Numerical Methods for Partial Differential Equations”, Academic Press (1996)

### Other software refences

- C++: B. Stroustrup, ”C++ Programming Language”, Addison-Wesley (1997).

- [OpenMPI](https://www.open-mpi.org/)

- [HDF5](https://support.hdfgroup.org/HDF5/)

