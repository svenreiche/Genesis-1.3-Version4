# GENESIS 1.3 - Main Input Deck

## Main Input File

Unlike older version, which supported only a single namelist as input, the latest version splits up the namelist and processes them sequentially as they appear in the file. Therefore the order is important, e.g. you cannot define a time window after the electron beam or the field has been already generated. On the other hand the code supports multiple tracking in one run.
For that the particle and field distribution is kept in memory and then reused in the next tracking. If a fresh bunch is required the old has to be discarded and then a new one has to be generated. In that sense it is possible to combine steady-state and time-dependent run in one input deck by the sequence: generate field, generate beam, track, define time window, generate beam, generate field, and track.

The namelists themselves always start with the ampersand `&` and the name of the name list and are terminated with `&end`. Each line consists out of an assignment of a supported variable, e.g `npart = 2048` in the namelist `&setup`. Empty lines and lines, starting with the pound symbol `#` are ignored. If the same variable occurs more than once the latest occurrence is used for the calculation. If Genesis encounters a variable name, which it does not recognize it stops execution and prints out a list of the supported variables.

The following describes all supported namelist with their variables, including its required format and default value.

## Supported Namelists
  - [setup](#setup)
  - [alter_setup](#alter_setup)
  - [lattice](#lattice)
  - [time](#time)
  - [profiles](#profiles)
    - [profile_const](#profile_const)
    - [profile_gauss](#profile_gauss)
    - [profile_step](#profile_step)
    - [profile_polynom](#profile_polynom)
    - [profile_file](#profile_file)
  - [sequence](#sequences)
    - [sequence_const](#sequence_const)
    - [sequence_polynom](#sequence_polynom)
    - [sequence_power](#psequence_power)
    - [sequence_random](#profile_random)
  - [beam](#beam)
  - [field](#field)
  - [importdistribution](#importdistribution)
  - [importbeam](#importbeam)
  - [importfield](#importfield)
  - [importtransformation](#importtransformation)
  - [efield](#efield)
  - [sponrad](#sponrad)
  - [wake](#wake)
  - [sort](#sort)
  - [write](#write)
  - [track](#track)

[Back](#supported-namelists)
<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### setup

The namelist `setup` is a mandatory namelist and should be the first in the input deck. It contains the basic parameters to control the simulations. It can only be called once. If the user want to change some parameter the namelist `alter_setup` should be used.

- `rootname` (*string, \<empty>*): The basic string, with which all output files will start, unless the output filename is directly overwritten (see `write` namelist)
- `outputdir` (*string, \<empty>*): Output directory name.
- `lattice` (*string, \<empty>*): The name of the file which contains the undulator lattice description. This can also include some relative paths if the lattice file is not in the same directory as the input file.
- `beamline` (*string, \<empty>*): The name of the beamline, which has to be defined within the lattice file. For more information on the lattice file, see the next chapter.
- `gamma0` (*double, 11350.3*): The reference energy in unites of the electron rest mass. This is the reference energy which is used in the code at various place, mostly in the calculation of the matching condition, the reference focusing strength of quadrupoles and undulator as well as the default value if an electron distribution is generated.
- `lambda0` (*double, 1e-10*): The reference wavelength in meter, which is used as the wavelength in steady-state simulation or for defining the sample distance in time-dependent runs. It also acts as the default value when field distributions are generated.
- `delz` (*double, 0.015*): Preferred integration stepsize in meter. Note that this is not a strict value because Genesis tries to optimized the stepsize according to the elements it can resolve. E.g. if an undulator is 1.99 m long but the preferred stepsize is 2 cm than it uses a stepsize which is the closest to preserve the number of integration step. In this case the preferred stepsize gives 99.5 steps which is than rounded to 100 and thus resulting in an actual stepsize of 1.99 cm. Note that outside of the undulator, which are free drifts for the radiation field, Genesis progresses the electron beam and radiation field in larger steps, namely one step per resolved element (drift, quadrupole, phase shifter).
- `seed` (*int, 123456789*): Seed to initialize the random number generator, which is used for shot noise calculation and undulator lattice errors, though it is recommended that the random number generator seed is redefined explicitly for undulator errors in its corresponding namelist.
- `npart` (*int, 8192*): Number of macro particles per slice. Note that the number must be a multiple of the used bins `nbins` otherwise Genesis will exit with an error. If one-for-one simulations are used, this parameter has no meaning.
- `nbins` (*int, 4*): Number of macro particles, which are grouped into beamlets for gener ating the correct shot noise. For one-for-one simulations this parameter has no meaning
- `one4one` (*bool, false*): Flag to enable or disable resolving each electron in the simulation. This is mandatory for certain features, such as sorting or slicing of particle distributions. If set to `true` other parameters such as `npart` and `nbins` are obsolete and do not need to be defined. It is recommended to estimate the number of electrons, which are generated in the simulations, because this can easily required memory beyond what is available on the computer.
- `shotnoise` (*bool, true*): Flag to enable the calculation of shotnoise per each slice during generation of the electron distribution. It is recommended to set the value to `false` for steady-state or scan simulations.
- `beam_global_stat` (*bool, false*): Flag to enable extra output of beam parameters of the entire bunch, such as energy, energy spread etc. The data are placed in the HDF group ”Global” within the group ”Beam” of the output file
- `field_global_stat` (*bool, false*): Flag for the field output, similar to `beam_global_stat`.
- `exclude_spatial_output` (*bool, false*): Flag to suppress the datasets in the output file for the x- and y-position and size (both Beam and Field) and px- and py-position (Beam only). This might be useful to reduce the file size of the output file, if these datasets are not needed for the post-processing
- `exclude_fft_output` (*bool, false*): Flag to suppress the datasets in the output file for the field divergence and pointing. Since it also disable the FFT calculation of the 2D wavefronts it speeds up the execution time slightly. If the code has been compiled without the support of the FFTW library this parametr has no effect.
- `exclude_intensity_output` (*bool, false*): Flag to suppress the datasets for the near and farfield intensity and phase for the radiation field. If excluded the output file size becomes smaller but no post-processing calculation of the spectra is possible.
- `exclude_energy_output` (*bool, false*): Flag to suppress the datasets in the output file for the mean energy and energy spread of the electron beam.
- `exclude_aux_output` (*bool, false*): Flag to suppress the auxiliary datasets in the output file. In the moment it is the long-range longitudinal electric field as seen by the electrons.
- `exclude_current_output` (*bool, true*): Flag to reduce the size of the current dataset for the electron beam. Under most circumstances the current profile is constant and only the initial current profile is written out. However, simulation with one-4-one set to `true` and sorting events the current profile might change. Example are ESASE/HGHG schemes. By setting the flag to false the current profile is written out at each output step similar to radiation power and bunching profile.
- `exclude_field_dump` (*bool, false*): Exclude the field dump to `.fld.h5`.
- `write_meta_file` (*bool, false*): Write a metadata file.
- `semaphore_file_name` (*string, <derived from rootname>*): Providing a file name for the semaphore file always switches on writing the "done" semaphore file, overriding 'write_semaphore_file' flag. This allows to switch on semaphore functionality just by specifying corresponding command line argument -- no modification of G4 input file needed.
- `write_semaphore_file` (*bool, false*): Write a semaphore file when the simulation has completed.
- `write_semaphore_file_done` (*bool, false*): Alias for `write_semaphore_file`. This takes precedence over `write_semaphore_file` if both are specified.
- `write_semaphore_file_started` (*bool, false*): Write a semaphore file at startup, after the setup block is parsed.

[Back](#supported-namelists)


<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### alter_setup

A namelist to change some parameters within the simulation, which have been defined alread by the `setup`-namelist. The change values are stored in the setup module so that for another invocation of alter_setup some defaults values are use which have been defined in the preceding call of alter_setup

- `rootname` (*string, <taken from SETUP>*): The basic string, with which all output files will start, unless the output filename is directly overwritten (see `write`-namelist)
- `beamline` (*string, \<empty>*): The name of the beamline, which has to be defined within the lattice file. This way another beamline can be selected in the case the simulation has multiple stages
- `delz` (*double, <taken from SETUP>*): Preferred integration stepsize in meter. Note that this is not a strict value because Genesis tries to optimized the stepsize according to the elements it can resolve. E.g. if an undulator is 1.99 m long but the preferred stepsize is 2 cm than it uses a stepsize which is the closes to preserve the number of integration step. In this case the preferred stepsize gives 99.5 steps which is than rounded to 100 and thus resulting in an actual stepsize of 1.99 cm. Note that outside of the undulator Genesis, which are free drifts for the radiation field, it progress the electron beam and radiation field in larger steps, namely one step per resolved element (drift, quadrupole, phase shifter).
- `harmonic` (*int, 1*): If the value is not 1 than a harmonic conversion is done. This has several consequences. The reference wavelength in `setup` is divided by the harmonic number, the sample rate in `time` is multiplied by the harmonic number, the ponderomotive phases of all macro particles are scaled with the harmonic number, all radiation fields, which are not identical to the harmonic numbers are deleted, while an existing harmonic field is changed to be at the fundamental wavelength
- `subharmonic` (*int, 1*): If the value is not 1 than a down conversion is done. It is similar to the action of `harmonics` but in the opposite directions. For the radiation field all field definitions are deleted except for the fundamental, which is converted to a harmonic. In this case the fundamental field needs to be defined before another tracking is called.
- `resample` (*bool, false*): If this is set to true and only if one-for-one simulations are used the harmonic and subharmonic conversion can re-sample to the new wavelength. In the case of up-conversion the slices are split and the total number of slices increases. Same with the radiation field. An previously existing harmonic field, which is now becoming the fundamental, is interpolated between the existing sample points (still needs to be implemented). If a new field is generated it has automatically the new number of slices. If also prevents that the sample rate is changed by remaining unchanged.
- `disable` (*bool, false*): Disable non-matching radiation harmonic.


[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### lattice

This namelist is used to change the raw lattice from the lattice file, such as generating errors in the position of the elements. The namelist can be defined several times to add more than one error source to the lattice.

- `zmatch` (*double, 0*): If the position within the undulator in meter is non-zero than Genesis tries to calculate the matched optics function for a periodic solution. In the case that it cannot find a solution than it will report it. Found solution will also be the default values for a succeeding beam generation, so that no explicit optical functions need to be defined any longer. If the lattice is highly non-periodic it is recommended
    to find the matching condition with an external program such as MAdX.
- `element` (*string, \<empty>*): Name of the element type, which will be changed, e.g. Undulator if undulator modules are altered. Only the first 4 letters need to be defined. If there is no match, e.g. due to a type, nothing will be changed. It acts rather as a filter than a mandatory element. Elements of the type `MARKER` are not supported.
- `field` (*string, \<empty>*): attribute name for a given element. The names are the same as in the definition of the lattice file. The field acts as a filter again. With non-matching events nothing will be changed.
- `value` (*double, 0, or sequence label*): The new value. If a reference to a sequence is used, values can be different depending on how many elements are changed. For a double the value would be the same for all elements affected.
- `instance` (*integer, 0*): The instances of affected elements. If a positive value is given, than only that element is changed, where its occurence matches the number. E.g. for a value of 3 only the third element is selected. For a value of 0 all elements are changed. The ability to change more than one but less than all is currently not supported.
- `add` (*bool, true*): If the value is `true`, the changes are added to the existing value. For a value of `false`, the old values are overwritten.
- `resolvePeriod` (*bool, false*): currently unused.

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### time

This namelist defines the time window/range for simulation with more than just one slice.
For reference the complementary axis of the undulator axis, which is normally the position in the time frame, is expressed in a position `s`. Normally everything is aligned to the origins = 0, in particular when external distributions are imported. Note that for parallel execution the number of slices per core must be the same for an efficient writing of the output files. Therefore Genesis extends the time-window to symmetrize the number of slices per core by extending it towards larger values of `s`.
As an example, with `XLAMDS=1e-6` and a length `SLEN=20e-6` a call of Genesis with 24 cores would generate a time-window of 24 microns because each core would have one slice, while 15 cores would expand it to 30 microns with 2 slices per core each.

This module defines also scans in either field or beam parameters if the corresponding flag is set. Technically it generates the beam and field as for time-dependence but disables slippage during simulations. That way the radiation field is kept in the same slice, acting as steady-state simulations.

- `s0` (*double, 0*): Starting point of the time-window in meters.
- `slen` (*double, 0*): Length of the time window in meters. Note that for parallel jobs this might be adjusted towards larger values.
- `sample` (*int, 1*): Sample rate in units of the reference wavelength from thesetup namelist, so that the number of slices is given by `SLEN / LAMBDA0 /SAMPLE` after `SLEN` has been adjusted to fit the MPI size.
- `time` (*bool, true*): Flag to indicate time-dependent run. Note that time-dependent simulations are enabled already by using this namelist. This flag has the functionality to differentiate between time-dependent run and scans, which disable the slippage in the tracking. To restrict the simulation to steady-state the `time` namelist has to be omitted from the input deck.

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### profiles

Profiles are defining a dependence on the position in the time frame, which then can be used to defined the behavior of various parameters such as slice emittance, energy spread etc. All profiles are in common that they have to have a label with which they are referred to in other namelists. To indicate the reference to a profile and thus allows to distinguish Genesis between normal numbers the label name must have the additional character ’@’ in front. E.g. if the name of the label is `energy` then in the beam name list it is referred to by `gamma = @energy`.

#### profile_const

- `label` (*string, \<empty>*): Name of the profile, which is used to refer to it in later calls of namelists
- `c0` (*double, 0*): constant value to be used.

#### profile_gauss

- `label`(*string, \<empty>*): Name of the profile, which is used to refer to it in later calls of namelists
- `c0` (*double, 0*): Maximum function value of the Gaussian distribution
- `s0` (*double, 0*): Center point of the Gaussian distribution
- `sig` (*double, 0*): Standard deviation of the Gaussian distribution

#### profile_step

- `label` (*string, \<empty>*): Name of the profile, which is used to refer to it in later calls of namelists
- `c0` (*double, 0*): Constant term
- `s_start` (*double, 0*): Starting point of the step function
- `s_end` (*double, 0*): Ending point of the step function

#### profile_polynom

- `label` (*string, \<empty>*): Name of the profile, which is used to refer to it in later calls of namelists
- `c0`(*double, 0*): Constant term
- `c1`(*double, 0*): Term proportional to s
- `c2`(*double, 0*): Term proportional to s^2
- `c3`(*double, 0*): Term proportional to s^3
- `c4`(*double, 0*): Term proportional to s^4

#### profile_file

- `label`(*string, \<empty>*): Name of the profile, which is used to refer to it in later calls of namelists
- `xdata` (*string, \<empty>*): Points to a dataset in an HDF5 file to define the `s`-position for the look-up table. The format is `filename/group1/.../groupn/datasetname`, where the naming of groups is not required if the dataset is at root level of the HDF file
- `ydata` (*string, \<empty>*): Same as y data but for the function values of the look-up table.
- `isTime` (*bool, false*): If true the `s`-position is a time variable and therefore multiplied with the speed of light `c` to get the position in meters.
- `reverse`(*bool, false*): if true the order in the look-up table is reverse. This is sometimes needed because time and spatial coordinates differ sometimes by a minus sign.
- `autoassign`(*bool, false*): use the HDF5 file from `xdata` (TODO more details).

[Back](#supported-namelists)


<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### sequences

Sequences act similar to profiles that they define a series of values which can be used to control dependence of certain parameters. While for profiles these are mostly properties of the electorn beam and radiation field, sequences are used for
occurences of beamline elements for specific parameters. Examples are applying a taper gradiant to the undulator offset or add random offsets to the quadrupole positions.
All sequences have in common that they have to have a label with which they are referred to in the lattice.
To indicate the reference to a sequence and thus allows to distinguish Genesis between normal numbers
the label name must have the additional character ’@’ in front. E.g.
if the name of the label is `taper` then in the parameter list of the modulator element
it is referred to by `aw = @taper`. A sequence is evoked for each occurence of the reference. As an example
in a line which is 6 times a Fodo Cell which each space between the quadrupoles hosting an undulator, a reference by the undulator field will evoke the first 12 elements of the sequence.

#### sequence_const

- `label` (*string, \<empty>*): Name of the sequence, which is used to refer to it in the lattice
- `c0` (*double, 0*): constant value to be used.

#### sequence_polynom

- `label` (*string, \<empty>*): Name of the sequence, which is used to refer to it in the lattice
- `c0`(*double, 0*): Constant term
- `c1`(*double, 0*): Term proportional to s
- `c2`(*double, 0*): Term proportional to s^2
- `c3`(*double, 0*): Term proportional to s^3
- `c4`(*double, 0*): Term proportional to s^4

#### sequence_power

- `label`(*string, \<empty>*): Name of the sequence, which is used to refer to it in the lattice
- `c0` (*double, 0*): Constant term
- `dc` (*double, 0*): Term scaling the growing power series before added to the constant term
- `alpha` (*double, 0*): power of the series
- `n0` (*integer, 1*): starting index of power growth. Otherwise the sequence uses only the constant term

#### sequence_random

- `label`(*string, \<empty>*): Name of the sequence, which is used to refer to it in the lattice
- `c0` (*double, 0*): Mean value
- `dc` (*double, 0*): Amplitude of the error, either the standard division for normal distribution or the min and max value for uniform distribution.
- `seed` (*integer, 100*): seed for the random number generator
- `normal` (*bool, true*): Flag for Gaussian distribution. If set to false a uniform distribution is used.


[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### beam

This namelist initiates the generation of the particle distribution to be kept in memory. Any time-dependence has to be defined before calling this namelist.

- `gamma` (*double, gamma0 or profile label*): Mean energy in units of the electron rest mass. If default value is given by the reference energy from the `setup`-namelist.
- `delgam` (*double, 0 or profile label*): RMS energy spread in units of the electron rest mass.
- `current` (*double, 1000 or profile label*): Current in Amperes.
- `ex` (*double, 0.3e-6 or profile label*): Normalized emittance in $x$ in units of meters
- `ey` (*double, 0.3e-6 or profile label*): Normalized emittance in $y$ in units of meters
- `betax` (*double, 15 or matched value or profile label*): Initial beta-function in $x$ in meters. If the matched command has been invoked before the default values are set to the results.
- `betay` (*double, 15 or matched value or profile label*): Initial beta-function in $y$ in meters. If the matched command has been invoked before the default values are set to the results.
- `alphax` (*double, 0 or matched value or profile label*): Initial alpha-function in $x$. If the matched command has been invoked before the default values are set to the results.
- `alphay` (*double, 0 or matched value or profile label*): Initial alpha-function in $y$. If the matched command has been invoked before the default values are set to the results.
- `xcenter` (*double, 0 or profile label*): Initial centroid position in $x$ in meter.
- `ycenter` (*double, 0 or profile label*): Initial centroid position in $y$ in meter.
- `pxcenter` (*double, 0 or profile label*): Initial centroid momentum in $x$ in units of $\gamma \beta_x$.
- `pycenter` (*double, 0 or profile label*): Initial centroid momentum in $y$ in units $\gamma \beta_y$.
- `bunch` (*double, 0 or profile label*): Initial bunching value
- `bunchphase` (*double, 0 or profile label*): Initial phase of the bunching
- `emod` (*double, 0 or profile label*): Initial energy modulation in units of the electron rest mass. This modulation is on the scale of the reference wavelength
- `emodphase` (*double, 0 or profile label*): Initial phase of the energy modulation

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### field

This namelist initiate the generation of the field distribution. It differs in one point from the generation of the beam. It can be called multiple times. If the variable `accumulate` is set to true, it does not delete the previous distribution but adds up the wavefronts. That way higher mode content in either spatial and time direction can be created.

- `lambda` (*double, lambda0*): Central frequency of the radiation mode. The default value is the reference wavelength from the `setup`-namelist.
- `power` (*double, 0 or profile label*): Radiation power in Watts
- `phase` (*double, 0 or profile label*): radiation phase in rads. Note that a linear profile results in a shift in the radiation wavelength, which is also the method if for the variable `lambda` a different value than the reference wavelength is used. In case of conflicts the profile for the phase definition has priority.
- `waist_pos` (*double, 0 or profile label*): Position where the focal point is located relative to the undulator entrance. Negative values place it before, resulting in a diverging radiation field.
- `waist_size` (*double, 100e-9 or profile label*): Waist size according to the definition of $w_0$ according to Siegman’s ’Laser’ handbook
- `xcenter` (*double, 0*): Center position in $x$ in meter of the Gauss-Hermite mode
- `ycenter` (*double, 0*): Center position in $y$ in meter of the Gauss-Hermite mode
- `xangle` (*double, 0*): Injection angle in $x$ in rad of the Gauss-Hermite mode
- `yangle` (*double, 0*): Injection angle in $y$ in rad of the Gauss-Hermite mode
- `dgrid` (*double, 1e-3 or by existing field*): Grid extension from the center to one edge. The whole grid is twice as large with 0 as the center position
- `ngrid` (*int, 151 or by existing field*): Number of grid points in one dimension. This value should be odd to enforce a grid point directly on axis. Otherwise the convergence in the simulations could be worse.
- `harm` (*int, 1*): Harmonic number of the radiation field with respect to the reference wavelength.
- `nx` (*int, 0*): Mode number in $x$ of the Gauss-Hermite mode
- `ny` (*int, 0*): Mode number in $y$ of the Gauss-Hermite mode
- `accumulate` (*bool, false*): If set the generated field is added to an existing field instead of overwriting it.

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### importdistribution

This namelist controls the import of an external distribution which are generated from Elegant. The file has to be in HDF5 format. In the distribution is a shell script to convert an Elegant sdds-output file into the HDF5 format. The distribution has to provide all 6 dimensions while the charge is supplied in this namelist. When imported the longitudinal position is changed so that the last particles is at $s=0$ micron.

Genesis will import the distribution and evaluate a part of the distribution in terms of emittance, twiss function etc. The part is defined by the parameters `eval_start` and `eval_end`. They are given in relative terms with respect to the total length of the distribution. E.g. values of 0.25 and 0.75 for these two parameters, respectively will use the central 50% of the bunch to determine the twiss function.

The beam can be match to new twiss parameters defined explicitly in the input deck or by the default values. Note that the matching region is defined by the evaluation range. So it is possible to match a certain slice in the distribution and the rest of the distribution will follow that distribution.

The flag `settimewindow` will update the length of the timewindow. But it is required that a time window is enabled in the input deck first. Also in this case import distribution needs to be called before the field definition.

One important parameter is slicewidth, with is the relative length of a slice in the distribution, for which the contained particles are used to recreate the internal distribution. It is a rolling window averaging of the particle distribution. Here a compromise has to be found between sufficient resolution to retrieve the time-dependent variation in current, emittance etc and having sufficient particles for a valid reconstruction.
A crude rule of thumb is that in average there should be at least 1000 particle for the recreation. As an example: a distribution of 200000 particles would allow for 200 slices with about 1000 particles in it.
These 200 slices correspond to a `slicewidth`  of 0.005.


- `file` (*string, \<empty>*): The file name of the distribution, including possible relative directories.
- `charge` (*double, 0*): Total charge of the distribution to calculate the current and individual charge per macro particle.
- `slicewidth` (*double, 0.01*): the fraction in length of the distribution which is used for reconstruction. E.g if the length is 10 micron and slic ewidth 0.02 then the reconstruction at the positions $s= 4\,\mu m$ is using those particles in the distribution, which are located in the slice from $3.9\, \mu m$ to $4.1\,\mu m$.
- `center` (*bool, false*): If set to true the particle distribution is recentered in transverse position, momenta and energy.
- `gamma0` (*double, gamma0 from setup*): If centering is enabled, new center in energy in units of electron rest mass.
- `x0` (*double, 0*): If centering is enabled, new center in $x$ in meter.
- `y0` (*double, 0*): If centering is enabled, new center in $y$ in meter.
- `px0` (*double, 0*): If centering is enabled, new mean momentum in $x$ in $\gamma \beta_x$.
- `py0` (*double, 0*): If centering is enabled, new mean momentum in y in $\gamma \beta_y$.
- `match` (*bool, false*): If set to `true`, the particle distribution is matched to new optical function values.
- `betax` (*double, 15 or matched value*): If matching is enabled, new beta function in $x$ in meters.
- `betay` (*double, 15 or matched value*): If matching is enabled, new beta function in $y$ in meters.
- `alphax`(*double, 0 or matched value*): If matching is enabled, new alpha function in $x$.
- `alphay`(*double, 0 or matched value*): If matching is enabled, new alpha function in $y$.
- `eval_start` (*double, 0*): evaluation start.
- `eval_end` (*double, 1*): evaluation end.
- `settimewindow` (*bool, false*): set time window.
- `align` (*int, 0*): currently unused.
- `align_start` (*double, 0*): currently unused.
- `align_end` (*double, 1*): currently unused.


This HDF5 file `file` should have the following datasets:
- `x` is the horizontal coordinate in meters
- `y` is the vertical coordinate in meters
- `xp` = $p_x$/$p_z$ is the dimensionless trace space horizontal momentum
- `yp` = $p_y$/$p_z$ is the dimensionless trace space vertical momentum
- `t` is the time in seconds
- `p`  = relativistic $\gamma*\beta$, which is the total momentum divided by $m_ec$


[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### importbeam

The modules controls the import of a Genesis 1.3 particle file to replace the internal generation of the particle distribution (note that the module `beam` should not be called). The routine defines also the parameter for a time-dependent run if the `time`-namelist hasn’t been defined yet.

- `file` (*string, \<empty>*): File name of a hdf5 complient datafile to contain the slice-wise particle distribution. It has to follow the internal Genesis 1.3 syntax.
- `time` (*bool, `true`*): If the time window hasn’t be defined it allows to run Genesis with the imported distribution in scan mode, when set to `false`. This would disable all slippage and long-range collective effects in the simulation

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### importfield

The modules controls the import of a Genesis 1.3 field file to replace the internal generation of the field distribution (note that the module `field` should only be called afterwards with the `accumulate`-option enabled). The routine defines also the parameter for a time-dependent run if the `time`-namelist hasn’t been defined yet.

- `file` (*string, \<empty>*): File name of a hdf5 compliant datafile to contain the slice-wise particle distribution. It has to follow the internal Genesis 1.3 syntax.
- `harmonic` (*int, 1*) defines the harmonic for the given Genesis run.
- `time` (*bool, true*): If the time window hasn’t be defined it allows to run Genesis with the imported distribution in scan mode, when set to `false`. This would disable all slippage and long-range collective effects in the simulation

[Back](#supported-namelists)


### importtransformation

Once an electron distribution is generated the namelist can be used to manipulate the distribution by shifting the particle by the vector dr or applying the transport matrix R. The applied transformation is `r1 = R*r0+dr`, where
`r0` is the initial particle vector and `r1` the final one. The transformation assumes the standard 6D vector of `(x,x',y,y',s,delta)`.
the supplied vector and matrix must have the corresponding shape (6 or 6x6).
The user can supply more than one vector or matrix, e.g. sampling at various positions `s`. Then the transformation used interpolated values.
Note that in the case of transport matrices and interpolated matrix does not preserve the emittance. In this case a high sample rate should be supplied to reduce this effect to a minimum.
Genesis will check the shape of the transport vector and matrices. If the rank is higher than needed (e.g. 2x6x6 for a transport matrix) then it assumes the first index
refers to the sample along the `s`-axis. In this case the sample distance `slen` should be also specified. In the case that `n=1` or `slen=0` only a global transformation is applied.

- `file` (*string, \<empty>*): File name of a hdf5 compliant datafile to contain the vector and matrix informations
- `vector` (*string, \<empty>*) Name of the dataset which contains the vector information. The shape must be either (6) or (n,6)
- `matrix` (*string, \<empty>*) Name of the dataset which contains the matrix information. The shape must be either (6,6) or (n,6,6)
- `slen` (*double, 0*): The length in meters between adjacent sample points (n>1), needed for the interpolation. If the value is zero only a global transformation is applied using the first entry.
[Back](#supported-namelists)



<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### efield

This namelist controls the long and short range space charge fields. The long range corresponds to any length scale longer than the slice length of the simulation, while the short range is on the resonant wavelength scale. Numerically they are treated differently.
The calculation for the short range is done on a radial-azimuthal grid, centered to the centroid position of the electron slice, while the long range is the sum of the space charge field in the rest frame where each slice is treated as a uniform disk.

- `longrange` (*bool, false*): Flag to enable the calculation of the long range space charge field.
- `rmax` (*double, 0*): Size of radial grid in meters. If the beam size gets larger than the grid the size is automatically adjusted to the maximum radius of the electrons with an additional 50% extension. When the mesh size is adjusted a message will be printed on screen.
- `nz` (*int, 0*): Number of longitudinal Fourier component of the short range space charge field. Note that this should be not in conflict with the beamlet size.
- `nphi` (*int, 0*): Number of azimuthal modes in the calculation of the short range space charge field.
- `ngrid` (*int, 100*): Number of grid points of the radial grid for the short range space charge field.

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### sponrad

This enables the effect of spontaneous radiation outside of the frequency band of the FEL simulation.

- `seed` (*int, 1234*): Seed for random number generator to model the quantum fluctuation of hard photons.
- `doLoss` (*bool, false*): If set to `true`, electrons will loose energy due to the emission of spontaneous radiation within the undulator
- `doSpread` (*bool, false*): If set to `true`, the energy spread will increase due to the fluctuation in the emission of hard photons of the spontaneous radiation.

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### wake

Genesis supports the calculation of three types of wakefields by specifying the typical input parameters (e.g. gap length for the geometric wakefield). It first solves the single particle wake and then convolutes with the current distribution. Therefore it follows the change in the wakepotential if a chirped beams undergoes a compression in a chicane. In addition an external loss factor can be supplied, which can also refer to a profile. In this case it is treated as the full wake and subtracted from the particle energy directly.

*Note that this functionality hasn't been fully tested yet or optimized for rapid calculation*

- `loss` (*double, 0 or profile label*): Loss in $eV/m$. This is a global loss function (in particular if a profile is defined). Its function values V(s) remains unchanged even if the current profile changes
- `radius` (*double, 2.5e-3*): Radius of the aperture if it is a round chanber or half the distance in the case of two parallel plates.
- `roundpipe` (*bool, true*): Flag to indicate the shape of the transverse cross-section of the aperture. If set to `true`, a round aperture is assumed, otherwise the model has two parallel plates.
- `conductivity` (*double, 0*): Conductivity of the vacuum material for the resistive wall wakefield function
- `relaxation` (*double, 0*): Relaxation distance (aka the mean free path of the electron in the vacuum material) for the resistive wall wakefields
- `material` (*string, \<empty>*): String literal to define conductivity and relaxation distance for either copper or aluminum by using the two character label ’CU’ or ’AL’ repectively. This overwrites also any explicit definition of the conductivity and relaxation value.
- `gap` (*double, 0*): Length in mm of a longitudinal gap in the aperture, exciting geometric wakes.
- `lgap` (*double, 1.0*): Effective length over which a single gap is applied. E.g. if there is a periodicity of 4.5 m at which there is always the same gap in the aperture for the geometrice wakes, then this value should be put to 4.5 m.
- `hrough` (*double, 0*): Amplitude in meters of a sinusoidal corrugation, modeling the effect of surface roughness wakes.
- `lrough` (*double, 1*): period lengthin meters of the sinusoidal corrugation of the surface roughness model.
- `transient` (*bool, false*): If set to `true`, Genesis includes the catch-up length of the origin of the wakefield to the particle effects. E.g. particles do not see immediatly the wake from those closer ahead of them than those further away. The catch-up distance is the distance in the undulator added to the starting position `ztrans`. If set to false the steady-state model is used, effectively setting `ztrans` to infinity. Enabling transient calculation will update the wakefield at each integration step, which can slow down the calculations.
- `ztrans` (*double, 0*): Reference location of the first source of the wake fields. A positive value means that the condition for wakes (e.g. a small aperture in the vacuum chamber) has already started and there has been already some length to establish the wakes. For a value of zero the source is right at the undulator start, while a negative value prevents any wake, till the interation position has passed that point.

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### sort

An empty namelist with no variables. It initiates the sorting and redistribution of particles only if one-for-one simulation is enabled. Note that harmonic conversion will automatically invokes sorting and therefore does not need to be called explicitly.

[Back](#supported-namelists)

### write

With this name list the field or particle distributions are dumped. The placeholder character `@` can be used to refer to the rootname of the simulation run, e.g. `field = @.final`

- `field` (*string,\<empty>*): if a filename is defined, Genesis writes out the field distribution of all harmonics. The harmonics are indicated by the suffix ’.hxxx.’ where xxx is the harmonic number. The filename gets the extension.fld.h5 automatically
- `beam` (*string, \<empty>*): if a filename is defined, Genesis writes out the particle distribution. The filename gets the `extension.par.h5` automatically
- `stride` (*integer,1*): For values larger than 1 the amount of particles written to the file is reduced by only writing each *stride*th particle to the dump file.

[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

### track

This namelist initiate the actually tracking through the undulator and then writing out the results. Normally all parameter should be defined before or defined in the lattice but the namelist allows some ’last minute’ change of the behavior of the code

- `zstop` (*double, 1e9*): If `zstop` is shorter than the lattice length the tracking stops at the specified position.
- `output_step` (*int, 1*): Defines the number of integration steps before the particle and field distribution is analyzed for output.
- `field_dump_step` (*int, 0*): Defines the number of integration steps before a field dump is written. Be careful because for time-dependent simulation it can generate many large output files.
- `beam_dump_step` (*int, 0*): Defines the number of integration steps before a particle dump is written. Be careful because for time-dependent simulation it can generate many large output files.
-  `sort_step` (*int,0*): Defines the number of steps of integration before the particle distribution is sorted. Works only for one-4-one simulations.
- `s0` (*double, <taken from TIME module>*): Option to override the default time window start from the TIME module.
- `slen` (*double, <taken from TIME module>*): Option to override the default time window length from the TIME module.
- `field_dump_at_undexit` (*bool, false*): Field dumps at the exit of the undulator (one dump for each undulator in the expanded lattice).
- `bunchharm` (*int, 1*): Bunching harmonic output setting. Must be >= 1.


[Back](#supported-namelists)

<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

## Internal namelists

The following namelists were not previously documented and have been extracted
from the source code.  These may not have the same level of support when
compared to the above section.

### alter_field

Field manipulator (TODO).

Note that the namelist `field_manipulator` is deprecated and will be removed in
the future. Use `alter_field` instead.

- `harm` (*int, 1*): harmonic
- `scale_power` (*double, 1.0*): power scaling factor
- `spp_l` (*double, 0.0*): TODO
- `spp_nsect` (*int, 0*): TODO
- `spp_phi0` (*double, 0.0*): TODO

### profile_file_multi

Generates profile objects `<label_prefix>.gamma`, `<label_prefix>.delgam`,
`<label_prefix>.current`, etc., each one corresponding to one `&profile_file`.

- `file` (*string, <empty>*): HDF5 filename.
- `label_prefix` (*string, \<empty>*): prefix for each object.
- `xdata` (*string, \<empty>*): Points to a dataset in an HDF5 file to define the `s`-position for the look-up table. The format is `filename/group1/.../groupn/datasetname`, where the naming of groups is not required if the dataset is at root level of the HDF file
- `ydata` (*string, \<empty>*): Same as y data but for the function values of the look-up table.
- `isTime` (*bool, false*): If true the `s`-position is a time variable and therefore multiplied with the speed of light `c` to get the position in meters.
- `reverse`(*bool, false*): if true the order in the look-up table is reverse. This is sometimes needed because time and spatial coordinates differ sometimes by a minus sign.

### sequence_list

A sequence of values given as a string.

- `label` (*string, \<empty>*): label for the sequence.
- `val` (*double, [<empty>]*): list of values.
- `default` (*double, 0*): default value to use for out-of-bound indices.

### sequence_filelist

A sequence list with data in a file.

- `label` (*string, \<empty>*): label for the sequence.
- `file` (*string, \<empty>*): filename to load the sequence from with one line per value.

### add_plugin_fielddiag

Optionally-compiled plugin for field diagnostics.

### add_plugin_beamdiag

Optionally-compiled plugin for beam diagnostics.
