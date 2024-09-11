## Example 3 : Time-dependent Simulation

*All files for running the example are found in the subdirectory examples/Example3-TimeDependent of the source code distribution*

Time-Dependent differs from steady-state simulation by tracking multiple slices of the electron bunch through the undulator. This
allows to resolve a band of frequencies, depending on the sample rate of the electron slices. 
For convention a set of slices is defined in the bunch-frame, which moves with the speed of the electrons. Therefore mostly the electrons are
staying at the same location while the radiation slips forward, namely by one wavelength per undulator period. In Genesis the coordinate 
in the bunch frame is given by **s**. This shouldn't be mixed with **z** which is the position in the undulator where the bunch frame is evaluated. 
In the underlying model the variable **z** acts as the independent variable of a proposed Hamiltonian 

This example is based on the steady-state input deck of the first example. To run time-dependent simulation the lattice file does not need to be changed. 

#### Input File

The definition of the beam-frame is done with the &time namelist.

```asm
&time
slen = 20e-6
sample = 3
&end
```

Here the total length is 20 microns and with a sample rate of 3 times the reference wavelength which is defined as **lambda0** in the main setup namelist.
The sample rate cannot be smaller than 1 but can be increased as long as the sample rate is still smaller than the cooperation length. The choice of the sample rate also defines
the spectral bandwidth according to Nyquist theorem. For the maximum sample rate the frequency is +/-50% around the central frequency. With larger sample rate the bandwidth is restricted
to +/- 50/**sample** percent.  Due to the underlying solver the FEL should be well within the given bandwidth of the simulation. Where this is valid or not, will not be checked by Genesis (Example would be
simulation with a very large energy chirp in the electron beam or seed signal well different from the central frequency). On the other hand the resonant frequency of the FEL does not need to be
exact the same as the central frequency.

There is also one possible pitfall that the desired integration stepsize is not in conflict with the sample rate. The given example requests an integration step size of 45 mm, which corresponds to
three undulator periods. This means that after one single integration step the field is pushed three wavelength. If the sample rate is smaller, e.g. 1 than after one integration step the field is jumping over electron slices without interaction.
As an example for sample = 1 the field in slice is copied from slice i to slice i+3.  There is no interaction with slice i+1 and i+2. As a consequence the simulation acts as three independent simulation 
with a sample rated of 3 each.  For this reason the example uses a sample rate of 3.

To illustrate time-varying input parameter, we set up the simulation for a Gaussian current profile and a linear chirp.
For that some profiles are defined with 
```asm
&profile_gauss
label=current
c0 = 2500
s0 = 10e-6
sig = 6e-6
&end

&profile_polynom
label=energy
c0=11347.
c1=1e7.
&end
```
The upper is the current profile with an rms length of 6 microns. Since the beam frame is defined for 20 microns, starting at s=0 micron. The peak of the Gaussian is placed at s=10 microns with the **s0** parameter.

The lower namelist is the energy chirp of the electron bunch. That namelist supports higher polynomial order but the corresponding coefficent are defaulting to zero.
Here only the constant and linear term is defined.

To assign these profiles to the beam parameters they are referred to by their name:

```asm
&beam
current=@current
gamma = @energy
...
&end
```

Time-simulation using much more resources than steady-state simulation and sufficient memory much be provided for storing the entire beam and field in memory.
Also the output files can be very large. Therefore it is better to reduced them.
One method is to only generate output after a given count of integration step. This is done in the track namelist with:

```asm
&track
output_step=3
&end
```
In this case after each 3rd integration step there will be an output

The other methode is to exclude certain fields in the output in combination for a collective output of all slices, such as the mean energy in the electron beam.
This is controlled in the setup namelist:

```asm
&setup
...
field_global_stat = true
beam_global_stat = true
exclude_spatial_output = true
exclude_fft_output = true
&end
```

**exclude_spatial_output** supresses any spatial information (cente rposition and rms sizes) while **exclude__fft_output** affects the field divergent. It has also the 
benefit that the simulation runs a little bit faster since the 2D FFT of the field wave-fronts can be CPU intensive

Global stats parameters can be fined in the group **Global** in the correpsonding field or beam group.


The last step to run a SASE simulation to set the radiaiton power to zero and set the shotnoise parameter in the setup namelist to one.


## Output
Since the simulation is tracking thousands of slice instead of a single slice for steady-state simulation, Genesis should be started with the mpirun command.

A typical launching command is ```mpirun -np 100 genesis4 Example3.in``` and it will take some time.
the output is similar to the steady-state simulation except some informations are added on the number of slices and the resulting beam-frame size.
Since the MPI execution of Genesis requires a symmetric number of slices per node it might not be able to resolve the requested size of the beam-frame.
To solve this Genesis will add that many slices to the last node that it holds the same amount as the others. This causes an extension of the time-window.
This depends on the number of modes chosen.

In addition the shotnoise depends on the number of modes, where 1000 slices are resolved by 20 or 50 nodes. This will be fixed in a future release for a better reporducibility of the simulation results.


The python script ```Example3.py``` will parse the output file and can be used as a template on how read the output file.
To run the script python should have the packages **matplotlib** and **h5py** installed.
The generated figure should be identical to the following plots.

#### Time-dependent beam parameters
![Plot1](Plots/Figure_1.png)

To verify the requested shape the beam current and the energy at the simulation start is plotted.

#### Pulse Energy and Far Field Intensity

![Plots2](Plots/Figure_2.png)

The Global group in the Field group holds the evolution of the pulse energy and the on-axis far field intensity.
Since Genesis allows for a lot of higher transverse modes to be emitted they are included in the calculation of the pulse energy for
SASE simulation. This generates the bump in the semi-log plot. Often it is hard to see if there is an exponential gain.
For that it is convenient to plot the on-axis far-field intensity. Since it has no contribution from higher transverse modes it is a "cleaner" signal for the FEL process.

#### Radiation Profile

![Plots3](Plots/Figure_3.png)

This is also visible in the profiles of the radiation power and farfield intensity. In the start-up the power follows quite closely
the current profile since it acts more like spontaneous radiation which is proportional to the current. The plot is taken after the first undulator module.

![Plots4](Plots/Figure_4.png)

Further downstream the power profile has more fluctuation similar to the farfield intensity. 

![Plots5](Plots/Figure_5.png)

A zoom in the core part of the beam-frame shows that close at saturation the profile of power and far-field intensity has roughly the same shape.
The slight shift in the power is still causes from some parasitic higher modes and indicates that the transverse coherence is not fully at 100%.

![Plots6](Plots/Figure_6.png)

A more compact plot showing the evolution is when for each step in z the profile is normalized to its mean values. Otherwise one could not see the start-up regime with the exponential growth of 
the radiation power. In this 2D plot the spiky structure becomes apparent after about half the undulator length. These spikes more forward in the beam-frame with a group velocity which is less than the spead of light.
At saturation this changes (the slope of the spikes is slightly changed). Also after saturation the power in the head and tail can reach saturation and the pulse gets longer.

![Plots7](Plots/Figure_7.png)

From the farfield intensity and phase the spectrum can be calculated. In this case the spectrum is wider than expected since the beam has an intrinsic energy chirp and thus boradening the emitted frequency spectrum
