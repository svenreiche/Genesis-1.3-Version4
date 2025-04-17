# GENESIS 1.3 - Lattice File

## The Lattice File

In comparison to previous versions of Genesis the definition of the undulator lattice is now completely defined in a lattice file. The main input file just refers to it and might change some parameters, such as adding undulator errors etc. 
Another major change is that the rigid scheme of fixed integration steps is now broken up. The lattice does not need to be aligned to a single integration step. Instead the lattice can be closer now to the real lattice.

The format of the lattice file is similar to other codes such as ELEGANT or MADX, which makes it somehow easier to convert between the formats. In general the syntax is the following:
```
label: element type={parameter=value[,···]};
```
Following beamline elements are currently supported: undulator, quadrupole, drift, corrector, chicane, phaseshifterandmarker. For the parsing of the elements Genesis only consider the first 4 letters. Therefore, `undulator` and `unduas` element name are both valid. This applies only for the elements. Other tags, such as labels and parameter have to be exact. Upper or lower cases are ignored because all letters are converted to lower case prior to parsing.


## Supported Lattice Elements
  - [undulator](#undulator)
  - [drift](#drift)
  - [quadrupole](#quadrupole)
  - [corrector](#corrector)
  - [chicane](#chicane)
  - [phaseshifter](#phaseshifter)
  - [marker](#marker)
  - [sequence](#sequence)
  - [line](#line)
  

Labels are used to identify elements and are referred to in the line element. More information is given at the end of this section.

### undulator

- `aw` (*double, 0*): The dimensionless rms undulator parameter. For planar undulator this value is smaller by a factor $1 / \sqrt{2}$ than its K-value, while for helical undulator rms and peak values are identical.
- `lambdau` (*double, 0, [m]*): Undulator period length in meter. Default is 0 m.
- `nwig` (*int, 0*): Number of periods.
- `helical` (*bool, false*): Boolean flag whether the undulator is planar or helical. A planar undulator has helical=`false`. Note that setting it to `true`, does not change the roll-off parameters for focusing. To be consistent they have to be set directly.
- `kx` (*double, 0*): Roll-off parameter of the quadratic term of the undulator field in x. It is normalized with respect to $k_u^2$. 
- `ky` (*double, 1*): Roll-off parameter of the quadratic term of the undulator field in y. 
- `ax` (*double, 0, [m]*): Offset of the undulator module in $x$ in meter. 
- `ay` (*double, 0, [m]*): Offset of the undulator module in $y$ in meter. 
- `gradx` (*double, 0*): Relative transverse gradient of undulator field in $x$ $\equiv (1/a_w) \partial a_w/\partial x$.
- `grady` (*double, 0*): Relative transverse gradient of undulator field in $y$ $\equiv (1/a_w) \partial a_w/\partial y$.

[Back](#supported-lattice-elements)

### drift

- `l` (*double, 0, [m]*): Length of the drift in meter.


[Back](#supported-lattice-elements)

### quadrupole

- `l` (*double, 0, [m]*): Length of the quadrupole in meter.
- `k1` (*double, 0, [1/$m^2$]*): Normalized focusing strength in 1/m^2.
- `dx` (*double, 0, [m]*): Offset in $x$ in meter.
- `dy` (*double, 0, [m]*): Offset in $y$ in meter.

[Back](#supported-lattice-elements)

### corrector

- `l` (*double, 0, [m]*): Length of the corrector in meter.
- `cx` (*double, 0, [rad]*): Kick angle in $x$ in units of radians.
- `cy` (*double, 0, [rad]*): Kick angle in $y$ in units of radians.

[Back](#supported-lattice-elements)

### chicane

- `l` (*double, 0, [m]*): Length of the chicane, which consists out of 4 dipoles without focusing. The first and last are placed at the beginning and end of the reserved space. The inner ones are defined by the drift length in between. Any remaining distance, namely the length subtracted by 4 times the dipole length and twice the drift length are placed between the second and third dipole.
- `lb` (*double, 0, [m]*): Length of an individual dipole in meter. 
- `ld` (*double, 0, [m]*): Drift between the outer and inner dipoles, projected onto the undulator axis. The actual path length is longer by the factor $1/\cos\theta$, where $\theta$ is the bending angle of an individual dipole.
- `delay` (*double, 0, [m]*): Path length difference between the straight path and the actual trajectory in meters. Genesis 1.3 calculates the bending angle internally starting from this value. $R_{56} = 2$`delay`.

[Back](#supported-lattice-elements)

### phaseshifter

- `l` (*double, 0, [m]*): Length of the phase shifter in meter.
- `phi` (*double, 0, [rad]*): Change in the ponderomotive phase of the electrons in units of rad. Note that Genesis 1.3 is doing an autophasing, so that the electrons at reference energy are not changing in ponderomotive phase in drifts.

[Back](#supported-lattice-elements)

### marker

- `dumpfield` (*int, 0/1*): A non-zero value enforces the dump of the field distribution of this zero length element.
- `dumpbeam` (*int, 0/1*): A non-zero value enforces the dump of the particle distribution.
- `sort` (*int, 0/1*): A non-zero value enforces the sorting of particles, if one-for-one simulations are enabled.
- `stop` (*int, 0/1*): A non-zero value stops the execution of the tracking module. Note that the output file still contains the full length with zeros as output for those integration steps which are no further calculated.

[Back](#supported-lattice-elements)

### sequence
This element provides the functionality of sequences in the same way as in the main input deck.
Note that it has nothing to do with beamline sequences used in tracking programs, e.g. MadX. 
See the manual for the main input deck for the explicit argument requirements for each sequence.
The label of the sequence is used as the reference to the sequence.

- `type` (*string, empty*) - String defining the type of sequence. It must match of the possibilities given in the input deck, e.g. `const` for a constant sequence
- `...` Same arguments as in the namelists, except that the field `label` is automatically defined by the label field in the lattice file.

An example for a sequence definition is:
```
VAL: SEQUENCE = {type = const, c0 = 3.5};
UND: UNDULATOR = {lambdau=0.015000, nwig=266, aw=@val, helical= True};
```
Note that occurences of sequences are parsed first before evaluating individual elements. Thus the sequence definition can occur at any location in the lattice file.

### line

The line uses the labels to define the order of the elements in the beam line. Elements can be referred to multiple times. Each time the line element creates an individual copy so that they are not treated as a single entity. That allows that in a later stage the elements can be changed differently, e.g. when errors in the alignment are generated. The line element can also contain other line elements. Genesis 1.3 is unrolling the nested references, however stops after an iteration depth of 10. This catches the error that one line element uses another and vice versa.

Elements can be duplicated by a preceding multiplication. The following is an example
```
FODO: Line = {F,DRI,D,DRI};
LAT: Line = {6*FODO};
```
There is one difference between Genesis and Elegant/Madx. It allows a superposition of secondary elements (everything except of undulator modules) so that correction and focusing can be superimposed onto the undulator field. An example is
```
UND: Undulator = {aw=1, lambdau=0.02, nwig=100};
FODO: Line = {UND,F@0.2,D@1.2};
```
The postfix \@ with the position in meter places the element directly there. Note that it also resets the counter. For a 10 cm long quadrupole the next position after ```F@0.2``` would be 0.3 m. To avoid overlapping with succeeding elements it is better to define the pointer with a marker element, e.g. the upper example can be improved with:
```
UND: Undulator = {aw=1, lambdau=0.02, nwig=100};
FODO: Line = {UND,F@0.2,D@1.2,M@2.0};
LAT: Line = {6*FODO};
```
If the absolute position places the element outside of an existing element than the missing space is padded with a drift section.

[Back](#supported-lattice-elements)


<div style="page-break-after: always; visibility: hidden"> \pagebreak </div>

