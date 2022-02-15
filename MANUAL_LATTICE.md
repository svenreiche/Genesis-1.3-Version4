# GENESIS 1.3 - Manual
  
** Lattice Elements **
  - Undulator
  - Drift
  - Quadrupole.
  - Corrector
  - Chicane.
  - Phaseshifter
  - Marker
  - Line


## The Lattice File

In comparison to previous versions of Genesis the definition of the undulator lattice is now
completely defined in a lattice file. The main input file just refers to it and might change
some parameters, such as adding undulator errors etc. Another major change is that the
rigid scheme of fixed integration steps is now broken up. The lattice does not need to be
aligned to a single integration step. Instead the lattice can be closer now to the real lattice.

The format of the lattice file is similar to other codes such as ELEGANT or MADX, which
makes it somehow easier to convert between the formats. In general the syntax is the
following:

label: element type={parameter=value[,···]};

Following beamline elements are currently supported: undulator, quadrupole, drift,
corrector, chicane, phaseshifterandmarker. For the parsing of the elements Genesis
only consider the first 4 letters. Thereforeundulatorandunduas element name are both
valid. This applies only for the elements. Other tags, such as labels and parameter have to
be exact. Upper or lower cases are ignored because all letters are converted to lower case
prior to parsing.

Labels are used to identify elements and are referred to in the line element. More information
is given at the end of this section.

## Undulator

- aw: The dimensionless rms undulator parameter. For planar undulator this value is
    smaller by the factor 1/

### √

```
2 than its K-value, while for helical undulator rms and peak
values are identical. Default value is zero.
```
- lambdau: Undulator period length in meter. Default is 0 m.
- nwig: Number of periods
- helical: Boolean flag whether the undulator is planar or helical. Default value is
    false, indicating a planar undulator. Note that setting it to true, does not change the
    roll-off parameters. To be consistent they have to be set directly.
- kx: Roll-off parameter of the quadratic term of the undulator field in x. There are
    normalized with respect toku^2. Default is 0
- ky: Roll-off parameter of the quadratic term of the undulator field in y. Default is 1

### 25


26 GENESIS 1.3 – Manual

- ax: Offset of the undulator module in x in meter. Default is 0.
- ay: Offset of the undulator module in y in meter. Default is 0.
- gradx: Relative transverse gradient of undulator field in x≡(1/aw)∂aw/∂x. Default
    is 0.
- grady: Relative transverse gradient of undulator field in y≡(1/aw)∂aw/∂y. Default
    is 0.

## Drift

- l: Length of the drift in meter. Default value is 0 m.

## Quadrupole.

- l: Length of the quadrupole in meter. Default value is 0 m.
- k1: Normalized focusing strength in 1/m^2. Default value is 0
- dx: Offset in x in meter. Default is zero.
- dy: Offset in y in meter. Default is zero.

## Corrector

- l: Length of the corrector in meter. Default value is 0 m.
- cx: Kick angle in x in units ofγβx. Default value is 0.
- cy: Kick angle in y in units ofγβy. Default value is 0.

## Chicane.

- l: Length of the chicane, which consists out of 4 dipoles without focusing. The first
    and last are placed at the beginning and end of the reserved space. The inner ones
    are defined by the drift length in between. Any remaining distance, namely the length
    subtracted by 4 times the dipole length and twice the drift length are placed between
    the second and third dipole. Default value is 0 m.
- lb: Length of an individual dipole in meter. Default value is 0 m


The Lattice File 27

- ld: Drift between the outer and inner dipoles, projected onto the undulator axis. The
    actual path length is longer by the factor 1/cosθ, whereθis the bending angle of an
    individual dipole.
- delay: Path length difference between the straight path and the actual trajectory in
    meters. From this value the bending angle is calculated internally by Genesis. Default
    delay is 0 m.

## Phaseshifter

- l: Length of the phase shifter in meter. Default value is 0 m.
- phi: Change in the ponderomotive phase of the electrons in units of rad. Default value
    is 0. Note that Genesis is doing an autophasing, so that the electrons at reference
    energy are not changing in ponderomotive phase in drifts.

## Marker

- dumpfield: A non-zero value enforces the dump of the field distribution of this zero
    length element.
- dumpbeam: A non-zero value enforces the dump of the particle distribution.
- sort: A non-zero value enforces the sorting of particles, if one-for-one simulations are
    enabled.
- stop: A non-zero value stops the execution of the tracking module. Note that the
    output file still contains the full length with zeros as output for those integration steps
    which are no further calculated.

## Line

The line uses the labels to define the order of the elements in the beam line. Elements can
be referred to multiple times. Each time the line element creates an individual copy so that
they are not treated as a single entity. That allows that in a later stage the elements can be
changed differently, e.g. when errors in the alignment are generated. The line element can
also contain other line elements. Genesis is unrolling the nested references, however stops
after an iteration depth of 10. This catches the error that one line element uses another and
vice versa.

Elements can be duplicated by a preceding multiplication. The following is an example


28 GENESIS 1.3 – Manual

FODO: Line = {F,DRI,D,DRI};
LAT: Line = {6*FODO};

There is one difference between Genesis and Elegant/Madx. It allows a superposition of
secondary elements (everything except of undulator modules) so that correction and focusing
can be superimposed onto the undulator field. An example is

UND: Undulator = {aw=1, lambdau=0.02, nwig=100};
FODO: Line = {UND,F@0.2,D@1.2};

The postfix @ with the position in meter places the element directly there. Note that it also
resets the counter. For a 10 cm long quadrupole the next position afterF@0.2would be 0.3
m. To avoid overlapping with succeeding elements it is better to define the pointer with a
marker element, e.g. the upper example can be improved with:

UND: Undulator = {aw=1, lambdau=0.02, nwig=100};
FODO: Line = {UND,F@0.2,D@1.2,M@2.0};
LAT: Line = {6*FODO};

If the absolute position places the element outside of an existing element than the missing
space is padded with a drift section.




