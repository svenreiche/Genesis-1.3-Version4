#
# GENESIS makefile
#
# libraries
#
# on some systems -lmpi_cxx is needed to build
LIB= -lm -lstdc++ -lmpi_cxx -lfftw3 
#LIB= -lm -lstdc++ -lfftw3 
#
INCLUDE=-I./include
#
# cpp-macros - to enable the FFTW library
#
DMACRO = -DFFTW
#
#  compilers
#
VPATH = src/Core src/IO src/Lattice src/Util src/Main src/Loading
CCOMPILER = h5pcc
#
#  flags
#
# FLAGS = -O2
FLAGS = -g -O0 -Wall
#
#  executable name
#
EXECUTABLE = gencore
#
# profiling
#
#SCOREP = /cluster/home/sreiche/profiling/scorep/bin/scorep
#SCORDEF = -DSCOREP
#
# targets
#
OBJECTS = Sorting.o BesselJ.o Inverfc.o Hammerslay.o RandomU.o GaussHermite.o StringProcessing.o Track.o Setup.o AlterSetup.o Time.o Wake.o Parser.o Dump.o SponRad.o EField.o LoadBeam.o ImportBeam.o LoadField.o ImportField.o Profile.o Series.o ShotNoise.o QuietLoading.o Optics.o Lattice.o LatticeElements.o LatticeParser.o AlterLattice.o Gencore.o TrackBeam.o Control.o Field.o FieldSolver.o EFieldSolver.o Incoherent.o Collective.o Beam.o BeamDiag.o BeamDiag_Demo.o BeamSolver.o Undulator.o HDF5base.o readBeamHDF5.o writeBeamHDF5.o readFieldHDF5.o writeFieldHDF5.o SDDSBeam.o Output.o GenMain.o VersionInfo.o

.PHONY: genesis genesisexecutable clean install beta

genesis:	$(OBJECTS) build_info.o
	ar -cvq libgenesis13.a $(OBJECTS) build_info.o
	mv libgenesis13.a ./lib
	$(SCOREP) $(CCOMPILER) src/Main/mainwrap.cpp -o $(EXECUTABLE) $(INCLUDE) -lgenesis13 -Llib $(LIB)

genesisexecutable:	$(OBJECTS)
	$(SCOREP) $(CCOMPILER)  -o $(EXECUTABLE) $(OBJECTS) $(LIB)



.cpp.o:
	$(SCOREP) $(CCOMPILER) $(FLAGS) -c $(DMACRO) $(SCORDEF) $(INCLUDE) $<

### rules for build info
build_info.o: build_info.c
	$(SCOREP) $(CCOMPILER) $(FLAGS) $(INCLUDE) -c $<
build_info.c: FORCE
	rm -f build_info.c
	./build_info.sh
# so-called "force target", see GNU make manual "rules without recipes or prerequisites"
FORCE:

### end rules for build info

clean:
	rm -f src/Core/*~
	rm -f src/IO/*~
	rm -f src/Lattice/*~
	rm -f src/Loading/*~
	rm -f src/Main/*~
	rm -f src/Util/*~
	rm -f include/*~
	rm -f *.o
	rm -f lib/*.a
	rm -f build_info.o build_info.c
	rm -f $(EXECUTABLE)

install:
	cp ./$(EXECUTABLE) ~/bin/genesis4

beta:
	cp ./$(EXECUTABLE) ~/bin/$(EXECUTABLE).beta
