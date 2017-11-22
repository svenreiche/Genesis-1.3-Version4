#
# GENESIS makefile
#
# libraries
#
#LIB= -lm -lmpi_cxx -lz -lgfortran -lmpe -L/gpfs/home/reiche/Apps/MPE/lib -llmpe 
#LIB= -lm -lmpi_cxx -lz -lgfortran  -lmpe -L/gpfs/home/reiche/Apps/MPE/lib
#LIB= -lm -lmpi_cxx -lz 
LIB= -lm -lstdc++ -lmpi_cxx
#
#INCLUDE=-I./include -I/gpfs/home/reiche/Apps/MPE/include
INCLUDE=-I./include
#
# cpp-macros
#
DMACRO = -DMPEEE
#
#  compilers
#
VPATH = src/Core src/IO src/Lattice src/Util src/Main src/Loading
CCOMPILER = h5pcc
FCOMPILER = gfortran
#
#  executable name 
#
EXECUTABLE = gencore
#
# targets
#
OBJECTS = Sorting.o BesselJ.o Inverfc.o Hammerslay.o RandomU.o GaussHermite.o StringProcessing.o Track.o Setup.o AlterSetup.o Time.o Parser.o Dump.o SponRad.o EField.o LoadBeam.o ImportBeam.o LoadField.o ImportField.o Profile.o ShotNoise.o QuietLoading.o Optics.o Lattice.o LatticeElements.o LatticeParser.o AlterLattice.o Gencore.o TrackBeam.o Control.o Field.o FieldSolver.o EFieldSolver.o Incoherent.o Beam.o BeamSolver.o Undulator.o HDF5base.o readBeamHDF5.o writeBeamHDF5.o readFieldHDF5.o writeFieldHDF5.o SDDSBeam.o Output.o  MPEProfiling.o main.o

genesis:	$(OBJECTS)
	$(CCOMPILER)  -o $(EXECUTABLE) $(OBJECTS) $(LIB)



.cpp.o:
	$(CCOMPILER) -O -c $(DMACRO) $(INCLUDE) $<

.f.o:
	$(FCOMPILER) -O  -c $<

clean:
	rm -f *~
	rm -f *.o
	rm $(EXECUTABLE)

install:
	cp ./$(EXECUTABLE) ~/bin/genesis4

beta:
	cp ./$(EXECUTABLE) ~/bin/$(EXECUTABLE).beta











