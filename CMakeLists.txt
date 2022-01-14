cmake_minimum_required(VERSION 3.10)
project(Genesis_1_3_Version4 CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENTIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#SET(CMAKE_C_COMPILER h5pcc)
#SET(CMAKE_CXX_COMPILER h5pcc)



# include mpi
#find_package(MPI REQUIRED)
#if (MPI_FOUND)
#    Message(STATUS "MPI Found")
#else()
#    Message(STATUS "MPI not Found")
#endif()

#include hdf5
#find_package(HDF5 REQUIRED COMPONENTS C HL)
#if (HDF5_FOUND)
#    Message(STATUS "HDF5 Found")
#else()
#    Message(STATUS "HDF5 not Found")
#endif()
#if (HDF5_IS_PARALLEL)
#    Message(STATUS "Parallel HDF5 supported")
#else()
#    Message(STATUS "No parallel HDF5 Support")
#endif()
#SET(HDF5_PREFER_PARALLEL TRUE)

# include fftw as optional library
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW QUIET fftw3 IMPORTED_TARGET)
if (TARGET PkgConfig::FFTW)
    message(STATUS "FFTW found")
else()
    message(STATUS "FFTW not found")
    message(STATUS "Excluding FFTW specific source code")
endif()


add_custom_command(
    OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/include/version.h
    ALL
    COMMAND
    ${CMAKE_COMMAND} -D TARGET_DIR=${CMAKE_CURRENT_BINARY_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/git-hash.cmake
    WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(
    get_git_hash
    ALL
    DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/include/version.h
   
)



# set compilers
#SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_CXX_COMPILER mpicxx)

# include directories
include_directories(PkgConfig::FFTW)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${HDF5_C_INCLUDE_DIRS})
include_directories(include)



add_library(genesis13 STATIC
        src/Util/BesselJ.cpp
        src/Util/GaussHermite.cpp
        src/Util/Hammerslay.cpp
        src/Util/Inverfc.cpp
        src/Util/RandomU.cpp
        src/Util/Sorting.cpp
        src/Util/StringProcessing.cpp
        src/Lattice/AlterLattice.cpp
        src/Lattice/Lattice.cpp
        src/Lattice/LatticeElements.cpp
        src/Lattice/LatticeParser.cpp
        src/Lattice/Optics.cpp
        src/Core/Beam.cpp
        src/Core/BeamSolver.cpp
        src/Core/Collective.cpp
        src/Core/Control.cpp
        src/Core/EFieldSolver.cpp
        src/Core/Field.cpp
        src/Core/FieldSolver.cpp
        src/Core/Gencore.cpp
        src/Core/Incoherent.cpp
        src/Core/TrackBeam.cpp
        src/Core/Undulator.cpp
        src/Core/Diagnostic.cpp
        src/Core/DiagnosticUser.cpp
        src/IO/HDF5base.cpp
        src/IO/Output.cpp
        src/IO/readBeamHDF5.cpp
        src/IO/readFieldHDF5.cpp
        src/IO/SDDSBeam.cpp
        src/IO/writeBeamHDF5.cpp
        src/IO/writeFieldHDF5.cpp
        src/Loading/ImportBeam.cpp
        src/Loading/ImportField.cpp
        src/Loading/LoadBeam.cpp
        src/Loading/LoadField.cpp
        src/Loading/Profile.cpp
        src/Loading/QuietLoading.cpp
        src/Loading/Series.cpp
        src/Loading/ShotNoise.cpp
        src/Main/AlterSetup.cpp
        src/Main/Dump.cpp
        src/Main/EField.cpp
        src/Main/GenMain.cpp
        src/Main/mainwrap.cpp
        src/Main/Parser.cpp
        src/Main/Setup.cpp
        src/Main/SponRad.cpp
        src/Main/Time.cpp
        src/Main/Track.cpp
        src/Main/Wake.cpp
        )

target_include_directories(genesis13
	PRIVATE
	${CMAKE_CURRENT_BINARY_DIR}/include
)

if (TARGET PkgConfig::FFTW)
  target_compile_definitions(genesis13 PRIVATE FFTW=1)
endif()


add_dependencies(genesis13 get_git_hash)

add_executable(genesis4 src/Main/mainwrap.cpp)
if (TARGET PkgConfig::FFTW)
    target_link_libraries(genesis4 genesis13 PkgConfig::FFTW)
else()
    target_link_libraries(genesis4 genesis13)     
endif()
