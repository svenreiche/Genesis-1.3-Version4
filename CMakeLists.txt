cmake_minimum_required(VERSION 3.16)
project(Genesis_1_3_Version4 CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENTIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# check for compiler wrapper
find_program(H5PCC h5pcc)
if (H5PCC)
    Message(STATUS "Wrapper H5PCC found")
    SET(CMAKE_C_COMPILER h5pcc)
    SET(CMAKE_CXX_COMPILER h5pcc)
else()
    Message(STATUS "Wrapper H5PCC not found")
    Message(STATUS "Add and link explicitly to HDF5")
    SET(CMAKE_C_COMPILER mpicc)
    SET(CMAKE_CXX_COMPILER mpicxx)
endif()

# include mpi

find_package(MPI REQUIRED)
if (MPI_FOUND)
    Message(STATUS "MPI Found")
else()
    Message(STATUS "MPI not Found")
endif()

#include hdf5
find_package(HDF5 REQUIRED C)
if (HDF5_FOUND)
    Message(STATUS "HDF5 Found")
else()
    Message(STATUS "HDF5 not Found")
endif()
if (HDF5_IS_PARALLEL)
    Message(STATUS "Parallel HDF5 supported")
else()
    Message(STATUS "No parallel HDF5 Support")
endif()
SET(HDF5_PREFER_PARALLEL TRUE)


# include fftw as optional library
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW QUIET fftw3 IMPORTED_TARGET)
if (TARGET PkgConfig::FFTW)
    message(STATUS "FFTW found")
else()
    message(STATUS "FFTW not found")
    message(STATUS "Excluding FFTW specific source code")
endif()




# generation of the version include file
SET(GIT_HASH "unknown")
find_package(Git QUIET)
if (GIT_FOUND)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H
            OUTPUT_VARIABLE GIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            WORKING_DIRECTORY
               ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
message(STATUS "Git hash is ${GIT_HASH}")
set(_GENESIS_VERSION_MAJOR "4")
set(_GENESIS_VERSION_MINOR "6")
set(_GENESIS_VERSION_REV "1")
set(_GENESIS_VERSION_BETA "true")
execute_process(
        COMMAND
            whoami
        TIMEOUT
            1
        OUTPUT_VARIABLE
            _user_name
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(TIMESTAMP _config_time "%Y-%m-%d %H:%M:%S [UTC]" UTC)
configure_file(include/version.h.in include/version.h @ONLY)



# include directories
include_directories(PkgConfig::FFTW)
include_directories(${MPI_INCLUDE_PATH})
Message(STATUS "Path for MPI include files: ${MPI_INCLUDE_PATH}")
include_directories(${HDF5_C_INCLUDE_DIRS})
Message(STATUS "Path for HDF5 include files: ${HDF5_C_INCLUDE_DIR}")
include_directories(include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

# list of source files

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

if (TARGET PkgConfig::FFTW)
  target_compile_definitions(genesis13 PRIVATE FFTW=1)
endif()


add_executable(genesis4 src/Main/mainwrap.cpp)
#target_link_libraries(genesis4 genesis13 MPI::MPI_C HDF5::HDF5)
if (TARGET PkgConfig::FFTW)
    target_link_libraries(genesis4 genesis13 MPI::MPI_C PkgConfig::FFTW)
#    target_link_libraries(genesis4 genesis13 MPI::MPI_C HDF5::HDF5 PkgConfig::FFTW)
else()
    target_link_libraries(genesis4 genesis13 MPI::MPI_C)
#    target_link_libraries(genesis4 genesis13 MPI::MPI_C HDF5::HDF5)
endif()
