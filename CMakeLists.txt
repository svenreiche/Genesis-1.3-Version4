cmake_minimum_required(VERSION 3.10)
project(Genesis_1_3_Version4 CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENTIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



# container for needed librarires
list(APPEND _libraries)


# check for existing compiler
#find_program(H5PCC h5pcc)
#if (H5PCCccc)
    # set compilers
#    Message(STATUS "hdf5 + mpi wrapper found - No additional libraries needed")

 #   list(APPEND _libraries MPI_CXX)
#    SET(CMAKE_C_COMPILER h5pcc)
#    SET(CMAKE_CXX_COMPILER h5pcc)

#else()
#    #SET(CMAKE_C_COMPILER mpicc)
#    #SET(CMAKE_CXX_COMPILER mpicxx)
#    # check for MPI
#    find_package(MPI REQUIRED)
#    Message(STATUS "MPI Found")
#    # and check for HDF5
#    SET(HDF5_PREFER_PARALLEL TRUE)
#    find_package(HDF5 REQUIRED COMPONENTS C HL)
#    Message(STATUS "HDF5 Found")

#    list(APPEND _libraries MPI::MPI_C HDF5::HDF5)
#    # add the include paths for both libraries
#    include_directories(${MPI_INCLUDE_PATH})
#    include_directories(${HDF5_C_INCLUDE_DIRS})
#    if (HDF5_IS_PARALLEL)   # this should be required - needs a better check
#        Message(STATUS "Parallel HDF5 supported")
#    else()
#        Message(STATUS "No parallel HDF5 Support")
#    endif()
#endif()

# include fftw as optional library
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW QUIET fftw3 IMPORTED_TARGET)
if (TARGET PkgConfig::FFTW)
    message(STATUS "FFTW found")
    list(APPEND _libraries PkgConfig::FFTW)
    include_directories(PkgConfig::FFTW)
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



# include directories
include_directories(include)


# some status print out before compilation
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
MESSAGE( STATUS "CMAKE_C_FLAGS: " ${CMAKE_C_FLAGS} )
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )
Message(STATUS "Used libraries: ${_libraries}")

# objects for the main part - maybe needs a better name
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


# add optional macro definition

if (TARGET PkgConfig::FFTW)
    message(STATUS "Defining the macro FFTW in the source code")
    target_compile_definitions(genesis13 PRIVATE FFTW=1)
endif()


add_dependencies(genesis13 get_git_hash)

add_executable(genesis4 src/Main/mainwrap.cpp)

target_link_libraries(genesis4 genesis13 ${_libraries})

