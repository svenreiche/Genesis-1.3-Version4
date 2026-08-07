# Every subdirectory containing a run.sh becomes one CTest test. The cases are
# driven by shell scripts rather than by CMake so that they can also be run from
# a source tree without a build directory, or against an already installed
# executable, which is what a packaging recipe needs.

find_program(BASH_EXECUTABLE bash)
find_package(Python3 COMPONENTS Interpreter)

set(GENESIS_TEST_RANKS "1 2 3"
    CACHE STRING "MPI rank counts that the test suite compares against each other")

if(NOT BASH_EXECUTABLE)
    message(STATUS "bash not found, the test suite will not be registered")
    return()
endif()

file(GLOB test_runners ${CMAKE_CURRENT_SOURCE_DIR}/*/run.sh)

foreach(runner ${test_runners})
    get_filename_component(case_dir ${runner} DIRECTORY)
    get_filename_component(case_name ${case_dir} NAME)

    add_test(NAME ${case_name}
             COMMAND ${BASH_EXECUTABLE} ${runner}
             WORKING_DIRECTORY ${case_dir})

    # A case exits with 77 when it cannot run for want of MPI or of a Python
    # interpreter with h5py. That is reported as a skip rather than a failure.
    set_tests_properties(${case_name} PROPERTIES
        SKIP_RETURN_CODE 77
        ENVIRONMENT
            "GENESIS4=$<TARGET_FILE:genesis4>;MPIEXEC=${MPIEXEC_EXECUTABLE};PYTHON=${Python3_EXECUTABLE};RANKS=${GENESIS_TEST_RANKS};WORKDIR=${CMAKE_CURRENT_BINARY_DIR}/work")
endforeach()
