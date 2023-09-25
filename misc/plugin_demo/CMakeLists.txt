# BEGIN: Copied from top-level CMakeLists.txt

cmake_minimum_required(VERSION 3.12)
# project(Genesis_1_3_Version4 CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENTIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)
# set(CMAKE_INSTALL_RPATH "/usr/local/lib")
set(CMAKE_CXX_FLAGS '-lstdc++')

# developer options
set(CMAKE_VERBOSE_MAKEFILE ON)
add_compile_options(-Wall -Wextra)

# OpenMPI: preprocessor macro to disable C++ interface (then we don't have to link libmpi_cxx)
add_compile_definitions(OMPI_SKIP_MPICXX)

# END

project(demo VERSION 1.0.0 DESCRIPTION "diagdemo")

add_library(${PROJECT_NAME} SHARED
	DiagFieldPowerDemo.cc
	DiagUtil.cc
)

