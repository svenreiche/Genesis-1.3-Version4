project(g4_tb)

add_executable(g4_tb g4_tb.cpp)

### copied from top-level CMakeLists.txt file ###
# The -DFFTW=1 is cruical: otherwise the storage layout of class
# DiagnosticFieldBase is different in g4_tb and the code blocks
# compiled for the "GENESIS 1.3" binary (such as
# src/Core/DiagnosticHook.cpp and src/Core/DiagnosticHookLI.cpp)
if (TARGET PkgConfig::FFTW)
    message(STATUS "Defining the macro FFTW in the source code")
    target_compile_definitions(g4_tb PRIVATE FFTW=1)
endif()
if(USE_DPI)
   target_compile_definitions(g4_tb PRIVATE USE_DPI)
endif()
### end of copied definitions ###


target_link_libraries(g4_tb genesis13 ${_libraries})
