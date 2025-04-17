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
set(_GENESIS_VERSION_REV "8")
set(_GENESIS_VERSION_BETA "false")
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
configure_file(${CMAKE_CURRENT_LIST_DIR}/include/version.h.in ${TARGET_DIR}/include/version.h @ONLY)
