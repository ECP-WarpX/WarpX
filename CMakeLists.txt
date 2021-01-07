# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.15.0)
project(WarpX VERSION 0.20.5)

include(${WarpX_SOURCE_DIR}/cmake/WarpXFunctions.cmake)

# In-source tree builds are messy and can screw up the build system.
# Avoid building at least in the same dir as the root dir:
if(CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Building in-source is not supported! "
            "Create a build directory and remove "
            "${CMAKE_SOURCE_DIR}/CMakeCache.txt ${CMAKE_SOURCE_DIR}/CMakeFiles/")
endif()


# CMake policies ##############################################################
#
# CMake 3.18+: CMAKE_CUDA_ARCHITECTURES
# https://cmake.org/cmake/help/latest/policy/CMP0104.html
if(POLICY CMP0104)
    cmake_policy(SET CMP0104 OLD)
endif()


# CCache Support ##############################################################
#
# this is an optional tool that stores compiled object files; allows fast
# re-builds even with "make clean" in between. Mainly used to store AMReX
# objects
set_ccache()


# Output Directories ##########################################################
#
# temporary build directories
set_default_build_dirs()

# install directories
set_default_install_dirs()


# Options and Variants ########################################################
#
option(WarpX_APP           "Build the WarpX executable application"     ON)
option(WarpX_ASCENT        "Ascent in situ diagnostics"                 OFF)
option(WarpX_LIB           "Build WarpX as a shared library"            OFF)
option(WarpX_MPI           "Multi-node support (message-passing)"       ON)
option(WarpX_OPENPMD       "openPMD I/O (HDF5, ADIOS)"                  OFF)
option(WarpX_PSATD         "spectral solver support"                    OFF)
option(WarpX_QED           "PICSAR QED (requires Boost and PICSAR)"     OFF)
# TODO: sensei, legacy hdf5?

set(WarpX_DIMS_VALUES 2 3 RZ)
set(WarpX_DIMS 3 CACHE STRING "Simulation dimensionality (2/3/RZ)")
set_property(CACHE WarpX_DIMS PROPERTY STRINGS ${WarpX_DIMS_VALUES})
if(NOT WarpX_DIMS IN_LIST WarpX_DIMS_VALUES)
    message(FATAL_ERROR "WarpX_DIMS (${WarpX_DIMS}) must be one of ${WarpX_DIMS_VALUES}")
endif()

set(WarpX_PRECISION_VALUES SINGLE DOUBLE)
set(WarpX_PRECISION DOUBLE CACHE STRING "Floating point precision (SINGLE/DOUBLE)")
set_property(CACHE WarpX_PRECISION PROPERTY STRINGS ${WarpX_PRECISION_VALUES})
if(NOT WarpX_PRECISION IN_LIST WarpX_PRECISION_VALUES)
    message(FATAL_ERROR "WarpX_PRECISION (${WarpX_PRECISION}) must be one of ${WarpX_PRECISION_VALUES}")
endif()

set(WarpX_COMPUTE_VALUES NOACC OMP CUDA SYCL HIP)
set(WarpX_COMPUTE OMP CACHE STRING "On-node, accelerated computing backend (NOACC/OMP/CUDA/SYCL/HIP)")
set_property(CACHE WarpX_COMPUTE PROPERTY STRINGS ${WarpX_COMPUTE_VALUES})
if(NOT WarpX_COMPUTE IN_LIST WarpX_COMPUTE_VALUES)
    message(FATAL_ERROR "WarpX_COMPUTE (${WarpX_COMPUTE}) must be one of ${WarpX_COMPUTE_VALUES}")
endif()

option(WarpX_MPI_THREAD_MULTIPLE "MPI thread-multiple support, i.e. for async_io" ON)
mark_as_advanced(WarpX_MPI_THREAD_MULTIPLE)

set(WarpX_PARSER_DEPTH 24 CACHE STRING
    "Maximum parser depth for input file functions")
mark_as_advanced(WarpX_PARSER_DEPTH)

option(WarpX_amrex_internal                    "Download & build AMReX" ON)

# change the default build type to RelWithDebInfo (or Release) instead of Debug
set_default_build_type("RelWithDebInfo")

# this defined the variable BUILD_TESTING which is ON by default
#include(CTest)


# Dependencies ################################################################
#

# AMReX
#   builds AMReX from source (default) or finds an existing install
include(${WarpX_SOURCE_DIR}/cmake/dependencies/AMReX.cmake)
#   suppress warnings in AMReX headers (use -isystem instead of -I)
make_third_party_includes_system(AMReX::amrex AMReX)

# PICSAR
#   builds PICSAR from source
include(${WarpX_SOURCE_DIR}/cmake/dependencies/PICSAR.cmake)

# openPMD
#   builds openPMD-api from source (default) or finds an existing install
include(${WarpX_SOURCE_DIR}/cmake/dependencies/openPMD.cmake)

# PSATD
if(WarpX_PSATD)
    # FFTW (non-GPU) and cuFFT (GPU)
    if(NOT WarpX_COMPUTE STREQUAL CUDA)
        find_package(PkgConfig REQUIRED QUIET)
        if(WarpX_PRECISION STREQUAL "DOUBLE")
            pkg_check_modules(fftw3 REQUIRED IMPORTED_TARGET fftw3)
        else()
            pkg_check_modules(fftw3f REQUIRED IMPORTED_TARGET fftw3f)
        endif()
    endif()
    # BLASPP and LAPACKPP
    if(WarpX_DIMS STREQUAL RZ)
        # FIXME send a couple of upstream PRs to those repos...
        find_package(blaspp REQUIRED)
        find_package(lapackpp REQUIRED)
        find_package(OpenMP REQUIRED)  # pulled by the two above
    endif()
endif()


# Targets #####################################################################
#
if(NOT WarpX_APP AND NOT WarpX_LIB)
    message(FATAL_ERROR "Need to build at least WarpX app or "
                        "library/Python bindings")
endif()

# collect all objects for compilation
add_library(WarpX OBJECT)
set(_ALL_TARGETS WarpX)

# executable application
#   note: we currently avoid a dependency on a core library
#         for simpler usage, but could make this an option
if(WarpX_APP)
    add_executable(app)
    add_executable(WarpX::app ALIAS app)
    target_link_libraries(app PRIVATE WarpX)
    set(_BUILDINFO_SRC app)
    list(APPEND _ALL_TARGETS app)
endif()

# link into a shared library
if(WarpX_LIB)
    add_library(shared MODULE)
    add_library(WarpX::shared ALIAS shared)
    target_link_libraries(shared PUBLIC WarpX)
    set(_BUILDINFO_SRC shared)
    list(APPEND _ALL_TARGETS shared)

    set_target_properties(WarpX shared PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()

# own headers
target_include_directories(WarpX PUBLIC
    $<BUILD_INTERFACE:${WarpX_SOURCE_DIR}/Source>
)

# if we include <AMReX_buildInfo.H> we will need to call:
include(AMReXBuildInfo)
generate_buildinfo(${_BUILDINFO_SRC} "${WarpX_SOURCE_DIR}")
target_link_libraries(WarpX PRIVATE buildInfo::${_BUILDINFO_SRC})
unset(_BUILDINFO_SRC)

# add sources
target_sources(WarpX PRIVATE Source/WarpX.cpp)
if(WarpX_APP)
    target_sources(app PRIVATE Source/main.cpp)
endif()

add_subdirectory(Source/BoundaryConditions)
add_subdirectory(Source/Diagnostics)
add_subdirectory(Source/Evolve)
add_subdirectory(Source/FieldSolver)
add_subdirectory(Source/Filter)
add_subdirectory(Source/Initialization)
add_subdirectory(Source/Laser)
add_subdirectory(Source/Parallelization)
add_subdirectory(Source/Parser)
add_subdirectory(Source/Particles)
add_subdirectory(Source/Python)
add_subdirectory(Source/Utils)

# C++ properties: at least a C++14 capable compiler is needed
foreach(warpx_tgt IN LISTS _ALL_TARGETS)
    target_compile_features(${warpx_tgt} PUBLIC cxx_std_14)
endforeach()
set_target_properties(${_ALL_TARGETS} PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

# link dependencies
target_link_libraries(WarpX PUBLIC WarpX::thirdparty::AMReX)

if(WarpX_PSATD)
    if(WarpX_COMPUTE STREQUAL CUDA)
        # CUDA_ADD_CUFFT_TO_TARGET(WarpX)
        target_link_libraries(WarpX PUBLIC cufft)
    else()
        if(WarpX_PRECISION STREQUAL "DOUBLE")
            target_link_libraries(WarpX PUBLIC PkgConfig::fftw3)
        else()
            target_link_libraries(WarpX PUBLIC PkgConfig::fftw3f)
        endif()
    endif()
    if(WarpX_DIMS STREQUAL RZ)
        # FIXME send a couple of upstream PRs to those repos...
        target_link_libraries(WarpX PUBLIC blaspp)
        target_link_libraries(WarpX PUBLIC lapackpp)
    endif()
endif()

if(WarpX_OPENPMD)
    target_compile_definitions(WarpX PUBLIC WARPX_USE_OPENPMD)
    target_link_libraries(WarpX PUBLIC openPMD::openPMD)
endif()

if(WarpX_QED)
    target_compile_definitions(WarpX PUBLIC WARPX_QED)
    # FIXME WARPX_QED_TABLE_GEN
    target_link_libraries(WarpX PUBLIC pxrmp_qed)
endif()

# AMReX helper function: propagate CUDA specific target & source properties
if(WarpX_COMPUTE STREQUAL CUDA)
    setup_target_for_cuda_compilation(WarpX)
    foreach(warpx_tgt IN LISTS _ALL_TARGETS)
        setup_target_for_cuda_compilation(${warpx_tgt})
    endforeach()
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
        foreach(warpx_tgt IN LISTS _ALL_TARGETS)
            target_compile_features(${warpx_tgt} PUBLIC cuda_std_14)
        endforeach()
        set_target_properties(${_ALL_TARGETS} PROPERTIES
            CUDA_EXTENSIONS OFF
            CUDA_STANDARD_REQUIRED ON
        )
    endif()
endif()

# fancy binary name for build variants
set_warpx_binary_name()


# Defines #####################################################################
#
if(WarpX_DIMS STREQUAL 3)
    target_compile_definitions(WarpX PUBLIC WARPX_DIM_3D)
elseif(WarpX_DIMS STREQUAL 2)
    target_compile_definitions(WarpX PUBLIC WARPX_DIM_XZ)
elseif(WarpX_DIMS STREQUAL RZ)
    target_compile_definitions(WarpX PUBLIC WARPX_DIM_RZ)
endif()

if(WarpX_OPENPMD)
    target_compile_definitions(WarpX PUBLIC WARPX_USE_OPENPMD)
endif()

if(WarpX_QED)
    target_compile_definitions(WarpX PUBLIC WARPX_QED)
    if(WarpX_QED_TABLE_GEN)
        target_compile_definitions(WarpX PUBLIC WarpX_QED_TABLE_GEN)
    endif()
endif()

if(WarpX_PSATD)
    target_compile_definitions(WarpX PUBLIC WARPX_USE_PSATD)
endif()

target_compile_definitions(WarpX PUBLIC
    WARPX_PARSER_DEPTH=${WarpX_PARSER_DEPTH})

# <cmath>: M_PI
# <unistd.h>: parser
if(WIN32)
    target_compile_definitions(WarpX PRIVATE _USE_MATH_DEFINES)
    target_compile_definitions(WarpX PRIVATE YY_NO_UNISTD_H)
endif()


# Warnings ####################################################################
#
set_cxx_warnings()


# Generate Configuration and .pc Files ########################################
#
# these files are used if WarpX is installed and picked up by a downstream
# project (not needed yet)

#include(CMakePackageConfigHelpers)
#write_basic_package_version_file("WarpXConfigVersion.cmake"
#    VERSION ${WarpX_VERSION}
#    COMPATIBILITY SameMajorVersion
#)


# Installs ####################################################################
#
# headers, libraries and executables
set(WarpX_INSTALL_TARGET_NAMES)
if(WarpX_APP)
    list(APPEND WarpX_INSTALL_TARGET_NAMES app)
endif()
if(WarpX_LIB)
    list(APPEND WarpX_INSTALL_TARGET_NAMES shared)
endif()

install(TARGETS ${WarpX_INSTALL_TARGET_NAMES}
    EXPORT WarpXTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# CMake package file for find_package(WarpX::WarpX) in depending projects
#install(EXPORT WarpXTargets
#    FILE WarpXTargets.cmake
#    NAMESPACE WarpX::
#    DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
#)
#install(
#    FILES
#        ${WarpX_BINARY_DIR}/WarpXConfig.cmake
#        ${WarpX_BINARY_DIR}/WarpXConfigVersion.cmake
#    DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
#)


# Tests #######################################################################
#

#if(BUILD_TESTING)
#    enable_testing()
#
#    add_test(...)
#endif()


# Status Summary for Build Options ############################################
#
warpx_print_summary()
