# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.14.0)
project(WarpX VERSION 0.20.5)

include(${WarpX_SOURCE_DIR}/cmake/WarpXFunctions.cmake)

# In-source tree builds are messy and can screw up the build system.
# Avoid building at least in the same dir as the root dir:
if(CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Building in-source is not supported! "
            "Create a build directory and remove "
            "${CMAKE_SOURCE_DIR}/CMakeCache.txt ${CMAKE_SOURCE_DIR}/CMakeFiles/")
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
option(WarpX_ASCENT        "Ascent in situ diagnostics"                 OFF)
option(WarpX_MPI           "Multi-node support (message-passing)"       ON)
option(WarpX_OPENPMD       "openPMD I/O (HDF5, ADIOS)"                  OFF)
option(WarpX_PSATD         "spectral solver support"                    OFF)
option(WarpX_QED           "PICSAR QED (requires Boost and PICSAR)"     OFF)
# TODO: python, sensei, legacy hdf5?

set(WarpX_DIMS_VALUES 2 3 RZ)
set(WarpX_DIMS 3 CACHE STRING "Simulation dimensionality (2/3/RZ)")
set_property(CACHE WarpX_DIMS PROPERTY STRINGS ${WarpX_DIMS_VALUES})
if(NOT WarpX_DIMS IN_LIST WarpX_DIMS_VALUES)
    message(FATAL_ERROR "WarpX_DIMS (${WarpX_DIMS}) must be one of ${WarpX_DIMS_VALUES}")
endif()

set(WarpX_PRECISION_VALUES double single)
set(WarpX_PRECISION double CACHE STRING "Floating point precision (single/double)")
set_property(CACHE WarpX_PRECISION PROPERTY STRINGS ${WarpX_PRECISION_VALUES})
if(NOT WarpX_PRECISION IN_LIST WarpX_PRECISION_VALUES)
    message(FATAL_ERROR "WarpX_PRECISION (${WarpX_PRECISION}) must be one of ${WarpX_PRECISION_VALUES}")
endif()

set(WarpX_COMPUTE_VALUES NOACC OMP CUDA DPCPP) # HIP
set(WarpX_COMPUTE OMP CACHE STRING "On-node, accelerated computing backend (NOACC/OMP/CUDA/DPCPP)")
set_property(CACHE WarpX_COMPUTE PROPERTY STRINGS ${WarpX_COMPUTE_VALUES})
if(NOT WarpX_COMPUTE IN_LIST WarpX_COMPUTE_VALUES)
    message(FATAL_ERROR "WarpX_PRECISION (${WarpX_COMPUTE}) must be one of ${WarpX_COMPUTE_VALUES}")
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
    if(NOT ENABLE_CUDA)
        find_package(PkgConfig REQUIRED QUIET)
        if(WarpX_PRECISION STREQUAL "double")
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
# executable
add_executable(WarpX)
add_executable(WarpX::WarpX ALIAS WarpX)

# own headers
target_include_directories(WarpX PRIVATE
    $<BUILD_INTERFACE:${WarpX_SOURCE_DIR}/Source>
)

# if we include <AMReX_buildInfo.H> we will need to call:
include(AMReXBuildInfo)
generate_buildinfo(WarpX "${WarpX_SOURCE_DIR}")

# add sources
target_sources(WarpX
  PRIVATE
    Source/main.cpp
    Source/WarpX.cpp
)
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
#add_subdirectory(Source/Python)
add_subdirectory(Source/Utils)

# C++ properties: at least a C++14 capable compiler is needed
target_compile_features(WarpX PUBLIC cxx_std_14)
set_target_properties(WarpX PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

# link dependencies
target_link_libraries(WarpX PUBLIC WarpX::thirdparty::AMReX)

if(WarpX_PSATD)
    if(ENABLE_CUDA)
        # CUDA_ADD_CUFFT_TO_TARGET(WarpX)
        target_link_libraries(WarpX PUBLIC cufft)
    else()
        if(WarpX_PRECISION STREQUAL "double")
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
if(ENABLE_CUDA)
    setup_target_for_cuda_compilation(WarpX)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
        target_compile_features(WarpX PUBLIC cuda_std_14)
        set_target_properties(WarpX PROPERTIES
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
    target_compile_definitions(WarpX PRIVATE WARPX_DIM_3D)
elseif(WarpX_DIMS STREQUAL 2)
    target_compile_definitions(WarpX PRIVATE WARPX_DIM_XZ)
elseif(WarpX_DIMS STREQUAL RZ)
    target_compile_definitions(WarpX PRIVATE WARPX_DIM_RZ)
endif()

if(WarpX_OPENPMD)
    target_compile_definitions(WarpX PRIVATE WARPX_USE_OPENPMD)
endif()

if(WarpX_QED)
    target_compile_definitions(WarpX PRIVATE WARPX_QED)
    if(WarpX_QED_TABLE_GEN)
        target_compile_definitions(WarpX PRIVATE WarpX_QED_TABLE_GEN)
    endif()
endif()

if(WarpX_PSATD)
    target_compile_definitions(WarpX PRIVATE WARPX_USE_PSATD)
endif()

target_compile_definitions(WarpX PRIVATE
    WARPX_PARSER_DEPTH=${WarpX_PARSER_DEPTH})


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
set(WarpX_INSTALL_TARGET_NAMES WarpX)

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
