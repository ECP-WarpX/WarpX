# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.18.0)
project(WarpX VERSION 22.05)

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
# AMReX 21.06+ supports CUDA_ARCHITECTURES with CMake 3.20+
# CMake 3.18+: CMAKE_CUDA_ARCHITECTURES
# https://cmake.org/cmake/help/latest/policy/CMP0104.html
if(POLICY CMP0104)
    cmake_policy(SET CMP0104 OLD)
endif()

# We use simple syntax in cmake_dependent_option, so we are compatible with the
# extended syntax in CMake 3.22+
# https://cmake.org/cmake/help/v3.22/policy/CMP0127.html
if(POLICY CMP0127)
    cmake_policy(SET CMP0127 NEW)
endif()


# C++ Standard in Superbuilds #################################################
#
# This is the easiest way to push up a C++17 requirement for AMReX, PICSAR and
# openPMD-api until they increase their requirement.
set_cxx17_superbuild()


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
include(CMakeDependentOption)
option(WarpX_APP           "Build the WarpX executable application"     ON)
option(WarpX_ASCENT        "Ascent in situ diagnostics"                 OFF)
option(WarpX_EB            "Embedded boundary support"                  OFF)
cmake_dependent_option(WarpX_GPUCLOCK
                           "Add GPU kernel timers (cost function)"      ON
                           "WarpX_COMPUTE STREQUAL CUDA OR WarpX_COMPUTE STREQUAL HIP" OFF)
option(WarpX_LIB           "Build WarpX as a shared library"            OFF)
option(WarpX_MPI           "Multi-node support (message-passing)"       ON)
option(WarpX_OPENPMD       "openPMD I/O (HDF5, ADIOS)"                  ON)
option(WarpX_PSATD         "spectral solver support"                    OFF)
option(WarpX_SENSEI        "SENSEI in situ diagnostics"                 OFF)
option(WarpX_QED           "QED support (requires PICSAR)"                    ON)
option(WarpX_QED_TABLE_GEN "QED table generation (requires PICSAR and Boost)" OFF)

set(WarpX_DIMS_VALUES 1 2 3 RZ)
set(WarpX_DIMS 3 CACHE STRING "Simulation dimensionality (1/2/3/RZ)")
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

option(WarpX_amrex_internal                    "Download & build AMReX" ON)

# change the default build type to Release (or RelWithDebInfo) instead of Debug
set_default_build_type("Release")

# Option to enable interprocedural optimization
# (also know as "link-time optimization" or "whole program optimization")
option(WarpX_IPO                                "Compile WarpX with interprocedural optimization (will take more time)" OFF)

# note: we could skip this if we solely build WarpX_APP, but if we build a
# shared WarpX library or a third party, like ImpactX, uses ablastr in a
# shared library (e.g., for Python bindings), then we need relocatable code.
option(ABLASTR_POSITION_INDEPENDENT_CODE
    "Build ABLASTR with position independent code" ${WarpX_LIB})
mark_as_advanced(ABLASTR_POSITION_INDEPENDENT_CODE)

# this defined the variable BUILD_TESTING which is ON by default
#include(CTest)


# Dependencies ################################################################
#

# AMReX
#   builds AMReX from source (default) or finds an existing install
include(${WarpX_SOURCE_DIR}/cmake/dependencies/AMReX.cmake)
#   suppress warnings in AMReX headers (use -isystem instead of -I)
warpx_make_third_party_includes_system(AMReX::amrex AMReX)

# PICSAR
#   builds PICSAR from source
include(${WarpX_SOURCE_DIR}/cmake/dependencies/PICSAR.cmake)

# openPMD
#   builds openPMD-api from source (default) or finds an existing install
include(${WarpX_SOURCE_DIR}/cmake/dependencies/openPMD.cmake)

# PSATD
include(${WarpX_SOURCE_DIR}/cmake/dependencies/FFT.cmake)
if(WarpX_PSATD)
    # BLASPP and LAPACKPP
    if(WarpX_DIMS STREQUAL RZ)
        find_package(blaspp CONFIG REQUIRED)
        find_package(lapackpp CONFIG REQUIRED)
        find_package(OpenMP REQUIRED)  # pulled by the two above
    endif()
endif()


# Targets #####################################################################
#
# collect all objects for compilation
add_library(WarpX OBJECT)
add_library(ablastr)

# ABLASTR library
set(_BUILDINFO_SRC ablastr)
set(_ALL_TARGETS WarpX ablastr)
add_library(WarpX::ablastr ALIAS ablastr)

# executable application
#   note: we currently avoid a dependency on a core library
#         for simpler usage, but could make this an option
if(WarpX_APP)
    add_executable(app)
    add_executable(WarpX::app ALIAS app)
    target_link_libraries(app PRIVATE WarpX ablastr)
    set(_BUILDINFO_SRC app)
    list(APPEND _ALL_TARGETS app)
endif()

# link into a shared library
if(WarpX_LIB)
    add_library(shared MODULE)
    add_library(WarpX::shared ALIAS shared)
    target_link_libraries(shared PUBLIC WarpX ablastr)
    set(_BUILDINFO_SRC shared)
    list(APPEND _ALL_TARGETS shared)

    set_target_properties(WarpX shared PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
    set(ABLASTR_POSITION_INDEPENDENT_CODE ON CACHE BOOL
        "Build ABLASTR with position independent code" FORCE)
endif()

# ABLASTR library (static or shared)
set_target_properties(ablastr PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)
if(ABLASTR_POSITION_INDEPENDENT_CODE)
    set_target_properties(ablastr PROPERTIES
        POSITION_INDEPENDENT_CODE ON
    )
endif()

# own headers
target_include_directories(WarpX PUBLIC
    $<BUILD_INTERFACE:${WarpX_SOURCE_DIR}/Source>
    $<BUILD_INTERFACE:${WarpX_BINARY_DIR}/Source>
)
target_include_directories(ablastr PUBLIC
    # future: own directory root
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

add_subdirectory(Source/ablastr)
add_subdirectory(Source/BoundaryConditions)
add_subdirectory(Source/Diagnostics)
add_subdirectory(Source/EmbeddedBoundary)
add_subdirectory(Source/Evolve)
add_subdirectory(Source/FieldSolver)
add_subdirectory(Source/Filter)
add_subdirectory(Source/Initialization)
add_subdirectory(Source/Laser)
add_subdirectory(Source/Parallelization)
add_subdirectory(Source/Particles)
add_subdirectory(Source/Python)
add_subdirectory(Source/Utils)

# C++ properties: at least a C++17 capable compiler is needed
foreach(warpx_tgt IN LISTS _ALL_TARGETS)
    target_compile_features(${warpx_tgt} PUBLIC cxx_std_17)
endforeach()
set_target_properties(${_ALL_TARGETS} PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

# Interprocedural optimization
if(WarpX_IPO)
    enable_IPO("${_ALL_TARGETS}")
endif()

# link dependencies
target_link_libraries(ablastr PUBLIC WarpX::thirdparty::AMReX)
target_link_libraries(WarpX PUBLIC ablastr)

if(WarpX_PSATD)
    target_link_libraries(ablastr PUBLIC WarpX::thirdparty::FFT)
    if(WarpX_DIMS STREQUAL RZ)
        target_link_libraries(ablastr PUBLIC blaspp)
        target_link_libraries(ablastr PUBLIC lapackpp)
    endif()
endif()

if(WarpX_OPENPMD)
    target_compile_definitions(ablastr PUBLIC WARPX_USE_OPENPMD)
    target_link_libraries(ablastr PUBLIC openPMD::openPMD)
endif()

if(WarpX_QED)
    target_compile_definitions(ablastr PUBLIC WARPX_QED)
    if(WarpX_QED_TABLE_GEN)
        target_compile_definitions(ablastr PUBLIC WARPX_QED_TABLE_GEN)
    endif()
    target_link_libraries(ablastr PUBLIC PXRMP_QED::PXRMP_QED)
endif()

# AMReX helper function: propagate CUDA specific target & source properties
if(WarpX_COMPUTE STREQUAL CUDA)
    foreach(warpx_tgt IN LISTS _ALL_TARGETS)
        setup_target_for_cuda_compilation(${warpx_tgt})
    endforeach()
    foreach(warpx_tgt IN LISTS _ALL_TARGETS)
        target_compile_features(${warpx_tgt} PUBLIC cuda_std_17)
    endforeach()
    set_target_properties(${_ALL_TARGETS} PROPERTIES
        CUDA_EXTENSIONS OFF
        CUDA_STANDARD_REQUIRED ON
    )
endif()

# avoid building all object files if we are only used as ABLASTR library
if(NOT WarpX_APP AND NOT WarpX_LIB)
    set_target_properties(WarpX PROPERTIES
        EXCLUDE_FROM_ALL 1
        EXCLUDE_FROM_DEFAULT_BUILD 1
    )
endif()

# fancy binary name for build variants
set_warpx_binary_name()


# Defines #####################################################################
#
if(WarpX_DIMS STREQUAL 3)
    target_compile_definitions(ablastr PUBLIC WARPX_DIM_3D WARPX_ZINDEX=2)
elseif(WarpX_DIMS STREQUAL 2)
    target_compile_definitions(ablastr PUBLIC WARPX_DIM_XZ WARPX_ZINDEX=1)
elseif(WarpX_DIMS STREQUAL 1)
    target_compile_definitions(ablastr PUBLIC WARPX_DIM_1D_Z WARPX_ZINDEX=0)
elseif(WarpX_DIMS STREQUAL RZ)
    target_compile_definitions(ablastr PUBLIC WARPX_DIM_RZ WARPX_ZINDEX=1)
endif()

if(WarpX_GPUCLOCK)
    target_compile_definitions(ablastr PUBLIC WARPX_USE_GPUCLOCK)
endif()

if(WarpX_OPENPMD)
    target_compile_definitions(ablastr PUBLIC WARPX_USE_OPENPMD)
endif()

if(WarpX_QED)
    target_compile_definitions(ablastr PUBLIC WARPX_QED)
    if(WarpX_QED_TABLE_GEN)
        target_compile_definitions(ablastr PUBLIC WarpX_QED_TABLE_GEN)
    endif()
endif()

if(WarpX_PSATD)
    target_compile_definitions(ablastr PUBLIC WARPX_USE_PSATD)
endif()

# <cmath>: M_PI
if(WIN32)
    target_compile_definitions(ablastr PUBLIC _USE_MATH_DEFINES)
endif()


# Warnings ####################################################################
#
set_cxx_warnings()


# Generate Configuration and .pc Files ########################################
#
get_source_version(WarpX ${WarpX_SOURCE_DIR})
configure_file(
    ${WarpX_SOURCE_DIR}/Source/Utils/WarpXVersion.H.in
    ${WarpX_BINARY_DIR}/Source/Utils/WarpXVersion.H
    @ONLY
)

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

# simplified library alias
# this is currently expected by Python bindings
if(WarpX_LIB)
    if(WarpX_DIMS STREQUAL RZ)
        set(lib_suffix "rz")
    else()
        set(lib_suffix "${WarpX_DIMS}d")
    endif()
    if(WIN32)
        set(mod_ext "dll")
    else()
        set(mod_ext "so")
    endif()
    if(IS_ABSOLUTE ${CMAKE_INSTALL_LIBDIR})
        set(ABS_INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR})
    else()
        set(ABS_INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
    endif()
    install(CODE "file(CREATE_LINK
        $<TARGET_FILE_NAME:shared>
        ${ABS_INSTALL_LIB_DIR}/libwarpx.${lib_suffix}.${mod_ext}
        COPY_ON_ERROR SYMBOLIC)")
endif()

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


# pip helpers for the pywarpx package #########################################
#
if(WarpX_LIB)
    set(PYINSTALLOPTIONS "" CACHE STRING
        "Additional parameters to pass to `pip install`")

    # build the wheel by re-using the shared library we build
    add_custom_target(pip_wheel
        ${CMAKE_COMMAND} -E rm -f ${WarpX_BINARY_DIR}/pywarpx*whl
        COMMAND
            ${CMAKE_COMMAND} -E env PYWARPX_LIB_DIR=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
                python3 -m pip wheel -v --no-build-isolation ${WarpX_SOURCE_DIR}
        WORKING_DIRECTORY
            ${WarpX_BINARY_DIR}
        DEPENDS
            shared
    )

    # this will also upgrade/downgrade dependencies, e.g., when the version of picmistandard changes
    add_custom_target(pip_install_requirements
        python3 -m pip install ${PYINSTALLOPTIONS} -r ${WarpX_SOURCE_DIR}/requirements.txt
        WORKING_DIRECTORY
            ${WarpX_BINARY_DIR}
    )

    # We force-install because in development, it is likely that the version of
    # the package does not change, but the code did change. We need --no-deps,
    # because otherwise pip would also force reinstall all dependencies.
    add_custom_target(pip_install
        ${CMAKE_COMMAND} -E env WARPX_MPI=${WarpX_MPI}
            python3 -m pip install --force-reinstall --no-deps ${PYINSTALLOPTIONS} ${WarpX_BINARY_DIR}/pywarpx*whl
        WORKING_DIRECTORY
            ${WarpX_BINARY_DIR}
        DEPENDS
            shared pip_wheel pip_install_requirements
    )
endif()


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
