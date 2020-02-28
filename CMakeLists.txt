# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.14.0)
project(WarpX VERSION 0.20.2)

# list(APPEND CMAKE_MODULE_PATH "${WarpX_SOURCE_DIR}/cmake")


# Project structure ###########################################################
#
# temporary build directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    CACHE PATH "Build directory for archives")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    CACHE PATH "Build directory for libraries")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    CACHE PATH "Build directory for binaries")

# install directories
include(GNUInstallDirs)
set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/WarpX"
    CACHE PATH "CMake config package location for installed targets")
if(WIN32)
    set(CMAKE_INSTALL_LIBDIR Lib
        CACHE PATH "Object code libraries")
    set_property(CACHE CMAKE_INSTALL_CMAKEDIR PROPERTY VALUE "cmake")
endif()


# Options and Variants ########################################################
#
function(WarpX_option name description default)
    set(WarpX_USE_${name} ${default} CACHE STRING "${description}")
    set_property(CACHE WarpX_USE_${name} PROPERTY
        STRINGS "ON;TRUE;AUTO;OFF;FALSE"
    )
    if(WarpX_HAVE_${name})
        set(WarpX_HAVE_${name} TRUE)
    else()
        set(WarpX_HAVE_${name})
    endif()
    # list of all possible options
    set(WarpX_CONFIG_OPTIONS ${WarpX_CONFIG_OPTIONS} ${name} PARENT_SCOPE)
endfunction()

WarpX_option(MPI            "Parallel, Multi-Node for clusters"         AUTO)
WarpX_option(openPMD        "openPMD I/O"                               AUTO)
WarpX_option(PYTHON         "Enable Python bindings"                    AUTO)
# PSATD (fftw; RZ:+blas+lapack), PICSAR/QED, hdf5, sensei, conduit, ascent, BL_ASSERT et al. (ON)

option(WarpX_DIM             "Spatial dimensions"      3)
option(WarpX_USE_INTERNAL_AMReX  "Download & build AMReX"  ON)
option(WarpX_USE_INTERNAL_PICSAR "Download & build PICSAR" OFF)
# CUDA/HIP
# OpenMP
# SP or DP; SP or DP for particles
# tiny profiler

set(CMAKE_CONFIGURATION_TYPES "Release;Debug;MinSizeRel;RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the build type, e.g. Release or Debug." FORCE)
endif()

#include(CTest)
# automatically defines: BUILD_TESTING, default is ON

# Dependencies ################################################################
#
# external library: MPI (optional)
#if(WarpX_USE_MPI STREQUAL AUTO)
#    find_package(MPI)
#    if(MPI_FOUND)
#        set(WarpX_HAVE_MPI TRUE)
#    else()
#        set(WarpX_HAVE_MPI FALSE)
#    endif()
#elseif(WarpX_USE_MPI)
#    find_package(MPI REQUIRED)
#    set(WarpX_HAVE_MPI TRUE)
#else()
#    set(WarpX_HAVE_MPI FALSE)
#endif()

# external library: AMReX (required)

# ENABLE_MPI = ON
# ENABLE_OMP = OFF
# ENABLE_CUDA = OFF
# ENABLE_DP = ON
# ENABLE_EB = OFF
# ENABLE_FORTRAN_INTERFACES = OFF
# ENABLE_LINEAR_SOLVERS = ON
# ENABLE_AMRDATA = OFF
# ENABLE_PARTICLES = ON
# ENABLE_DP_PARTICLES = ON
# ENABLE_SENSEI_INSITU = OFF
# ENABLE_CONDUIT = OFF
# ENABLE_SUNDIALS = OFF
# ENABLE_HYPRE = OFF
# ENABLE_PETSC =
# ENABLE_FPE = OFF
# ENABLE_ASSERTIONS = OFF
# ENABLE_BASE_PROFILE = OFF
# ENABLE_TINY_PROFILE = OFF
# ENABLE_TRACE_PROFILE = OFF
# ENABLE_MEM_PROFILE = OFF
# ENABLE_COMM_PROFILE = OFF
# ENABLE_BACKTRACE = OFF
# ENABLE_PROFPARSER = OFF
# TP_PROFILE = None
# ENABLE_CUDA = OFF
# ENABLE_ACC = OFF

if(WarpX_USE_INTERNAL_AMReX)
    include(FetchContent)
    set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
    #set(BUILD_TESTING OFF)
    #set(BUILD_EXAMPLES OFF)
    set(ENABLE_TUTORIALS OFF CACHE INTERNAL "")
    set(ENABLE_PARTICLES ON CACHE INTERNAL "")
    set(ENABLE_TINY_PROFILE ON CACHE INTERNAL "")
    set(AMReX_DIM ${WarpX_DIM} CACHE INTERNAL "")
    #set(AMReX_ENABLE_PIC ON CACHE INTERNAL "")
    # TEMPORARY until fixed upstream
    set(ENABLE_FORTRAN_INTERFACES ON CACHE INTERNAL "")
    #set(ENABLE_FORTRAN ON CACHE INTERNAL "")
    FetchContent_Declare(AMReX
      GIT_REPOSITORY "https://github.com/AMReX-Codes/amrex.git"
      GIT_TAG        "development"
      #GIT_REPOSITORY "/home/axel/src/amrex/.git"
      #GIT_REPOSITORY "https://github.com/ax3l/amrex.git"
      #GIT_TAG        "topic-cudaTargetFlags"
      BUILD_IN_SOURCE 0
    )
    FetchContent_MakeAvailable(AMReX)

    # FIXME hackity-hack path
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_BINARY_DIR}/_deps/amrex-src/Tools/CMake")

    message(WARNING "VERSION: ${AMREX_GIT_VERSION}")
    # message(WARNING "VERSION: ${AMReX_VERSION}")  # FIXME missing upstream

    message(STATUS "AMReX: Using INTERNAL version ${AMReX_VERSION}")
else()
    find_package(AMReX 3.7.0 CONFIG REQUIRED)
    message(STATUS "AMReX: Found version ${AMReX_VERSION}")
endif()

if(ENABLE_CUDA)
    enable_language(CUDA)
    include(AMReX_SetupCUDA)
    # set_source_files_properties(main.cpp PROPERTIES LANGUAGE CUDA)
endif()

add_library(WarpX::thirdparty::AMReX INTERFACE IMPORTED)
target_link_libraries(WarpX::thirdparty::AMReX
    INTERFACE AMReX::amrex)

# external library: PICSAR (required)
# external library: FFTW (optional)
# external library: openPMD-api (optional)


# Targets #####################################################################
#
set(CORE_SOURCE
    Source/main.cpp
    Source/WarpX.cpp
)

#if(ENABLE_CUDA)
#    set_source_files_properties(Source/main.cpp Source/WarpX.cpp PROPERTIES LANGUAGE CUDA)
#endif()

# executable
add_executable(WarpX ${CORE_SOURCE})
add_executable(WarpX::WarpX ALIAS WarpX)

# own headers
target_include_directories(WarpX PUBLIC
    $<BUILD_INTERFACE:${WarpX_SOURCE_DIR}/Source>
#    $<BUILD_INTERFACE:${WarpX_BINARY_DIR}/Source>
#    $<INSTALL_INTERFACE:Source>
)

# more sources
add_subdirectory(Source/BoundaryConditions)
add_subdirectory(Source/Diagnostics)
add_subdirectory(Source/Evolve)
add_subdirectory(Source/FieldSolver)
add_subdirectory(Source/Filter)
#add_subdirectory(Source/FortranInterface)
add_subdirectory(Source/Initialization)
add_subdirectory(Source/Laser)
add_subdirectory(Source/Parallelization)
add_subdirectory(Source/Parser)
add_subdirectory(Source/Particles)
#add_subdirectory(Source/Python)
#add_subdirectory(Source/QED)
add_subdirectory(Source/Utils)

# we include <AMReX_buildInfo.H>
include(AMReXBuildInfo)  # FIXME this requires Python at build time...
# FIXME hackity-hack path
generate_buildinfo(WarpX "${PROJECT_BINARY_DIR}/_deps/amrex-src/")

# FIXME
#set_target_properties(WarpX
#   PROPERTIES
#     LINKER_LANGUAGE
#       Fortran
#)

# properties
target_compile_features(WarpX
    PUBLIC cxx_std_14
)
set_target_properties(WarpX PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

if(ENABLE_CUDA)
    set_cpp_sources_to_cuda_language( WarpX )
endif()


# FIXME hackity-hack path to PICSAR
target_include_directories(WarpX PRIVATE
    /home/axel/src/picsar/src/multi_physics/QED/src
)

#if(WarpX_HAVE_MPI)
#    # MPI targets: CMake 3.9+
#    # note: often the PUBLIC dependency to CXX is missing in C targets...
#    target_link_libraries(WarpX PUBLIC MPI::MPI_C MPI::MPI_CXX)
#endif()

# JSON Backend and User-Facing Runtime Options
#target_link_libraries(WarpX PRIVATE WarpX::thirdparty::AMReX)
#target_include_directories(WarpX SYSTEM PRIVATE
#    $<TARGET_PROPERTY:WarpX::thirdparty::AMReX,INTERFACE_INCLUDE_DIRECTORIES>)
#target_link_libraries(WarpX PRIVATE WarpX::thirdparty::AMReX)
# should not be needed... if we setup cuda via amrex helpers...
target_compile_options(
    WarpX::thirdparty::AMReX INTERFACE
      $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
      $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>)
target_link_libraries(WarpX PRIVATE WarpX::thirdparty::AMReX)

# HDF5 Backend
#if(WarpX_HAVE_HDF5)
#    target_link_libraries(WarpX PRIVATE ${HDF5_LIBRARIES})
#    target_include_directories(WarpX SYSTEM PRIVATE ${HDF5_INCLUDE_DIRS})
#    target_compile_definitions(WarpX PRIVATE ${HDF5_DEFINITIONS})
#endif()

# Defines #####################################################################
#
target_compile_definitions(WarpX PRIVATE -DWARPX_DIM_3D)
target_compile_definitions(WarpX PRIVATE -DDWARPX_DIM_XZ)

# DEFINES += -DWARPX_USE_PY
# DEFINES += -DWARPX_DIM_3D
# DEFINES += -DWARPX_DIM_RZ
# DEFINES += -DWARPX_DIM_XZ
# CXXFLAGS += -DWARPX_QED
# CXXFLAGS += -DWARPX_QED_TABLE_GEN
# DEFINES += -DWARPX_GIT_VERSION=\"$(WARPX_GIT_VERSION)\"
# DEFINES += -DWARPX_USE_OPENPMD
# DEFINES += -DWARPX_USE_PSATD
# DEFINES += -DWARPX_USE_PSATD_HYBRID
# DEFINES += -DWARPX_DO_ELECTROSTATIC
# DEFINES += -DWARPX_USE_HDF5


# Warnings ####################################################################
#
# TODO: LEGACY! Use CMake TOOLCHAINS instead!
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # list(APPEND CMAKE_CXX_FLAGS "-fsanitize=address") # address, memory, undefined
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    # set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
    # set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -fsanitize=address")

    # note: might still need a
    #   export LD_PRELOAD=libclang_rt.asan.so
    # or on Debian 9 with Clang 6.0
    #   export LD_PRELOAD=/usr/lib/llvm-6.0/lib/clang/6.0.0/lib/linux/libclang_rt.asan-x86_64.so:
    #                     /usr/lib/llvm-6.0/lib/clang/6.0.0/lib/linux/libclang_rt.ubsan_minimal-x86_64.so
    # at runtime when used with symbol-hidden code (e.g. pybind11 module)

    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wshadow -Woverloaded-virtual -Wextra-semi -Wunreachable-code")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wshadow -Woverloaded-virtual -Wunreachable-code")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # Warning C4503: "decorated name length exceeded, name was truncated"
    # Symbols longer than 4096 chars are truncated (and hashed instead)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4503")
    # Yes, you should build against the same C++ runtime and with same
    # configuration (Debug/Release). MSVC does inconvenient choices for their
    # developers, so be it. (Our Windows-users use conda-forge builds, which
    # are consistent.)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4251")
endif ()


# Generate Files with Configuration Options ###################################
#
# TODO configure a version.hpp
#configure_file(
#    ${WarpX_SOURCE_DIR}/include/WarpX/config.hpp.in
#    ${WarpX_BINARY_DIR}/include/WarpX/config.hpp
#    @ONLY
#)

#configure_file(
#    ${WarpX_SOURCE_DIR}/WarpXConfig.cmake.in
#    ${WarpX_BINARY_DIR}/WarpXConfig.cmake
#    @ONLY
#)

# get absolute paths to linked libraries
function(WarpXreclibs tgtname outname)
    get_target_property(PC_PRIVATE_LIBS_TGT ${tgtname} INTERFACE_LINK_LIBRARIES)
    foreach(PC_LIB IN LISTS PC_PRIVATE_LIBS_TGT)
       if(TARGET ${PC_LIB})
           WarpXreclibs(${PC_LIB} ${outname})
       else()
           if(PC_LIB)
               string(APPEND ${outname} " ${PC_LIB}")
           endif()
       endif()
    endforeach()
    set(${outname} ${${outname}} PARENT_SCOPE)
endfunction()

if(WarpX_HAVE_PKGCONFIG)
    WarpXreclibs(WarpX WarpX_PC_PRIVATE_LIBS)
    if(BUILD_SHARED_LIBS)
        set(WarpX_PC_STATIC false)
    else()
        set(WarpX_PC_STATIC true)
    endif()
    configure_file(
        ${WarpX_SOURCE_DIR}/WarpX.pc.in
        ${WarpX_BINARY_DIR}/WarpX.pc
        @ONLY
    )
endif()

include(CMakePackageConfigHelpers)
write_basic_package_version_file("WarpXConfigVersion.cmake"
    VERSION ${WarpX_VERSION}
    COMPATIBILITY SameMajorVersion
)


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
#if(WarpX_HAVE_PYTHON)
#    install(TARGETS WarpX.py
#        DESTINATION ${CMAKE_INSTALL_PYTHONDIR}
#    )
#endif()
#install(DIRECTORY "${WarpX_SOURCE_DIR}/include/WarpX"
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#    FILES_MATCHING PATTERN "*.hpp"
#)
#install(
#    FILES ${WarpX_BINARY_DIR}/include/WarpX/config.hpp
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/WarpX
#)

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
# pkg-config .pc file for depending legacy projects
#   This is for projects that do not use a build file generator, e.g.
#   because they compile manually on the command line or write their
#   Makefiles by hand.
if(WarpX_HAVE_PKGCONFIG)
    install(
        FILES       ${WarpX_BINARY_DIR}/WarpX.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    )
endif()


# Tests #######################################################################
#
#if(BUILD_TESTING)
#    enable_testing()
#
#    add_test(...)
#endif()


# Examples ####################################################################
#
if(BUILD_EXAMPLES)
    # ...
endif()


# Status Message for Build Options ############################################
#
message("")
message("WarpX build configuration:")
message("  Version: ${WarpX_VERSION}")
message("  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} "
                        "${CMAKE_CXX_COMPILER_VERSION} "
                        "${CMAKE_CXX_COMPILER_WRAPPER}")
message("    ${CMAKE_CXX_COMPILER}")
message("")
message("  Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message("        bin: ${CMAKE_INSTALL_BINDIR}")
message("        lib: ${CMAKE_INSTALL_LIBDIR}")
message("    include: ${CMAKE_INSTALL_INCLUDEDIR}")
message("      cmake: ${CMAKE_INSTALL_CMAKEDIR}")
if(WarpX_HAVE_PYTHON)
    message("     python: ${CMAKE_INSTALL_PYTHONDIR}")
endif()
message("")
message("  Build Type: ${CMAKE_BUILD_TYPE}")
#if(BUILD_SHARED_LIBS)
#    message("  Library: shared")
#else()
#    message("  Library: static")
#endif()
message("  Testing: ${BUILD_TESTING}")
#message("  Invasive Tests: ${WarpX_USE_INVASIVE_TESTS}")
#message("  Internal VERIFY: ${WarpX_USE_VERIFY}")
message("  Build Options:")

foreach(opt IN LISTS WarpX_CONFIG_OPTIONS)
  if(${WarpX_HAVE_${opt}})
    message("    ${opt}: ON")
  else()
    message("    ${opt}: OFF")
  endif()
endforeach()
message("    DIM: ${WarpX_DIM}")
message("")

