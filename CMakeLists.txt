# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.14.0)
project(WarpX VERSION 0.20.5)

include(${WarpX_SOURCE_DIR}/cmake/WarpXFunctions.cmake)


# Output Directories ##########################################################
#
# temporary build directories
set_default_build_dirs()

# install directories
set_default_install_dirs()


# Options and Variants ########################################################
#
warpx_option(QED                "PICSAR QED (requires Boost)"              OFF)
warpx_option(OPENPMD            "openPMD I/O"                              OFF)
#warpx_option(PYTHON             "Enable Python bindings"                  AUTO)
# PSATD (fftw; RZ:+blaspp+lapackpp), PICSAR/QED, hdf5, sensei, conduit, ascent, BL_ASSERT et al. (ON)

option(WarpX_amrex_internal    "Download & build AMReX" ON)
option(WarpX_picsar_internal   "Download & build PICSAR" ON)
option(WarpX_openpmd_internal  "Download & build openPMD-api" ON)

# change the default build type to RelWithDebInfo (or Release) instead of Debug
set_default_build_type("RelWithDebInfo")


# Dependencies ################################################################
#

# AMReX
#   builds AMReX from source (default) or finds an existing install
set(WarpX_amrex_repo "https://github.com/AMReX-Codes/amrex.git"
    CACHE STRING
    "Repository URI to pull and build AMReX from if(WarpX_amrex_internal)")
set(WarpX_amrex_branch "development"
    CACHE STRING
    "Repository branch for WarpX_amrex_repo if(WarpX_amrex_internal)")
include(${WarpX_SOURCE_DIR}/cmake/dependencies/AMReX.cmake)
#   suppress warnings in AMReX headers (use -isystem instead of -I)
make_third_party_includes_system(AMReX::amrex AMReX)

# PICSAR
#   builds PICSAR from source
set(WarpX_picsar_repo "https://bitbucket.org/berkeleylab/picsar.git"
    CACHE STRING
    "Repository URI to pull and build PICSAR from if(WarpX_picsar_internal)")
set(WarpX_picsar_branch "master"
    CACHE STRING
    "Repository branch for WarpX_picsar_repo if(WarpX_picsar_internal)")
include(${WarpX_SOURCE_DIR}/cmake/dependencies/PICSAR.cmake)

# openPMD
#   builds openPMD-api from source (default) or finds an existing install
set(WarpX_openpmd_repo "https://github.com/openPMD/openPMD-api.git"
    CACHE STRING
    "Repository URI to pull and build openPMD-api from if(WarpX_openpmd_internal)")
set(WarpX_openpmd_branch "dev"
    CACHE STRING
    "Repository branch for WarpX_openpmd_repo if(WarpX_openpmd_internal)")
include(${WarpX_SOURCE_DIR}/cmake/dependencies/openPMD.cmake)

# external library: FFTW (optional)
# external library: BLASPP LAPACKPP (optional)


# Targets #####################################################################
#
# executable
add_executable(WarpX)
add_executable(WarpX::WarpX ALIAS WarpX)

# own headers
target_include_directories(WarpX PUBLIC
    $<BUILD_INTERFACE:${WarpX_SOURCE_DIR}/Source>
    $<BUILD_INTERFACE:${WarpX_BINARY_DIR}/Source>
    $<INSTALL_INTERFACE:Source>
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
#add_subdirectory(Source/FortranInterface)
add_subdirectory(Source/Initialization)
add_subdirectory(Source/Laser)
add_subdirectory(Source/Parallelization)
add_subdirectory(Source/Parser)
add_subdirectory(Source/Particles)
#add_subdirectory(Source/Python)
add_subdirectory(Source/Utils)

# C++ properties: at least a C++14 capable compiler is needed
target_compile_features(WarpX
    PUBLIC cxx_std_14
)
set_target_properties(WarpX PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
)

# link dependencies
# FIXME use HiPACE::thirdparty::AMReX as soon as
#       missing _deps/fetchedamrex-build/mod_files is solved
target_link_libraries(WarpX PUBLIC AMReX::amrex)
if(WarpX_HAVE_OPENPMD)
    target_compile_definitions(WarpX PUBLIC WARPX_USE_OPENPMD)
    target_link_libraries(WarpX PUBLIC openPMD::openPMD)
endif()
if(WarpX_HAVE_QED)
    target_compile_definitions(WarpX PUBLIC WARPX_QED)
    # FIXME WARPX_QED_TABLE_GEN
    target_link_libraries(WarpX PUBLIC pxrmp_qed)
endif()

# AMReX helper function: propagate CUDA specific target & source properties
if(ENABLE_CUDA)
    setup_target_for_cuda_compilation(WarpX)
endif()


# Defines #####################################################################
#
target_compile_definitions(WarpX PRIVATE -DWARPX_DIM_3D)
target_compile_definitions(WarpX PRIVATE -DDWARPX_DIM_XZ)


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
# this defined the variable BUILD_TESTING which is ON by default
include(CTest)

#if(BUILD_TESTING)
#    enable_testing()
#
#    add_test(...)
#endif()


# Status Summary for Build Options ############################################
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
