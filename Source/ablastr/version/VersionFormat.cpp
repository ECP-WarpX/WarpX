/* Copyright 2022 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "VersionFormat.H"
#include "WarpX.H"

#include <AMReX_Print.H>

#ifdef __CUDACC_VER_MAJOR__
#   include <cuda.h>
#endif
#ifdef WARPX_USE_OPENPMD
#   include "openPMD/openPMD.hpp"
#endif
#if defined(AMREX_USE_MPI)
#   include <mpi.h>
#endif

#include <map>
#include <string>
#include <sstream>

// Stringizes arguments after they have been expanded.
#define ABLASTR_STRINGIZE(...) ABLASTR_STRINGIZE2(__VA_ARGS__)
#define ABLASTR_STRINGIZE2(...) #__VA_ARGS__


void
ablastr::version::printSoftwareVersions()
{
    std::string const versionNotFound("unknown");
    std::string const versionDisabled("disabled");

    std::stringstream warpxCompute;
#if defined(AMREX_USE_CUDA)
    warpxCompute << "CUDA";
#elif defined(AMREX_USE_HIP)
    warpxCompute << "HIP";
#elif defined(AMREX_USE_DPCPP)
    warpxCompute << "SYCL/DPC++";
#elif defined(AMREX_USE_OMP)
    warpxCompute << "OpenMP";
#else
    warpxCompute << "NOACC";
#endif

    std::stringstream picsar;
#ifdef PICSAR_GIT_VERSION
    picsar << std::string(PICSAR_GIT_VERSION);
#elif !defined(WARPX_QED)
    picsar << versionDisabled;
#endif
    if( picsar.str().empty() )
        picsar << versionNotFound;

    std::stringstream picsarTableGen;
#if defined(WARPX_QED_TABLE_GEN)
    picsarTableGen << " (with table generation)";
#endif

#if 0 && defined(WARPX_QED_TABLE_GEN)
    // commented, since only transitively from PICSAR
    std::stringstream boost;
    boost << int(BOOST_VERSION / 100000) << "." << int(BOOST_VERSION / 100 % 1000) << "."
          << int(BOOST_VERSION % 100);
#endif

    std::stringstream buildType;
#ifdef CMAKE_BUILD_TYPE
    buildType << ABLASTR_STRINGIZE(CMAKE_BUILD_TYPE);
#else
    buildType << versionNotFound;
#endif

    std::stringstream os;
#ifdef CMAKE_SYSTEM
    os << ABLASTR_STRINGIZE(CMAKE_SYSTEM);
#else
    os << versionNotFound;
#endif

    std::stringstream arch;
#ifdef CMAKE_SYSTEM_PROCESSOR
    arch << ABLASTR_STRINGIZE(CMAKE_SYSTEM_PROCESSOR);
#else
    arch << versionNotFound;
#endif

    std::stringstream cxx;
#ifdef CMAKE_CXX_COMPILER_ID
    cxx << ABLASTR_STRINGIZE(CMAKE_CXX_COMPILER_ID);
#else
    cxx << versionNotFound;
#endif

    std::stringstream cxxVersion;
#ifdef CMAKE_CXX_COMPILER_VERSION
    cxxVersion << ABLASTR_STRINGIZE(CMAKE_CXX_COMPILER_VERSION);
#else
    cxxVersion << versionNotFound;
#endif

    std::stringstream cmake;
#ifdef CMAKE_VERSION
    cmake << ABLASTR_STRINGIZE(CMAKE_VERSION);
#else
    cmake << versionNotFound;
#endif

    std::stringstream computeVersion;
#ifdef __CUDACC_VER_MAJOR__
    computeVersion << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#endif

#if defined(AMREX_USE_OMP)
    std::map< int, std::string > mapOpenMP = {
        {200505, "2.5"},
        {200805, "3.0"},
        {201107, "3.1"},
        {201307, "4.0"},
        {201511, "4.5"},
        {201811, "5.0"}
    };

    if (mapOpenMP.count(_OPENMP) == 1u)
        computeVersion << mapOpenMP[_OPENMP];
    else
        computeVersion << _OPENMP;
#endif

    // TODO: HIP, SYCL

    std::stringstream mpiStandard;
    std::stringstream mpiFlavor;
    std::stringstream mpiFlavorVersion;
#if defined(AMREX_USE_MPI)
    mpiStandard << MPI_VERSION << "." << MPI_SUBVERSION;
#   if defined(OMPI_MAJOR_VERSION)
    // includes derivatives such as Bull MPI, Sun, ...
    mpiFlavor << "OpenMPI";
    mpiFlavorVersion << OMPI_MAJOR_VERSION << "." << OMPI_MINOR_VERSION << "." << OMPI_RELEASE_VERSION;
#   elif defined(MPICH_VERSION)
    /* includes MPICH2 and MPICH3 and
     * derivates such as IBM, Cray, MS, Intel, MVAPICH(2), ... */
    mpiFlavor << "MPICH";
    mpiFlavorVersion << MPICH_VERSION;
#   else
    mpiFlavor << versionNotFound;
    mpiFlavorVersion << versionNotFound;
#   endif
#else
    mpiStandard << versionDisabled;
    mpiFlavor << versionDisabled;
    mpiFlavorVersion << versionDisabled;
#endif

#ifdef WARPX_USE_OPENPMD
    std::string openPMD = openPMD::getVersion();
    std::string openPMDstandard = openPMD::getStandard();
    std::string openPMDext;
    auto const openPMDextVec = openPMD::getFileExtensions();
    for (auto it = openPMDextVec.cbegin();
         it != openPMDextVec.cend(); it++)
    {
        openPMDext.append(*it);

        if (openPMDextVec.size() > 1)
            if (openPMDextVec.end() - 1 != it)
                openPMDext.append(",");
    }
#else
    std::string openPMD = versionDisabled;
    std::string openPMDstandard = versionDisabled;
    std::string openPMDext = versionDisabled;
#endif

    // TODO: FFTW
    // TODO: BLAS++/LAPACK++
    // TODO: Ascent
    // TODO: SENSEI

    // command-line output
    //   AMReX: automatically printed
    //   app: add manually before this call
    amrex::Print() << "Build options:" << "\n";
    amrex::Print() << "  Build type: " << buildType.str() << "\n";
    amrex::Print() << "  OS:         " << os.str() << "\n";
    amrex::Print() << "  arch:       " << arch.str() << "\n";
    amrex::Print() << "  CXX:        " << cxx.str() << " (" << cxxVersion.str() << ")" << "\n";
    amrex::Print() << "  compute:    " << warpxCompute.str() << " (" << computeVersion.str() << ")" << "\n";
    amrex::Print() << "Third party:" << std::endl;
    amrex::Print() << "  CMake:      " << cmake.str() << "\n";
    amrex::Print() << "  MPI:        " << std::endl
                   << "    standard: " << mpiStandard.str() << "\n"
                   << "    flavor:   " << mpiFlavor.str() << " (" << mpiFlavorVersion.str() << ")\n";
    amrex::Print() << "  PICSAR:     " << picsar.str() << picsarTableGen.str() << "\n";
    amrex::Print() << "  openPMD:    " << openPMD
                   << " (standard: " << openPMDstandard
                   << "; files: " << openPMDext << ")\n";
    // only transitive in PICSAR QED with tables
    // amrex::Print() << "  Boost:      " << boost.str() << "\n";
    // TODO: FFTW
    // TODO: BLAS++/LAPACK++
    // TODO: Ascent
    // TODO: SENSEI
}
