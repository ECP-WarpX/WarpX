/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyWarpX.H"

#include <WarpX.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#if defined(AMREX_DEBUG) || defined(DEBUG)
#   include <cstdio>
#endif
#include <string>


namespace py = pybind11;
//using namespace warpx;

namespace warpx {
    struct Config {};
}

void init_WarpX (py::module& m)
{
    py::class_<WarpX> warpx(m, "WarpX");
    warpx
        .def(py::init<>())

        // from AmrCore->AmrMesh
        .def("Geom",
            //[](ImpactX const & ix, int const lev) { return ix.Geom(lev); },
            py::overload_cast< int >(&WarpX::Geom, py::const_),
            py::arg("lev")
        )
        .def("DistributionMap",
            [](WarpX const & ix, int const lev) { return ix.DistributionMap(lev); },
            //py::overload_cast< int >(&ImpactX::DistributionMap, py::const_),
            py::arg("lev")
        )
        .def("boxArray",
            [](WarpX const & ix, int const lev) { return ix.boxArray(lev); },
            //py::overload_cast< int >(&ImpactX::boxArray, py::const_),
            py::arg("lev")
        )
        .def("multifab",
            [](WarpX const & ix, std::string const multifab_name) {
                if (ix.multifab_map.count(multifab_name) > 0) {
                    return ix.multifab_map.at(multifab_name);
                } else {
                    throw std::runtime_error("The MultiFab '" + multifab_name + "' is unknown or is not allocated!");
                }
            },
            py::arg("multifab_name"),
            py::return_value_policy::reference_internal,
            "Return MultiFabs by name, e.g., 'Efield_aux[x][l=0]', 'Efield_cp[x][l=0]', ..."
        )
    ;

    py::class_<warpx::Config>(m, "Config")
//        .def_property_readonly_static(
//            "warpx_version",
//            [](py::object) { return Version(); },
//            "WarpX version")
        .def_property_readonly_static(
            "have_mpi",
            [](py::object){
#ifdef AMREX_USE_MPI
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "have_gpu",
            [](py::object){
#ifdef AMREX_USE_GPU
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "have_omp",
            [](py::object){
#ifdef AMREX_USE_OMP
                return true;
#else
                return false;
#endif
            })
        .def_property_readonly_static(
            "gpu_backend",
            [](py::object){
#ifdef AMREX_USE_CUDA
                return "CUDA";
#elif defined(AMREX_USE_HIP)
                return "HIP";
#elif defined(AMREX_USE_DPCPP)
                return "SYCL";
#else
                return py::none();
#endif
            })
        ;
}
