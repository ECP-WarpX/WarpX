/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyWarpX.H"

#include <WarpX.H>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#define CONCAT_NAME(PRE, SUF) PRE ## SUF

// see: CMakeLists.txt, setup.py and __init__.py
#if defined(WARPX_DIM_1D_Z)
#  define PYWARPX_MODULE_NAME CONCAT_NAME(warpx_pybind_, 1d)
#elif defined(WARPX_DIM_XZ)
#  define PYWARPX_MODULE_NAME CONCAT_NAME(warpx_pybind_, 2d)
#elif defined(WARPX_DIM_RZ)
#  define PYWARPX_MODULE_NAME CONCAT_NAME(warpx_pybind_, rz)
#elif defined(WARPX_DIM_3D)
#  define PYWARPX_MODULE_NAME CONCAT_NAME(warpx_pybind_, 3d)
#endif

namespace py = pybind11;
//using namespace warpx;


// forward declarations of exposed classes
void init_WarpX(py::module&);

PYBIND11_MODULE(PYWARPX_MODULE_NAME, m) {
    // make sure AMReX types are known
    py::module::import("amrex");

    m.doc() = R"pbdoc(
            warpx_pybind
            --------------
            .. currentmodule:: warpx_pybind_(1d|2d|3d|rz)

            .. autosummary::
               :toctree: _generate
               WarpX
    )pbdoc";

    // note: order from parent to child classes
    init_WarpX(m);

    // API runtime version
    //   note PEP-440 syntax: x.y.zaN but x.y.z.devN
#ifdef PYIMPACTX_VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(PYIMPACTX_VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif

    // authors
    m.attr("__author__") =
        "Jean-Luc Vay, David P. Grote, Maxence Thevenet, Remi Lehe, Andrew Myers, Weiqun Zhang, Axel Huebl, et al.";

    // API runtime build-time feature variants
    // m.attr("variants") = warpx::getVariants();
    // TODO allow to query runtime versions of all dependencies

    // license SPDX identifier
    m.attr("__license__") = "BSD-3-Clause-LBNL";

    // TODO broken numpy if not at least v1.15.0: raise warning
    // auto numpy = py::module::import("numpy");
    // auto npversion = numpy.attr("__version__");
    // std::cout << "numpy version: " << py::str(npversion) << std::endl;
}
