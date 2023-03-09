/* Copyright 2021-2022 The ImpactX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "pyWarpX.H"
#include "WarpX_py.H"

#include <WarpX.H>
#include <Initialization/WarpXAMReXInit.H>
#include <Utils/WarpXUtil.H>

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

    m.def("amrex_init",
        [](const py::list args) {
            amrex::Vector<std::string> cargs;
            amrex::Vector<char*> argv;

            // Populate the "command line"
            for (const auto& v: args)
                cargs.push_back(v.cast<std::string>());
            for (auto& v: cargs)
                argv.push_back(&v[0]);
            int argc = argv.size();

            // note: +1 since there is an extra char-string array element,
            //       that ANSII C requires to be a simple NULL entry
            //       https://stackoverflow.com/a/39096006/2719194
            argv.push_back(NULL);
            char** tmp = argv.data();

            const bool build_parm_parse = (cargs.size() > 1);
            // TODO: handle version with MPI
            return warpx_amrex_init(argc, tmp, build_parm_parse);
        }, py::return_value_policy::reference,
        "Initialize AMReX library");
    m.def("amrex_finalize", [] () {amrex::Finalize();},
        "Close out the amrex related data");
    m.def("warpx_finalize", &WarpX::ResetInstance,
        "Close out the WarpX related data");
    m.def("convert_lab_params_to_boost",  &ConvertLabParamsToBoost,
        "Convert input parameters from the lab frame to the boosted frame");
    m.def("read_BC_params", &ReadBCParams,
        "Read the boundary condition parametes and check for consistency");
    m.def("check_gridding_for_RZ_spectral", &CheckGriddingForRZSpectral,
        "Ensure that the grid is setup appropriately with using the RZ spectral solver");

    // Expose the python callback function installation and removal functions
    m.def("add_python_callback", &InstallPythonCallback);
    m.def("remove_python_callback", &ClearPythonCallback);
    m.def("execute_python_callback", &ExecutePythonCallback, py::arg("name"));
}
