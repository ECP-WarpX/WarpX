/* Copyright 2024 The WarpX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "Python/pyWarpX.H"

#include <ablastr/fields/MultiFabRegister.H>


void init_MultiFabRegister (py::module & m)
{
    using ablastr::fields::MultiFabRegister;

    py::class_<ablastr::fields::MultiFabRegister>(m, "MultiFabRegister")
        .def("alloc_init",
             &MultiFabRegister::alloc_init
             // ... py::arg("name")
             // "..."
        )
        .def("alloc_like",
             &MultiFabRegister::alloc_like,
             py::arg("other_name"),
             py::arg("other_level")
             // "..."
        )
        .def("has",
             &MultiFabRegister::has,
             py::arg("name"),
             py::arg("level")
             // "..."
        )
        .def("list",
             &MultiFabRegister::list
             // "..."
        )
        .def("erase",
             &MultiFabRegister::erase,
             py::arg("name"),
             py::arg("level")
             // "..."
        )
    ;
}
