/* Copyright 2024 The WarpX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */
#include "Python/pyWarpX.H"

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>

void init_MultiFabRegister (py::module & m)
{
    using ablastr::fields::MultiFabRegister;

    py::class_<ablastr::fields::MultiFabRegister>(m, "MultiFabRegister")

        .def("alloc_init",
             py::overload_cast<
                 std::string,
                 int,
                 amrex::BoxArray const &,
                 amrex::DistributionMapping const &,
                 int,
                 amrex::IntVect const &,
                 std::optional<const amrex::Real>,
                 bool,
                 bool
             >(&MultiFabRegister::alloc_init),
             py::arg("name"),
             py::arg("level"),
             py::arg("ba"),
             py::arg("dm"),
             py::arg("ncomp"),
             py::arg("ngrow"),
             py::arg("initial_value"),
             py::arg("redistribute"),
             py::arg("redistribute_on_remake")
        )

        .def("alloc_init",
             py::overload_cast<
                 std::string,
                 ablastr::fields::Direction,
                 int,
                 amrex::BoxArray const &,
                 amrex::DistributionMapping const &,
                 int,
                 amrex::IntVect const &,
                 std::optional<const amrex::Real>,
                 bool,
                 bool
             >(&MultiFabRegister::alloc_init),
             py::arg("name"),
             py::arg("dir"),
             py::arg("level"),
             py::arg("ba"),
             py::arg("dm"),
             py::arg("ncomp"),
             py::arg("ngrow"),
             py::arg("initial_value"),
             py::arg("redistribute"),
             py::arg("redistribute_on_remake")
        )

        .def("alias_init",
             py::overload_cast<
                 std::string,
                 std::string,
                 int,
                 std::optional<const amrex::Real>
             >(&MultiFabRegister::alias_init)
             // ... py::arg("name")
             // "..."
        )

        .def("alias_init",
             py::overload_cast<
                 std::string,
                 std::string,
                 ablastr::fields::Direction,
                 int,
                 std::optional<const amrex::Real>
             >(&MultiFabRegister::alias_init)
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
             py::overload_cast<
                 std::string,
                 int
             >(&MultiFabRegister::has)
             //py::arg("name"),
             // "..."
        )

        .def("has",
             py::overload_cast<
                 std::string,
                 ablastr::fields::Direction,
                 int
             >(&MultiFabRegister::has)
             //py::arg("name"),
             // "..."
        )

        .def("get",
             py::overload_cast<
                 std::string,
                 int
             >(&MultiFabRegister::get)
             //py::arg("name"),
             // "..."
        )

        .def("get",
             py::overload_cast<
                 std::string,
                 ablastr::fields::Direction,
                 int
             >(&MultiFabRegister::get)
             //py::arg("name"),
             // "..."
        )

        //.def("list",
        //     &MultiFabRegister::list
        //     // "..."
        //)

        .def("erase",
             py::overload_cast<
                 std::string,
                 int
             >(&MultiFabRegister::erase)
             //py::arg("name"),
             // "..."
        )

        .def("erase",
             py::overload_cast<
                 std::string,
                 ablastr::fields::Direction,
                 int
             >(&MultiFabRegister::erase)
             //py::arg("name"),
             // "..."
        )
    ;
}
