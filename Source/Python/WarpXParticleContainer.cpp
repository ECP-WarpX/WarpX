/* Copyright 2021-2022 The WarpX Community
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "pyWarpX.H"

#include <Particles/WarpXParticleContainer.H>

namespace py = pybind11;

void init_WarpXParIter (py::module& m)
{
    py::class_<WarpXParIter, amrex::ParIter<0,0,PIdx::nattribs>
    > wpi (m, "WarpXParIter");
    wpi
        .def(py::init<WarpXParticleContainer&, int>())
        ;
}

void init_WarpXParticleContainer (py::module& m)
{
    py::class_<
        WarpXParticleContainer,
        amrex::ParticleContainer<0, 0, PIdx::nattribs, 0>
    > wpc (m, "WarpXParticleContainer");
    wpc
        .def("add_real_comp",
            [](WarpXParticleContainer& pc, const std::string& name, bool const comm) { pc.AddRealComp(name, comm); },
            py::arg("name"), py::arg("comm")
        )
        .def("add_n_particles",
            &WarpXParticleContainer::AddNParticles,
            py::arg("lev"), py::arg("n"),
            py::arg("x"), py::arg("y"), py::arg("z"),
            py::arg("ux"), py::arg("uy"), py::arg("uz"),
            py::arg("nattr_real"), py::arg("attr_real"),
            py::arg("nattr_int"), py::arg("attr_int"),
            py::arg("uniqueparticles"), py::arg("id")
        )
        .def("deposit_charge",
            static_cast<void (WarpXParticleContainer::*)(
                WarpXParIter&,
                amrex::Gpu::DeviceVector<amrex::Real> const &,
                const int * const, amrex::MultiFab*,
                const int, const long, const long,
                const int, const int, const int
            )>(&WarpXParticleContainer::DepositCharge),
            py::arg("pti"), py::arg("wp"), py::arg("ion_lev"),
            py::arg("rho"), py::arg("icomp"),
            py::arg("offset"), py::arg("np_to_depose"),
            py::arg("thread_num"), py::arg("lev"), py::arg("depos_lev")
        )
        ;
}
