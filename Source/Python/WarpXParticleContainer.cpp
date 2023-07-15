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
    py::class_<
        WarpXParIter, amrex::ParIter<0,0,PIdx::nattribs>
    >(m, "WarpXParIter")
        .def(py::init<amrex::ParIter<0,0,PIdx::nattribs>::ContainerType&, int>(),
            py::arg("particle_container"), py::arg("level"))
        .def(py::init<amrex::ParIter<0,0,PIdx::nattribs>::ContainerType&, int, amrex::MFItInfo&>(),
            py::arg("particle_container"), py::arg("level"),
            py::arg("info"))
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
        .def("num_real_comps", &WarpXParticleContainer::NumRealComps)
        .def("total_number_of_particles",
            &WarpXParticleContainer::TotalNumberOfParticles,
            py::arg("valid_particles_only"), py::arg("local")
        )
        // .def("deposit_charge",
        //     static_cast<void (WarpXParticleContainer::*)(
        //         WarpXParIter&,
        //         amrex::Gpu::DeviceVector<amrex::Real> const &,
        //         const int * const, amrex::MultiFab*,
        //         const int, const long, const long,
        //         const int, const int, const int
        //     )>(&WarpXParticleContainer::DepositCharge),
        //     py::arg("pti"), py::arg("wp"), py::arg("ion_lev"),
        //     py::arg("rho"), py::arg("icomp"),
        //     py::arg("offset"), py::arg("np_to_depose"),
        //     py::arg("thread_num"), py::arg("lev"), py::arg("depos_lev")
        // )
        .def("deposit_charge",
            [](WarpXParticleContainer& pc, WarpXParIter& pti,
            amrex::Gpu::DeviceVector<amrex::Real> const & wp,
            const int * const ion_lev,
            amrex::MultiFab* rho,
            const int icomp, const long offset, const long np_to_depose,
            const int thread_num, const int lev, const int depos_lev)
            {
                if (*ion_lev == -1)
                    pc.DepositCharge(
                        pti, wp, nullptr, rho, icomp, offset, np_to_depose,
                        thread_num, lev, depos_lev
                    );
                else
                    pc.DepositCharge(
                        pti, wp, ion_lev, rho, icomp, offset, np_to_depose,
                        thread_num, lev, depos_lev
                    );
            }
        )
    ;
}
