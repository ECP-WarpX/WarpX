/* Copyright 2021-2022 The WarpX Community
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */

#include <Particles/WarpXParticleContainer.H>
#include <Python/pyWarpX.H>

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
        .def("num_local_tiles_at_level",
            &WarpXParticleContainer::numLocalTilesAtLevel,
            py::arg("level")
        )
        .def("total_number_of_particles",
            &WarpXParticleContainer::TotalNumberOfParticles,
            py::arg("valid_particles_only"), py::arg("local")
        )
        .def("deposit_charge",
            [](WarpXParticleContainer& pc,
            amrex::MultiFab* rho, const int lev)
            {
                for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti)
                {
                    const long np = pti.numParticles();
                    auto& wp = pti.GetAttribs(PIdx::w);
                    pc.DepositCharge(pti, wp, nullptr, rho, 0, 0, np, 0, lev, lev);
                }
            },
            py::arg("rho"), py::arg("lev")
        )
    ;
}
