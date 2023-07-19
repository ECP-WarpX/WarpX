/* Copyright 2021-2022 The WarpX Community
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */

#include <Particles/MultiParticleContainer.H>
#include <Particles/WarpXParticleContainer.H>
#include <Python/pyWarpX.H>

namespace py = pybind11;

void init_MultiParticleContainer (py::module& m)
{
    py::class_<MultiParticleContainer> mpc(m, "MultiParticleContainer");
    mpc
        .def("get_particle_container_from_name",
            &MultiParticleContainer::GetParticleContainerFromName,
            py::arg("name"),
            py::return_value_policy::reference_internal
        )
    ;
}
