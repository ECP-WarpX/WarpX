/* Copyright 2021-2023 The WarpX Community
 *
 * Authors: Axel Huebl, Remi Lehe, Roelof Groenewald
 * License: BSD-3-Clause-LBNL
 */

#include <Particles/ParticleBoundaryBuffer.H>
#include <Particles/MultiParticleContainer.H>
#include <Python/pyWarpX.H>

namespace py = pybind11;


void init_ParticleBoundaryBuffer (py::module& m)
{
    py::class_<ParticleBoundaryBuffer>(m, "ParticleBoundaryBuffer")
        .def(py::init<>())
        .def("get_num_particles_in_container",
            [](ParticleBoundaryBuffer& pbc,
            const std::string species_name, int boundary, bool local)
            {
                return pbc.getNumParticlesInContainer(species_name, boundary, local);
            },
            py::arg("species_name"), py::arg("boundary"), py::arg("local")=true
        )
    ;
}
