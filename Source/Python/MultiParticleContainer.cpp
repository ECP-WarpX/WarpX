/* Copyright 2021-2022 The WarpX Community
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include "pyWarpX.H"

#include <Particles/MultiParticleContainer.H>

namespace py = pybind11;

void init_MultiParticleContainer (py::module& m)
{
    py::class_<MultiParticleContainer> mpc(m, "MultiParticleContainer");
}
