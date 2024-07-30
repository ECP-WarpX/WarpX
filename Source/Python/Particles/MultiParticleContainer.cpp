/* Copyright 2021-2022 The WarpX Community
 *
 * Authors: Axel Huebl, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */

#include "Python/pyWarpX.H"

#include <Particles/MultiParticleContainer.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>


void init_MultiParticleContainer (py::module& m)
{
    py::class_<MultiParticleContainer>(m, "MultiParticleContainer")
        .def("get_particle_container_from_name",
            &MultiParticleContainer::GetParticleContainerFromName,
            py::arg("name"),
            py::return_value_policy::reference_internal
        )

        .def("set_plasma_lens_strength",
             [](MultiParticleContainer& mpc,
             int i_lens,
             amrex::Real strength_Ex,
             amrex::Real strength_Ey,
             amrex::Real strength_Bx,
             amrex::Real strength_By
             ) {
                 mpc.h_repeated_plasma_lens_strengths_Ex.at(i_lens) = strength_Ex;
                 mpc.h_repeated_plasma_lens_strengths_Ey.at(i_lens) = strength_Ey;
                 mpc.h_repeated_plasma_lens_strengths_Bx.at(i_lens) = strength_Bx;
                 mpc.h_repeated_plasma_lens_strengths_By.at(i_lens) = strength_By;
                 amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                       mpc.h_repeated_plasma_lens_strengths_Ex.begin(), mpc.h_repeated_plasma_lens_strengths_Ex.end(),
                                       mpc.d_repeated_plasma_lens_strengths_Ex.begin());
                 amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                       mpc.h_repeated_plasma_lens_strengths_Ey.begin(), mpc.h_repeated_plasma_lens_strengths_Ey.end(),
                                       mpc.d_repeated_plasma_lens_strengths_Ey.begin());
                 amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                       mpc.h_repeated_plasma_lens_strengths_Bx.begin(), mpc.h_repeated_plasma_lens_strengths_Bx.end(),
                                       mpc.d_repeated_plasma_lens_strengths_Bx.begin());
                 amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                       mpc.h_repeated_plasma_lens_strengths_By.begin(), mpc.h_repeated_plasma_lens_strengths_By.end(),
                                       mpc.d_repeated_plasma_lens_strengths_By.begin());
                 amrex::Gpu::synchronize();
             },
             py::arg("i_lens"), py::arg("strength_Ex"), py::arg("strength_Ey"), py::arg("strength_Bx"), py::arg("strength_By"),
             R"pbdoc(Set the strength of the `i_lens`-th lens
Parameters
----------
i_lens: int
  Index of the lens to be modified
strength_Ex, strength_Ey, strength_Bx, strength_By: floats
  The electric and magnetic x- and y- focusing strengths of the lens)pbdoc"
        )
    ;
}
