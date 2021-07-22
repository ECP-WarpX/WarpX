/* Copyright 2020 Luca Fedeli, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_schwinger_process_wrapper_h_
#define WARPX_schwinger_process_wrapper_h_

#include "QedWrapperCommons.H"

#include <AMReX_Random.H>

#include <picsar_qed/physics/schwinger/schwinger_pair_engine_core.hpp>

#include <cmath>

/**
 * This function is a wrapper around a PICSAR function which calculates
 * the number of Schwinger pairs created at a given timestep and in a given
 * cell as a function of the EM field in that cell.
 *
 * @param[in] dV Volume of the cell.
 * @param[in] dt temporal step.
 * @param[in] Ex x-component of the electric field on the cell.
 * @param[in] Ey y-component of the electric field on the cell.
 * @param[in] Ez z-component of the electric field on the cell.
 * @param[in] Bx x-component of the magnetic field on the cell.
 * @param[in] By y-component of the magnetic field on the cell.
 * @param[in] Bz z-component of the magnetic field on the cell.
 * @param[in] PoissonToGaussianThreshold If the expected number of created
 *            pairs is below this parameter, a Poisson distribution is used
 *            to draw the number of created pairs. Otherwise a Gaussian
 *            distribution is used.
 * @return the number of pairs generated via the Schwinger process
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
getSchwingerProductionNumber (const amrex::Real dV, const amrex::Real dt,
              const amrex::ParticleReal Ex, const amrex::ParticleReal Ey, const amrex::ParticleReal Ez,
              const amrex::ParticleReal Bx, const amrex::ParticleReal By, const amrex::ParticleReal Bz,
              const amrex::ParticleReal PoissonToGaussianThreshold,
              amrex::RandomEngine const& engine)
{
    using namespace amrex;
    namespace pxr_p = picsar::multi_physics::phys;
    namespace pxr_sh = picsar::multi_physics::phys::schwinger;

    const auto expectedPairNumber =
        pxr_sh::expected_pair_number<amrex::Real, pxr_p::unit_system::SI>(
            Ex, Ey, Ez, Bx, By, Bz, dV, dt);

    if (expectedPairNumber <= PoissonToGaussianThreshold) {
        return amrex::RandomPoisson(expectedPairNumber, engine);
    } else {
        const auto numpairs =
            amrex::RandomNormal(expectedPairNumber, std::sqrt(expectedPairNumber), engine);
        return numpairs > 0._rt ? numpairs : 0._rt;
    }
}

#endif // WARPX_schwinger_process_wrapper_h_
