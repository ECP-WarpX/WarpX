/* Copyright 2021 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEIO_H_
#define PARTICLEIO_H_

#include "Particles/WarpXParticleContainer.H"

#include <AMReX_AmrParticles.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

enum struct ConvertDirection{WarpX_to_SI, SI_to_WarpX};

/** Convert particle momentum to/from SI
 *
 * Particle momentum is defined as gamma*velocity, which is neither
 * SI mass*gamma*velocity nor normalized gamma*velocity/c.
 * This converts momentum to SI units (or vice-versa) to write SI data
 * to file.
 * Photons are a special case, since particle momentum is defined as
 * (photon_energy/(m_e * c) ) * u, where u is the photon direction (a
 * unit vector).
 *
 * @tparam T_ParticleContainer a WarpX particle container or AmrParticleContainer
 * @param convert_direction convert to or from SI
 * @param pc the particle container to manipulate
 * @param mass the particle rest mass to use for conversion
 */
template< typename T_ParticleContainer >
void
particlesConvertUnits (ConvertDirection convert_direction, T_ParticleContainer* pc, amrex::ParticleReal const mass )
{
    using namespace amrex;

    // Compute conversion factor
    auto factor = 1_rt;

    if (convert_direction == ConvertDirection::WarpX_to_SI){
        factor = mass;
    } else if (convert_direction == ConvertDirection::SI_to_WarpX){
        factor = 1._rt/mass;
    }

    const int nLevels = pc->finestLevel();
    for (int lev=0; lev<=nLevels; lev++){
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIter pti(*pc, lev); pti.isValid(); ++pti)
        {
            // - momenta are stored as a struct of array, in `attribs`
            // The GetStructOfArrays is called directly since the convenience routine GetAttribs
            // is only available in WarpXParIter. ParIter is used here since the pc passed in
            // will sometimes be a PinnedMemoryParticleContainer (not derived from a WarpXParticleContainer).
            auto& attribs = pti.GetStructOfArrays().GetRealData();
            ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
            // Loop over the particles and convert momentum
            const long np = pti.numParticles();
            ParallelFor( np,
                         [=] AMREX_GPU_DEVICE (long i) {
                ux[i] *= factor;
                uy[i] *= factor;
                uz[i] *= factor;
            }
            );
        }
    }
}

#endif /* PARTICLEIO_H_ */
