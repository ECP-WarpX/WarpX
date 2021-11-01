/* Copyright 2021 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMPUTEPARTICLEDIAGFUNCTOR_H_
#define WARPX_COMPUTEPARTICLEDIAGFUNCTOR_H_

#include "Particles/WarpXParticleContainer.H"
#include <AMReX.H>
#include <AMReX_AmrParticles.H>

/**
 * \brief Functor to compute a diagnostic and store the result in existing ParticleContainer.
 */
class
ComputeParticleDiagFunctor
{
public:
    using ParticleContainer = typename amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;

    ComputeParticleDiagFunctor( ) {}
    /** Virtual Destructor to handle clean destruction of derived classes */
    virtual ~ComputeParticleDiagFunctor() = default;

    virtual void PrepareFunctorData ( int i_buffer, bool ZSliceInDomain,
                                      amrex::Real old_z_boost,
                                      amrex::Real current_z_boost, amrex::Real t_lab)
                                      {
                                          amrex::ignore_unused(i_buffer,
                                          ZSliceInDomain, old_z_boost,
                                          current_z_boost, t_lab);
                                      }
    /** Compute particle attributes and store the result in pc_dst particle container.
     *  \param[out] pc_dst  output particle container where the result is stored.
     *  \param[in] i_buffer snapshot index for which the particle buffer is processed
     */
    virtual void operator () (ParticleContainer& pc_dst, int &totalParticlesInBuffer, int i_buffer = 0) const = 0;
    virtual void InitData() {}
};

#endif // WARPX_COMPUTEPARTICLEDIAGFUNCTOR_H_
