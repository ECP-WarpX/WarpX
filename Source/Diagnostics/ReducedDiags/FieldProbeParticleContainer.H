/* Copyright 2021 Tiberius Rheaume, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FieldProbeParticleContainer_H_
#define WARPX_FieldProbeParticleContainer_H_

#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include <AMReX_BaseFwd.H>
#include <AMReX_AmrCoreFwd.H>

/**
 * This enumerated struct is used to index the field probe particle
 * values that are being stored as SoA data. Nattribs
 * is enumerated to give the number of attributes stored.
 */
struct ParticleVal
{
    enum
    {
        Ex = 0, Ey, Ez,
        Bx, By, Bz,
        S,
        nattribs
    };
};

/**
 * This class defines the FieldProbeParticleContainer
 * which is branched from the amrex::ParticleContainer.
 * nattribs tells the particle container to allot 7 SOA values.
 */
class FieldProbeParticleContainer
    : public amrex::ParticleContainer<0,0,ParticleVal::nattribs>
{
public:
    FieldProbeParticleContainer (amrex::AmrCore* amr_core);
    virtual ~FieldProbeParticleContainer() {}

    //! amrex iterator for our number of attributes
    using iterator = amrex::ParIter<0,0,ParticleVal::nattribs,0>;

    //similar to WarpXParticleContainer::AddNParticles but does not include v(x,y,z)
    void AddNParticles (int lev, int n, const amrex::ParticleReal* x, const amrex::ParticleReal* y, const amrex::ParticleReal* z);
};

#endif // WARPX_FieldProbeParticleContainer_H_
