/* Copyright 2021 Elisa Rheaume, Axel Huebl
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
#include <AMReX_Vector.H>

/**
 * This enumerated struct is used to index the field probe particle
 * values that are being stored as SoA data. Nattribs
 * is enumerated to give the number of attributes stored.
 */
struct FieldProbePIdx
{
    enum
    {
        Ex = 0, Ey, Ez,
        Bx, By, Bz,
        S, //!< the Poynting vector
        nattribs
    };
};

/**
 * This class defines the FieldProbeParticleContainer
 * which is branched from the amrex::ParticleContainer.
 * nattribs tells the particle container to allot 7 SOA values.
 */
class FieldProbeParticleContainer
    : public amrex::ParticleContainer<0, 0, FieldProbePIdx::nattribs>
{
public:
    FieldProbeParticleContainer (amrex::AmrCore* amr_core);
    virtual ~FieldProbeParticleContainer() {}

    //! amrex iterator for our number of attributes
    using iterator = amrex::ParIter<0, 0, FieldProbePIdx::nattribs, 0>;
    //! amrex iterator for our number of attributes (read-only)
    using const_iterator = amrex::ParConstIter<0, 0, FieldProbePIdx::nattribs, 0>;

    //! similar to WarpXParticleContainer::AddNParticles but does not include u(x,y,z)
    void AddNParticles (int lev, amrex::Vector<amrex::ParticleReal> const & x, amrex::Vector<amrex::ParticleReal> const & y, amrex::Vector<amrex::ParticleReal> const & z);
};

#endif // WARPX_FieldProbeParticleContainer_H_
