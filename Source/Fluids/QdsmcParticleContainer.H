/* Copyright 2024 Marco Acciarri (Helion Energy Inc.)
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_QdsmcParticleContainer_H_
#define WARPX_QdsmcParticleContainer_H_

#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include <AMReX_BaseFwd.H>
#include <AMReX_AmrCoreFwd.H>
#include <AMReX_Vector.H>

#include "Particles/ParticleBoundaries.H"
#include <ablastr/fields/MultiFabRegister.H>

/**
 * This enumerated struct is used to index the field probe particle
 * values that are being stored as SoA data. Nattribs
 * is enumerated to give the number of attributes stored.
 */
struct QdsmcPart
{
    enum
    {
#if !defined (WARPX_DIM_1D_Z)
        x,
#endif
#if defined (WARPX_DIM_3D)
        y,
#endif
        z,
        vx, vy, vz,
        entropy, np_real,

#ifdef WARPX_DIM_RZ
        theta,      ///< RZ needs all three position components
#endif
        nattribs //shoud be = ??
    };
};

/**
 * This class defines the QdsmcParticleContainer
 * which is branched from the amrex::ParticleContainer.
 */
class QdsmcParticleContainer
    : public amrex::ParticleContainerPureSoA<QdsmcPart::nattribs, 0>
{
public:
    static constexpr int NStructReal = 0;
    static constexpr int NStructInt = 0;
    static constexpr int NReal = QdsmcPart::nattribs;
    static constexpr int NInt = 0;

    QdsmcParticleContainer (amrex::AmrCore* amr_core);
    ~QdsmcParticleContainer() override = default;

    QdsmcParticleContainer ( QdsmcParticleContainer const &)             = delete;
    QdsmcParticleContainer& operator= ( QdsmcParticleContainer const & ) = delete;
    QdsmcParticleContainer ( QdsmcParticleContainer&& )                  = default;
    QdsmcParticleContainer& operator= ( QdsmcParticleContainer&& )       = default;

    //! amrex iterator for our number of attributes
    using iterator = amrex::ParIterSoA<QdsmcPart::nattribs, 0>;
    //! amrex iterator for our number of attributes (read-only)
    using const_iterator = amrex::ParConstIterSoA<QdsmcPart::nattribs, 0>;

    //! similar to WarpXParticleContainer::AddNParticles
    void AddNParticles (int lev, amrex::Long n,
                        amrex::Vector<amrex::ParticleReal> const & x, 
                        amrex::Vector<amrex::ParticleReal> const & y, 
                        amrex::Vector<amrex::ParticleReal> const & z,
                        amrex::Vector<amrex::ParticleReal> const & vx, 
                        amrex::Vector<amrex::ParticleReal> const & vy, 
                        amrex::Vector<amrex::ParticleReal> const & vz,
                        amrex::Vector<amrex::ParticleReal> const & entropy,
                        amrex::Vector<amrex::ParticleReal> const & np_real)
                        //amrex::Long id);

    // Function that Initializes the qdsmc particles positions
    void InitX();

    /* Function that Initializes the qdsmc particles velocities
    *  without including a thermal spread, this is only used by 
    *  the qdsmc solver for the electron temperature used in the 
    *  hybrid PIC model
    */
    void InitV_drift(ablastr::fields::VectorField& U);

    // Function that pushes particle positions by one time step
    void PushX(amrex::Real dt);

    /* Function that gathers a scalar multifab value to particles positions
    *  No interpolation needed since particles are initialized on each step at node positions
    */
    void Gather();

    // Function that deposits a scalar quantity from particles to a multifab.
    // Takes vector with a specific particle value (can be anything) and scalar multifab to deposit in
    //void DepositK(amrex::Vector<amrex::ParticleReal> const & p_val, ablastr::fields::ScalarField& scalar_field);

};

#endif // WARPX_QdsmcParticleContainer_H_
