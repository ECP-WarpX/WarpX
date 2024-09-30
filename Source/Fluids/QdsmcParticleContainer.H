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
        x, vx,
#endif
#if defined (WARPX_DIM_3D)
        y, vy,
#endif
        z, vz,
        KeNe,
#ifdef WARPX_DIM_RZ
        theta, vtheta      ///< RZ needs all three position components
#endif
        nattribs //shoud be = 2*dim + number of real attributes (KeNe for now)
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
    void AddNParticles (int lev, 
                        amrex::Vector<amrex::ParticleReal> const & x, 
                        amrex::Vector<amrex::ParticleReal> const & y, 
                        amrex::Vector<amrex::ParticleReal> const & z,
                        amrex::Vector<amrex::ParticleReal> const & vx, 
                        amrex::Vector<amrex::ParticleReal> const & vy, 
                        amrex::Vector<amrex::ParticleReal> const & vz,
                        amrex::Vector<amrex::ParticleReal> const & kN);

    // Function that Initializes the qdsmc particles positions
    void InitX();

    // Function that Initializes the qdsmc particles velocities
    void InitV();

    // Function that pushes particle positions by one time step
    void PushX(amrex::Real dt);

    // Function that gathers a scalar multifab value to particles positions
    // No interpolation needed since particles are initialized on each step at node positions
    void Gather();

    // Function that deposits a scalar quantity from particles to a multifab.
    // Takes vector with a specific particle value (can be anything) and scalar multifab to deposit in
    void Deposit(amrex::Vector<amrex::ParticleReal> const & p_val, ablastr::fields::ScalarField& scalar_field);




};

#endif // WARPX_QdsmcParticleContainer_H_
