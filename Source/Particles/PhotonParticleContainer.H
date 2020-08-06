/* Copyright 2019 Andrew Myers, David Grote, Luca Fedeli
 * Maxence Thevenet, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PhotonParticleContainer_H_
#define WARPX_PhotonParticleContainer_H_

#include "PhysicalParticleContainer.H"
#include <AMReX_Vector.H>

/**
 * Photon particles have no mass, they deposit no charge, and see specific QED
 * effects. For these reasons, they are stored in the separate particle
 * container PhotonParticleContainer, that inherts from
 * PhysicalParticleContainer. The particle pusher and current deposition, in
 * particular, are overriden in this container.
 */
class PhotonParticleContainer
    : public PhysicalParticleContainer
{
public:
    PhotonParticleContainer (amrex::AmrCore* amr_core,
                                    int ispecies,
                                    const std::string& name);
    virtual ~PhotonParticleContainer () {}

    virtual void InitData() override;

    virtual void Evolve (int lev,
                         const amrex::MultiFab& Ex,
                         const amrex::MultiFab& Ey,
                         const amrex::MultiFab& Ez,
                         const amrex::MultiFab& Bx,
                         const amrex::MultiFab& By,
                         const amrex::MultiFab& Bz,
                         const amrex::MultiFab& Ex_avg,
                         const amrex::MultiFab& Ey_avg,
                         const amrex::MultiFab& Ez_avg,
                         const amrex::MultiFab& Bx_avg,
                         const amrex::MultiFab& By_avg,
                         const amrex::MultiFab& Bz_avg,
                         amrex::MultiFab& jx,
                         amrex::MultiFab& jy,
                         amrex::MultiFab& jz,
                         amrex::MultiFab* cjx,
                         amrex::MultiFab* cjy,
                         amrex::MultiFab* cjz,
                         amrex::MultiFab* rho,
                         amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx,
                         const amrex::MultiFab* cEy,
                         const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx,
                         const amrex::MultiFab* cBy,
                         const amrex::MultiFab* cBz,
                         amrex::Real t,
                         amrex::Real dt,
                         DtType a_dt_type=DtType::Full) override;

    virtual void PushPX(WarpXParIter& pti,
                        amrex::FArrayBox const * exfab,
                        amrex::FArrayBox const * eyfab,
                        amrex::FArrayBox const * ezfab,
                        amrex::FArrayBox const * bxfab,
                        amrex::FArrayBox const * byfab,
                        amrex::FArrayBox const * bzfab,
                        const int ngE, const int /*e_is_nodal*/,
                        const long offset,
                        const long np_to_push,
                        int lev, int gather_lev,
                        amrex::Real dt, ScaleFields scaleFields,
                        DtType a_dt_type) override;

    // Do nothing
    virtual void PushP (int /*lev*/,
                        amrex::Real /*dt*/,
                        const amrex::MultiFab& /*Ex*/,
                        const amrex::MultiFab& /*Ey*/,
                        const amrex::MultiFab& /*Ez*/,
                        const amrex::MultiFab& /*Bx*/,
                        const amrex::MultiFab& /*By*/,
                        const amrex::MultiFab& /*Bz*/) override {}


    // DepositCurrent should do nothing for photons
    virtual void DepositCurrent(WarpXParIter& /*pti*/,
                                RealVector& /*wp*/,
                                RealVector& /*uxp*/,
                                RealVector& /*uyp*/,
                                RealVector& /*uzp*/,
                                const int * const /*ion_lev*/,
                                amrex::MultiFab* /*jx*/,
                                amrex::MultiFab* /*jy*/,
                                amrex::MultiFab* /*jz*/,
                                const long /*offset*/,
                                const long /*np_to_depose*/,
                                int /*thread_num*/,
                                int /*lev*/,
                                int /*depos_lev*/,
                                amrex::Real /*dt*/) override {}
};

#endif // #ifndef WARPX_PhotonParticleContainer_H_
