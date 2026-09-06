/* Copyright 2025 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_QED_EVOLVE_OPTICAL_DEPTH_H_
#define WARPX_QED_EVOLVE_OPTICAL_DEPTH_H_

#include "Particles/Gather/FieldGather.H"
#include "Particles/Gather/GetExternalFields.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "QEDInternals/QuantumSyncEngineWrapper.H"

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * \brief Functor advancing the quantum synchrotron optical depth of one particle.
 *
 * The optical depth is a Monte-Carlo quantity, decremented at the photon emission rate
 * dN/dt(energy, chi), and triggers an emission when it becomes negative. It is advanced
 * in its own pass rather than inside the particle pusher, so that the momentum push carries
 * no compile-time QED option and the implicit solver never has to roll it back.
 *
 * The fields are gathered at the particle position, as PhotonEmissionTransformFunc does.
 */
class QuantumSyncEvolveOpticalDepthFunc
{
public:
    /**
    * \brief Constructor of the QuantumSyncEvolveOpticalDepthFunc functor.
    *
    * @param[in] evolve_functor functor advancing the optical depth of a single particle
    * @param[in] ux,uy,uz pointers to the momentum components of the source species
    * @param[in] opt_depth pointer to the optical depth component of the source species
    * @param[in] dt the time step over which the optical depth is advanced
    * @param[in] a_pti particle iterator of the species undergoing quantum synchrotron
    * @param[in] lev the mesh-refinement level
    * @param[in] ngEB number of guard cells allocated for the E and B MultiFabs
    * @param[in] exfab,eyfab,ezfab FArrayBoxes of the electric field components
    * @param[in] bxfab,byfab,bzfab FArrayBoxes of the magnetic field components
    * @param[in] E_external_particle reference to the external E field
    * @param[in] B_external_particle reference to the external B field
    * @param[in] a_offset offset to apply to the particle indices
    */
    QuantumSyncEvolveOpticalDepthFunc (
        QuantumSynchrotronEvolveOpticalDepth evolve_functor,
        amrex::ParticleReal* AMREX_RESTRICT ux,
        amrex::ParticleReal* AMREX_RESTRICT uy,
        amrex::ParticleReal* AMREX_RESTRICT uz,
        amrex::ParticleReal* AMREX_RESTRICT opt_depth,
        amrex::Real dt,
        const WarpXParIter& a_pti, int lev, amrex::IntVect ngEB,
        amrex::FArrayBox const& exfab,
        amrex::FArrayBox const& eyfab,
        amrex::FArrayBox const& ezfab,
        amrex::FArrayBox const& bxfab,
        amrex::FArrayBox const& byfab,
        amrex::FArrayBox const& bzfab,
        amrex::Vector<amrex::ParticleReal>& E_external_particle,
        amrex::Vector<amrex::ParticleReal>& B_external_particle,
        int a_offset = 0);

    /**
    * \brief Functor call. Gathers the fields at the particle position and advances
    * its optical depth by one time step.
    *
    * @param[in] ip index of the particle
    */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void operator() (long ip) const noexcept
    {
        amrex::ParticleReal xp, yp, zp;
        m_get_position(ip, xp, yp, zp);

        amrex::ParticleReal ex = m_Ex_external_particle;
        amrex::ParticleReal ey = m_Ey_external_particle;
        amrex::ParticleReal ez = m_Ez_external_particle;
        amrex::ParticleReal bx = m_Bx_external_particle;
        amrex::ParticleReal by = m_By_external_particle;
        amrex::ParticleReal bz = m_Bz_external_particle;

        m_get_externalEB(ip, ex, ey, ez, bx, by, bz);

        doGatherShapeN(xp, yp, zp, ex, ey, ez, bx, by, bz,
                       m_ex_arr, m_ey_arr, m_ez_arr, m_bx_arr, m_by_arr, m_bz_arr,
                       m_ex_type, m_ey_type, m_ez_type, m_bx_type, m_by_type, m_bz_type,
                       m_dinv, m_xyzmin, m_lo, m_n_rz_azimuthal_modes,
                       m_nox, m_galerkin_interpolation);

        m_evolve_functor(m_ux[ip], m_uy[ip], m_uz[ip],
                         ex, ey, ez, bx, by, bz,
                         m_dt, m_opt_depth[ip]);
    }

private:
    QuantumSynchrotronEvolveOpticalDepth
        m_evolve_functor; /*!< A copy of the functor advancing the optical depth. It holds only
                               a few parameters and pointers to the lookup tables. */

    amrex::ParticleReal* AMREX_RESTRICT m_ux;
    amrex::ParticleReal* AMREX_RESTRICT m_uy;
    amrex::ParticleReal* AMREX_RESTRICT m_uz;
    amrex::ParticleReal* AMREX_RESTRICT m_opt_depth;
    amrex::Real m_dt;

    GetParticlePosition<PIdx> m_get_position;
    GetExternalEBField m_get_externalEB;
    amrex::ParticleReal m_Ex_external_particle;
    amrex::ParticleReal m_Ey_external_particle;
    amrex::ParticleReal m_Ez_external_particle;
    amrex::ParticleReal m_Bx_external_particle;
    amrex::ParticleReal m_By_external_particle;
    amrex::ParticleReal m_Bz_external_particle;

    amrex::Array4<const amrex::Real> m_ex_arr;
    amrex::Array4<const amrex::Real> m_ey_arr;
    amrex::Array4<const amrex::Real> m_ez_arr;
    amrex::Array4<const amrex::Real> m_bx_arr;
    amrex::Array4<const amrex::Real> m_by_arr;
    amrex::Array4<const amrex::Real> m_bz_arr;

    amrex::IndexType m_ex_type;
    amrex::IndexType m_ey_type;
    amrex::IndexType m_ez_type;
    amrex::IndexType m_bx_type;
    amrex::IndexType m_by_type;
    amrex::IndexType m_bz_type;

    amrex::XDim3 m_dinv;
    amrex::XDim3 m_xyzmin;

    bool m_galerkin_interpolation;
    int m_nox;
    int m_n_rz_azimuthal_modes;

    amrex::Dim3 m_lo;
};

#endif // WARPX_QED_EVOLVE_OPTICAL_DEPTH_H_
