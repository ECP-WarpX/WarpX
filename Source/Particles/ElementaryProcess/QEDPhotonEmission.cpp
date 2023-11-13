/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Particles/ElementaryProcess/QEDPhotonEmission.H"

#include "WarpX.H"

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IntVect.H>

#include <array>

PhotonEmissionTransformFunc::
PhotonEmissionTransformFunc (QuantumSynchrotronGetOpticalDepth opt_depth_functor,
                             int const opt_depth_runtime_comp,
                             QuantumSynchrotronPhotonEmission const emission_functor,
                             const WarpXParIter& a_pti, int lev, amrex::IntVect ngEB,
                             amrex::FArrayBox const& exfab,
                             amrex::FArrayBox const& eyfab,
                             amrex::FArrayBox const& ezfab,
                             amrex::FArrayBox const& bxfab,
                             amrex::FArrayBox const& byfab,
                             amrex::FArrayBox const& bzfab,
                             amrex::Vector<amrex::ParticleReal>& E_external_particle,
                             amrex::Vector<amrex::ParticleReal>& B_external_particle,
                             int a_offset):
    m_opt_depth_functor{opt_depth_functor},
    m_opt_depth_runtime_comp{opt_depth_runtime_comp},
    m_emission_functor{emission_functor},
    m_galerkin_interpolation{WarpX::galerkin_interpolation},
    m_nox{WarpX::nox},
    m_n_rz_azimuthal_modes{WarpX::n_rz_azimuthal_modes}
{

    using namespace amrex::literals;

    m_get_position  = GetParticlePosition<PIdx>(a_pti, a_offset);
    m_get_externalEB = GetExternalEBField(a_pti, a_offset);
    m_Ex_external_particle = E_external_particle[0];
    m_Ey_external_particle = E_external_particle[1];
    m_Ez_external_particle = E_external_particle[2];
    m_Bx_external_particle = B_external_particle[0];
    m_By_external_particle = B_external_particle[1];
    m_Bz_external_particle = B_external_particle[2];

    m_ex_arr = exfab.array();
    m_ey_arr = eyfab.array();
    m_ez_arr = ezfab.array();
    m_bx_arr = bxfab.array();
    m_by_arr = byfab.array();
    m_bz_arr = bzfab.array();

    m_ex_type = exfab.box().ixType();
    m_ey_type = eyfab.box().ixType();
    m_ez_type = ezfab.box().ixType();
    m_bx_type = bxfab.box().ixType();
    m_by_type = byfab.box().ixType();
    m_bz_type = bzfab.box().ixType();

    amrex::Box box = a_pti.tilebox();
    box.grow(ngEB);

    const std::array<amrex::Real,3>& dx = WarpX::CellSize(std::max(lev, 0));
    m_dx_arr = {dx[0], dx[1], dx[2]};

    // Lower corner of tile box physical domain (take into account Galilean shift)
    const std::array<amrex::Real, 3>& xyzmin = WarpX::LowerCorner(box, lev, 0._rt);
    m_xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};



    m_lo = amrex::lbound(box);
}
