/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Particles/ElementaryProcess/Ionization.H"

#include "WarpX.H"

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IntVect.H>

#include <algorithm>
#include <array>

IonizationFilterFunc::IonizationFilterFunc (const WarpXParIter& a_pti, int lev, amrex::IntVect ngE,
                                            amrex::FArrayBox const& exfab,
                                            amrex::FArrayBox const& eyfab,
                                            amrex::FArrayBox const& ezfab,
                                            amrex::FArrayBox const& bxfab,
                                            amrex::FArrayBox const& byfab,
                                            amrex::FArrayBox const& bzfab,
                                            amrex::Vector<amrex::Real> v_galilean,
                                            const amrex::Real* const AMREX_RESTRICT a_ionization_energies,
                                            const amrex::Real* const AMREX_RESTRICT a_adk_prefactor,
                                            const amrex::Real* const AMREX_RESTRICT a_adk_exp_prefactor,
                                            const amrex::Real* const AMREX_RESTRICT a_adk_power,
                                            int a_comp,
                                            int a_atomic_number,
                                            int a_offset) noexcept
{
    m_ionization_energies = a_ionization_energies;
    m_adk_prefactor = a_adk_prefactor;
    m_adk_exp_prefactor = a_adk_exp_prefactor;
    m_adk_power = a_adk_power;
    comp = a_comp;
    m_atomic_number = a_atomic_number;

    m_get_position  = GetParticlePosition(a_pti, a_offset);
    m_get_externalEB = GetExternalEBField(a_pti, a_offset);

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
    box.grow(ngE);

    const std::array<amrex::Real,3>& dx = WarpX::CellSize(std::max(lev, 0));
    m_dx_arr = {dx[0], dx[1], dx[2]};

    // Lower corner of tile box physical domain (take into account Galilean shift)
    amrex::Real cur_time = WarpX::GetInstance().gett_new(lev);
    const auto& time_of_last_gal_shift = WarpX::GetInstance().time_of_last_gal_shift;
    amrex::Real time_shift = (cur_time - time_of_last_gal_shift);
    amrex::Array<amrex::Real,3> galilean_shift = { v_galilean[0]*time_shift, v_galilean[1]*time_shift, v_galilean[2]*time_shift };
    const std::array<amrex::Real, 3>& xyzmin = WarpX::LowerCorner(box, galilean_shift, lev);
    m_xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};

    m_galerkin_interpolation = WarpX::galerkin_interpolation;
    m_nox = WarpX::nox;
    m_n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;

    m_lo = amrex::lbound(box);
}
