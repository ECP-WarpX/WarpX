/* Copyright 2019-2020
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WarpXPushFieldsEM_K_h
#define WarpXPushFieldsEM_K_h

#include "Utils/WarpXConst.H"
#include <AMReX.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_tilebox(
    amrex::Box& tb,
    const int zdir,
    const int nz_domain,
    const int tb_smallEnd_z,
    const int tb_bigEnd_z)
{
    using namespace amrex;

    const int nz_tile = tb_bigEnd_z;

    // Set the tileboxes so that they only cover the lower/upper half of the guard cells
    if (tb_smallEnd_z < 0)
    {
        const int nz_guard = -tb_smallEnd_z;
        tb.setBig(zdir, 0 - nz_guard/2 - 1);
    }
    else if (nz_tile > nz_domain)
    {
        const int nz_guard = nz_tile - nz_domain;
        tb.setSmall(zdir, nz_domain + nz_guard/2 + 1);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void damp_field(
    amrex::Array4<amrex::Real> const& mf_arr,
    const int i,
    const int j,
    const int k,
    const int icomp,
    const int nz_domain,
    const int tb_smallEnd_z,
    const int tb_bigEnd_z)
{
    using namespace amrex;

#if   (AMREX_SPACEDIM == 2)
    const int zidx = j;
#elif (AMREX_SPACEDIM == 3)
    const int zidx = k;
#endif

    const int nz_tile = tb_bigEnd_z;

    if (tb_smallEnd_z < 0)
    {
        // Apply damping factor in guards cells below the lower end of the domain
        const int nz_guard = -tb_smallEnd_z;

        const amrex::Real zcell = static_cast<amrex::Real>(zidx + nz_guard);

        const amrex::Real phase = MathConst::pi * zcell / nz_guard;
        const amrex::Real sin_phase = std::sin(phase);
        const amrex::Real damp_factor = sin_phase * sin_phase;

        mf_arr(i,j,k,icomp) *= damp_factor;
    }
    else if (nz_tile > nz_domain)
    {
        // Apply damping factor in guards cells above the upper end of the domain
        const int nz_guard = nz_tile - nz_domain;

        const amrex::Real zcell = static_cast<amrex::Real>(nz_tile - zidx);

        const amrex::Real phase = MathConst::pi * zcell / nz_guard;
        const amrex::Real sin_phase = std::sin(phase);
        const amrex::Real damp_factor = sin_phase * sin_phase;

        mf_arr(i,j,k,icomp) *= damp_factor;
    }
}

#endif
