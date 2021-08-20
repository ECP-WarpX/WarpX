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

/*
 * \brief Set a tilebox so that it only covers the lower/upper half of the guard cells.
 *
 * \param[in,out] tb tilebox to be modified
 * \param[in] dir direction where the tilebox smallEnd/bigEnd is modified
 * \param[in] n_domain number of cells in the whole simulation domain
 * \param[in] tb_smallEnd smallEnd of the tilebox
 * \param[in] tb_bigEnd bigEnd of the tilebox
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void constrain_tilebox_to_guards(
    amrex::Box& tb,
    const int dir,
    const int n_domain,
    const int tb_smallEnd,
    const int tb_bigEnd)
{
    using namespace amrex;

    const int n_tile = tb_bigEnd;

    if (tb_smallEnd < 0)
    {
        const int n_guard = -tb_smallEnd;
        tb.setBig(dir, 0 - n_guard/2 - 1);
    }
    else if (n_tile > n_domain)
    {
        const int n_guard = n_tile - n_domain;
        tb.setSmall(dir, n_domain + n_guard/2 + 1);
    }
}

/*
 * \brief Damp a given field in the guard cells along a given direction
 *
 * \param[in,out] mf_arr array that contains the field values to be damped
 * \oaram[in] i index along x
 * \oaram[in] j index along y (in 3D) or z (in 2D/RZ)
 * \oaram[in] k index along z (in 3D, \c k = 0 in 2D/RZ)
 * \param[in] icomp index along the fourth component of the array
 * \param]in] dir direction where the field will be damped
 * \param[in] n_domain number of cells in the whole simulation domain
 * \param[in] tb_smallEnd smallEnd of the current tilebox
 * \param[in] tb_bigEnd bigEnd of the current tilebox
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void damp_field_in_guards(
    amrex::Array4<amrex::Real> const& mf_arr,
    const int i,
    const int j,
    const int k,
    const int icomp,
    const int dir,
    const int n_domain,
    const int tb_smallEnd,
    const int tb_bigEnd)
{
    using namespace amrex;

    // dir = 0: idx = i (x)
    // dir = 1: idx = j (y in 3D, z in 2D/RZ)
    // dir = 2: idx = k (z in 3D)
    const int idx = ((dir == 0) ? i : ((dir == 1) ? j : k));

    const int n_tile = tb_bigEnd;

    if (tb_smallEnd < 0)
    {
        // Apply damping factor in guards cells below the lower end of the domain
        const int n_guard = -tb_smallEnd;

        const amrex::Real cell = static_cast<amrex::Real>(idx + n_guard);

        const amrex::Real phase = MathConst::pi * cell / n_guard;
        const amrex::Real sin_phase = std::sin(phase);
        const amrex::Real damp_factor = sin_phase * sin_phase;

        mf_arr(i,j,k,icomp) *= damp_factor;
    }
    else if (n_tile > n_domain)
    {
        // Apply damping factor in guards cells above the upper end of the domain
        const int n_guard = n_tile - n_domain;

        const amrex::Real cell = static_cast<amrex::Real>(n_tile - idx);

        const amrex::Real phase = MathConst::pi * cell / n_guard;
        const amrex::Real sin_phase = std::sin(phase);
        const amrex::Real damp_factor = sin_phase * sin_phase;

        mf_arr(i,j,k,icomp) *= damp_factor;
    }
}

#endif
