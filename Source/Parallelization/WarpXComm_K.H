/* Copyright 2019 Axel Huebl, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMM_K_H_
#define WARPX_COMM_K_H_

#include <AMReX.H>
#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp (int j, int k, int l,
                   amrex::Array4<amrex::Real      > const& arr_aux,
                   amrex::Array4<amrex::Real const> const& arr_fine,
                   amrex::Array4<amrex::Real const> const& arr_coarse,
                   const amrex::IntVect& arr_stag,
                   const amrex::IntVect& rr)
{
    using namespace amrex;

    // NOTE Indices (j,k,l) in the following refer to (x,z,-) in 2D and (x,y,z) in 3D

    // Refinement ratio
    const int rj = rr[0];
    const int rk = rr[1];
    const int rl = (AMREX_SPACEDIM == 2) ? 1 : rr[2];

    // Staggering (0: cell-centered; 1: nodal)
    const int sj = arr_stag[0];
    const int sk = arr_stag[1];
    const int sl = (AMREX_SPACEDIM == 2) ? 0 : arr_stag[2];

    // Number of points used for interpolation from coarse grid to fine grid
    const int nj = (sj == 0) ? 1 : 2;
    const int nk = (sk == 0) ? 1 : 2;
    const int nl = (sl == 0) ? 1 : 2;

    const int jc = amrex::coarsen(j, rj);
    const int kc = amrex::coarsen(k, rk);
    const int lc = amrex::coarsen(l, rl);

    amrex::Real wj;
    amrex::Real wk;
    amrex::Real wl;

    // Interpolate from coarse grid to fine grid using either 1 point with weight 1, if both grids
    // are cell-centered, or 2 points with weights depending on the distance, if both grids are nodal
    amrex::Real res = 0.0_rt;
    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                wj = (sj == 0) ? 1.0_rt : (rj - amrex::Math::abs(j - (jc + jj) * rj))
                                          / static_cast<amrex::Real>(rj);
                wk = (sk == 0) ? 1.0_rt : (rk - amrex::Math::abs(k - (kc + kk) * rk))
                                          / static_cast<amrex::Real>(rk);
                wl = (sl == 0) ? 1.0_rt : (rl - amrex::Math::abs(l - (lc + ll) * rl))
                                          / static_cast<amrex::Real>(rl);
                res += wj * wk * wl * arr_coarse(jc+jj,kc+kk,lc+ll);
            }
        }
    }
    arr_aux(j,k,l) = arr_fine(j,k,l) + res;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp (int j, int k, int l,
                   amrex::Array4<amrex::Real      > const& arr_aux,
                   amrex::Array4<amrex::Real const> const& arr_fine,
                   amrex::Array4<amrex::Real const> const& arr_coarse,
                   amrex::Array4<amrex::Real const> const& arr_tmp,
                   const amrex::IntVect& arr_fine_stag,
                   const amrex::IntVect& arr_coarse_stag,
                   const amrex::IntVect& rr)
{
    using namespace amrex;

    // NOTE Indices (j,k,l) in the following refer to (x,z,-) in 2D and (x,y,z) in 3D

    // Refinement ratio
    const int rj = rr[0];
    const int rk = rr[1];
    const int rl = (AMREX_SPACEDIM == 2) ? 1 : rr[2];

    // Staggering of fine array (0: cell-centered; 1: nodal)
    const int sj_fp = arr_fine_stag[0];
    const int sk_fp = arr_fine_stag[1];
    const int sl_fp = (AMREX_SPACEDIM == 2) ? 0 : arr_fine_stag[2];

    // Staggering of coarse array (0: cell-centered; 1: nodal)
    const int sj_cp = arr_coarse_stag[0];
    const int sk_cp = arr_coarse_stag[1];
    const int sl_cp = (AMREX_SPACEDIM == 2) ? 0 : arr_coarse_stag[2];

    // Number of points used for interpolation from coarse grid to fine grid
    int nj;
    int nk;
    int nl;

    int jc = amrex::coarsen(j, rj);
    int kc = amrex::coarsen(k, rk);
    int lc = amrex::coarsen(l, rl);

    amrex::Real tmp = 0.0_rt;
    amrex::Real fine = 0.0_rt;
    amrex::Real coarse = 0.0_rt;

    amrex::Real wj;
    amrex::Real wk;
    amrex::Real wl;

    // 1) Interpolation from coarse nodal to fine nodal

    nj = 2;
    nk = 2;
    nl = (AMREX_SPACEDIM == 2) ? 1 : 2;

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                wj = (rj - amrex::Math::abs(j - (jc + jj) * rj)) / static_cast<amrex::Real>(rj);
                wk = (rk - amrex::Math::abs(k - (kc + kk) * rk)) / static_cast<amrex::Real>(rk);
                wl = (AMREX_SPACEDIM == 2) ? 1.0_rt : (rl - amrex::Math::abs(l - (lc + ll) * rl))
                                                      / static_cast<amrex::Real>(rl);
                tmp += wj * wk * wl * arr_tmp(jc+jj,kc+kk,lc+ll);
            }
        }
    }

    // 2) Interpolation from coarse staggered to fine nodal

    nj = 2;
    nk = 2;
    nl = (AMREX_SPACEDIM == 2) ? 1 : 2;

    const int jn = (sj_cp == 1) ? j : j - rj / 2;
    const int kn = (sk_cp == 1) ? k : k - rk / 2;
    const int ln = (sl_cp == 1) ? l : l - rl / 2;

    jc = amrex::coarsen(jn, rj);
    kc = amrex::coarsen(kn, rk);
    lc = amrex::coarsen(ln, rl);

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                wj = (rj - amrex::Math::abs(jn - (jc + jj) * rj)) / static_cast<amrex::Real>(rj);
                wk = (rk - amrex::Math::abs(kn - (kc + kk) * rk)) / static_cast<amrex::Real>(rk);
                wl = (AMREX_SPACEDIM == 2) ? 1.0_rt : (rl - amrex::Math::abs(ln - (lc + ll) * rl))
                                                      / static_cast<amrex::Real>(rl);
                coarse += wj * wk * wl * arr_coarse(jc+jj,kc+kk,lc+ll);
            }
        }
    }

    // 3) Interpolation from fine staggered to fine nodal

    nj = (sj_fp == 0) ? 2 : 1;
    nk = (sk_fp == 0) ? 2 : 1;
    nl = (sl_fp == 0 && AMREX_SPACEDIM == 3) ? 2 : 1;

    const int jm = (sj_fp == 0) ? j-1 : j;
    const int km = (sk_fp == 0) ? k-1 : k;
    const int lm = (sl_fp == 0 && AMREX_SPACEDIM == 3) ? l-1 : l;

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                wj = 1.0_rt / static_cast<amrex::Real>(nj);
                wk = 1.0_rt / static_cast<amrex::Real>(nk);
                wl = 1.0_rt / static_cast<amrex::Real>(nl);
                fine += wj * wk * wl * arr_fine(jm+jj,km+kk,lm+ll);
            }
        }
    }

    // Final result
    arr_aux(j,k,l) = tmp + (fine - coarse);
}

/**
 * \brief Arbitrary-order interpolation function used to center a given MultiFab between two grids
 * with different staggerings. With the FDTD solver, this performs simple linear interpolation.
 * With the PSATD solver, this performs arbitrary-order interpolation based on the Fornberg coefficients.
 * The result is stored in the output array \c dst_arr.
 *
 * \param[in] j index along x of the output array
 * \param[in] k index along y (in 3D) or z (in 2D) of the output array
 * \param[in] l index along z (in 3D, \c l = 0 in 2D) of the output array
 * \param[in,out] dst_arr output array where interpolated values are stored
 * \param[in] src_arr input array storing the values used for interpolation
 * \param[in] dst_stag \c IndexType of the output array
 * \param[in] src_stag \c IndexType of the input array
 * \param[in] nox order of finite-order centering along x
 * \param[in] noy order of finite-order centering along y
 * \param[in] noz order of finite-order centering along z
 * \param[in] stencil_coeffs_x array of ordered Fornberg coefficients for finite-order centering stencil along x
 * \param[in] stencil_coeffs_y array of ordered Fornberg coefficients for finite-order centering stencil along y
 * \param[in] stencil_coeffs_z array of ordered Fornberg coefficients for finite-order centering stencil along z
 */
template< bool IS_PSATD = false >
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp (const int j,
                   const int k,
                   const int l,
                   amrex::Array4<amrex::Real      > const& dst_arr,
                   amrex::Array4<amrex::Real const> const& src_arr,
                   const amrex::IntVect& dst_stag,
                   const amrex::IntVect& src_stag,
                   const int nox = 2,
                   const int noy = 2,
                   const int noz = 2,
                   amrex::Real const* stencil_coeffs_x = nullptr,
                   amrex::Real const* stencil_coeffs_y = nullptr,
                   amrex::Real const* stencil_coeffs_z = nullptr)
{
    using namespace amrex;

    // Pad input array with zeros beyond ghost cells
    // for out-of-bound accesses due to large-stencil operations
    const auto src_arr_zeropad = [src_arr] (const int jj, const int kk, const int ll) noexcept
    {
        return src_arr.contains(jj,kk,ll) ? src_arr(jj,kk,ll) : 0.0_rt;
    };

    // Avoid compiler warnings
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(noy, stencil_coeffs_y);
#endif

    // Avoid compiler warnings
#ifndef WARPX_USE_PSATD
    amrex::ignore_unused(stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
#endif

    // If dst_nodal = true , we are centering from a staggered grid to a nodal grid
    // If dst_nodal = false, we are centering from a nodal grid to a staggered grid
    const bool dst_nodal = (dst_stag == amrex::IntVect::TheNodeVector());

    // See 1D examples below to understand the meaning of this integer shift
    const int shift = (dst_nodal) ? 0 : 1;

    // Staggering (s = 0 if cell-centered, s = 1 if nodal)
    const int sj = (dst_nodal) ? src_stag[0] : dst_stag[0];
    const int sk = (dst_nodal) ? src_stag[1] : dst_stag[1];
#if (AMREX_SPACEDIM == 3)
    const int sl = (dst_nodal) ? src_stag[2] : dst_stag[2];
#endif

    // Interpolate along j,k,l only if source MultiFab is staggered along j,k,l
    const bool interp_j = (sj == 0);
    const bool interp_k = (sk == 0);
#if (AMREX_SPACEDIM == 3)
    const bool interp_l = (sl == 0);
#endif

    const int noj = nox;
#if   (AMREX_SPACEDIM == 2)
    const int nok = noz;
#elif (AMREX_SPACEDIM == 3)
    const int nok = noy;
    const int nol = noz;
#endif

    // Additional normalization factor
    const amrex::Real wj = (interp_j) ? 0.5_rt : 1.0_rt;
    const amrex::Real wk = (interp_k) ? 0.5_rt : 1.0_rt;
#if   (AMREX_SPACEDIM == 2)
    constexpr amrex::Real wl = 1.0_rt;
#elif (AMREX_SPACEDIM == 3)
    const amrex::Real wl = (interp_l) ? 0.5_rt : 1.0_rt;
#endif

    // Min and max for interpolation loop along j
    const int jmin = (interp_j) ? j - noj/2 + shift     : j;
    const int jmax = (interp_j) ? j + noj/2 + shift - 1 : j;

    // Min and max for interpolation loop along k
    const int kmin = (interp_k) ? k - nok/2 + shift     : k;
    const int kmax = (interp_k) ? k + nok/2 + shift - 1 : k;

    // Min and max for interpolation loop along l
#if   (AMREX_SPACEDIM == 2)
    // l = 0 always
    const int lmin = l;
    const int lmax = l;
#elif (AMREX_SPACEDIM == 3)
    const int lmin = (interp_l) ? l - nol/2 + shift     : l;
    const int lmax = (interp_l) ? l + nol/2 + shift - 1 : l;
#endif

    // Number of interpolation points
    const int nj = jmax - jmin;
    const int nk = kmax - kmin;
    const int nl = lmax - lmin;

    // Example of 1D centering from nodal grid to nodal grid (simple copy):
    //
    //         j
    // --o-----o-----o--  result(j) = f(j)
    // --o-----o-----o--
    //  j-1    j    j+1
    //
    // Example of 1D linear centering from staggered grid to nodal grid:
    //
    //         j
    // --o-----o-----o--  result(j) = (f(j-1) + f(j)) / 2
    // -----x-----x-----
    //     j-1    j
    //
    // Example of 1D linear centering from nodal grid to staggered grid:
    // (note the shift of +1 in the indices with respect to the case above, see variable "shift")
    //
    //         j
    // --x-----x-----x--  result(j) = (f(j) + f(j+1)) / 2
    // -----o-----o-----
    //      j    j+1
    //
    // Example of 1D finite-order centering from staggered grid to nodal grid:
    //
    //                     j
    // --o-----o-----o-----o-----o-----o-----o--  result(j) = c_0 * (f(j-1) + f(j)  ) / 2
    // -----x-----x-----x-----x-----x-----x-----            + c_1 * (f(j-2) + f(j+1)) / 2
    //     j-3   j-2   j-1    j    j+1   j+2                + c_2 * (f(j-3) + f(j+2)) / 2
    //     c_2   c_1   c_0   c_0   c_1   c_2                + ...
    //
    // Example of 1D finite-order centering from nodal grid to staggered grid:
    // (note the shift of +1 in the indices with respect to the case above, see variable "shift")
    //
    //                     j
    // --x-----x-----x-----x-----x-----x-----x--  result(j) = c_0 * (f(j)   + f(j+1)) / 2
    // -----o-----o-----o-----o-----o-----o-----            + c_1 * (f(j-1) + f(j+2)) / 2
    //     j-2   j-1    j    j+1   j+2   j+3                + c_2 * (f(j-2) + f(j+3)) / 2
    //     c_2   c_1   c_0   c_0   c_1   c_2                + ...

    amrex::Real res = 0.0_rt;

    if( !IS_PSATD ) // FDTD (linear interpolation)
    {
        for (int ll = 0; ll <= nl; ll++)
        {
            for (int kk = 0; kk <= nk; kk++)
            {
                for (int jj = 0; jj <= nj; jj++)
                {
                    res += src_arr_zeropad(jmin+jj,kmin+kk,lmin+ll);
                }
            }
        }

    }

    else // PSATD (finite-order interpolation)
    {
        amrex::Real cj = 1.0_rt;
        amrex::Real ck = 1.0_rt;
        amrex::Real cl = 1.0_rt;

        amrex::Real const* scj = stencil_coeffs_x;
#if (AMREX_SPACEDIM == 2)
        amrex::Real const* sck = stencil_coeffs_z;
#elif (AMREX_SPACEDIM == 3)
        amrex::Real const* sck = stencil_coeffs_y;
        amrex::Real const* scl = stencil_coeffs_z;
#endif

        for (int ll = 0; ll <= nl; ll++)
        {
#if (AMREX_SPACEDIM == 3)
            if (interp_l) cl = scl[ll];
#endif
            for (int kk = 0; kk <= nk; kk++)
            {
                if (interp_k) ck = sck[kk];

                for (int jj = 0; jj <= nj; jj++)
                {
                    if (interp_j) cj = scj[jj];

                    res += cj * ck * cl * src_arr_zeropad(jmin+jj,kmin+kk,lmin+ll);
                }
            }
        }
    }

    dst_arr(j,k,l) = wj * wk * wl * res;
}

#endif
