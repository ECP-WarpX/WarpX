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

/**
 * \brief Interpolation function called within WarpX::UpdateAuxilaryDataSameType
 * with electromagnetic solver to interpolate data from the coarse and fine grids to the fine aux grid,
 * assuming that all grids have the same staggering (either collocated or staggered).
 *
 * \param[in] j index along x of the output array
 * \param[in] k index along y (in 3D) or z (in 2D) of the output array
 * \param[in] l index along z (in 3D, l=0 in 2D) of the output array
 * \param[in,out] arr_aux output array where interpolated values are stored
 * \param[in] arr_fine input fine-patch array storing the values to interpolate
 * \param[in] arr_coarse input coarse-patch array storing the values to interpolate
 * \param[in] arr_stag IndexType of the arrays
 * \param[in] rr mesh refinement ratios along each direction
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp (int j, int k, int l,
                   amrex::Array4<amrex::Real      > const& arr_aux,
                   amrex::Array4<amrex::Real const> const& arr_fine,
                   amrex::Array4<amrex::Real const> const& arr_coarse,
                   const amrex::IntVect& arr_stag,
                   const amrex::IntVect& rr)
{
    using namespace amrex;

    // Pad arr_coarse with zeros beyond ghost cells for out-of-bound accesses
    const auto arr_coarse_zeropad = [arr_coarse] (const int jj, const int kk, const int ll) noexcept
    {
        return arr_coarse.contains(jj,kk,ll) ? arr_coarse(jj,kk,ll) : 0.0_rt;
    };

    // NOTE Indices (j,k,l) in the following refer to (z,-,-) in 1D, (x,z,-) in 2D, and (x,y,z) in 3D

    // Refinement ratio
    const int rj = rr[0];
    const int rk = (AMREX_SPACEDIM == 1) ? 1 : rr[1];
    const int rl = (AMREX_SPACEDIM <= 2) ? 1 : rr[2];

    // Staggering (0: cell-centered; 1: nodal)
    const int sj = arr_stag[0];
    const int sk = (AMREX_SPACEDIM == 1) ? 0 : arr_stag[1];
    const int sl = (AMREX_SPACEDIM <= 2) ? 0 : arr_stag[2];

    // Number of points used for interpolation from coarse grid to fine grid
    const int nj = 2;
    const int nk = 2;
    const int nl = 2;

    const int jc = (sj == 0) ? amrex::coarsen(j - rj/2, rj) : amrex::coarsen(j, rj);
    const int kc = (sk == 0) ? amrex::coarsen(k - rk/2, rk) : amrex::coarsen(k, rk);
    const int lc = (sl == 0) ? amrex::coarsen(l - rl/2, rl) : amrex::coarsen(l, rl);

    // Interpolate from coarse grid to fine grid using 2 points
    // with weights depending on the distance, for both nodal and cell-centered grids
    const amrex::Real hj = (sj == 0) ? 0.5_rt : 0._rt;
    const amrex::Real hk = (sk == 0) ? 0.5_rt : 0._rt;
    const amrex::Real hl = (sl == 0) ? 0.5_rt : 0._rt;

    amrex::Real res = 0.0_rt;

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                const amrex::Real wj = (rj - amrex::Math::abs(j + hj - (jc + jj + hj) * rj)) / static_cast<amrex::Real>(rj);
                const amrex::Real wk = (rk - amrex::Math::abs(k + hk - (kc + kk + hk) * rk)) / static_cast<amrex::Real>(rk);
                const amrex::Real wl = (rl - amrex::Math::abs(l + hl - (lc + ll + hl) * rl)) / static_cast<amrex::Real>(rl);
                res += wj * wk * wl * arr_coarse_zeropad(jc+jj,kc+kk,lc+ll);
            }
        }
    }
    arr_aux(j,k,l) = arr_fine(j,k,l) + res;
}

/**
 * \brief Interpolation function called within WarpX::UpdateAuxilaryDataStagToNodal
 * to interpolate data from the coarse and fine grids to the fine aux grid,
 * with momentum-conserving field gathering, hence between grids with different staggering,
 * and assuming that the aux grid is collocated.
 *
 * \param[in] j index along x of the output array
 * \param[in] k index along y (in 3D) or z (in 2D) of the output array
 * \param[in] l index along z (in 3D, l=0 in 2D) of the output array
 * \param[in,out] arr_aux output array where interpolated values are stored
 * \param[in] arr_fine input fine-patch array storing the values to interpolate
 * \param[in] arr_coarse input coarse-patch array storing the values to interpolate
 * \param[in] arr_fine_stag IndexType of the fine-patch arrays
 * \param[in] arr_coarse_stag IndexType of the coarse-patch arrays
 * \param[in] rr mesh refinement ratios along each direction
 */
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

    // Pad input arrays with zeros beyond ghost cells
    // for out-of-bound accesses due to large-stencil operations
    const auto arr_fine_zeropad = [arr_fine] (const int jj, const int kk, const int ll) noexcept
    {
        return arr_fine.contains(jj,kk,ll) ? arr_fine(jj,kk,ll) : 0.0_rt;
    };
    const auto arr_coarse_zeropad = [arr_coarse] (const int jj, const int kk, const int ll) noexcept
    {
        return arr_coarse.contains(jj,kk,ll) ? arr_coarse(jj,kk,ll) : 0.0_rt;
    };
    const auto arr_tmp_zeropad = [arr_tmp] (const int jj, const int kk, const int ll) noexcept
    {
        return arr_tmp.contains(jj,kk,ll) ? arr_tmp(jj,kk,ll) : 0.0_rt;
    };

    // NOTE Indices (j,k,l) in the following refer to:
    // - (z,-,-) in 1D
    // - (x,z,-) in 2D
    // - (r,z,-) in RZ
    // - (x,y,z) in 3D

    // Refinement ratio
    const int rj = rr[0];
#if defined(WARPX_DIM_1D_Z)
    constexpr int rk = 1;
    constexpr int rl = 1;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const     int rk = rr[1];
    constexpr int rl = 1;
#else
    const int rk = rr[1];
    const int rl = rr[2];
#endif

    // Staggering of fine array (0: cell-centered; 1: nodal)
    const int sj_fp = arr_fine_stag[0];
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const int sk_fp = arr_fine_stag[1];
#elif defined(WARPX_DIM_3D)
    const int sk_fp = arr_fine_stag[1];
    const int sl_fp = arr_fine_stag[2];
#endif

    // Staggering of coarse array (0: cell-centered; 1: nodal)
    const int sj_cp = arr_coarse_stag[0];
#if defined(WARPX_DIM_1D_Z)
    constexpr int sk_cp = 0;
    constexpr int sl_cp = 0;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const     int sk_cp = arr_coarse_stag[1];
    constexpr int sl_cp = 0;
#else
    const int sk_cp = arr_coarse_stag[1];
    const int sl_cp = arr_coarse_stag[2];
#endif

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

    // 1) Interpolation from coarse nodal to fine nodal

    nj = 2;
#if defined(WARPX_DIM_1D_Z)
    nk = 1;
    nl = 1;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    nk = 2;
    nl = 1;
#else
    nk = 2;
    nl = 2;
#endif

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                auto c =  arr_tmp_zeropad(jc+jj,kc+kk,lc+ll);
                c *= (rj - amrex::Math::abs(j - (jc + jj) * rj)) / static_cast<amrex::Real>(rj);
#if (AMREX_SPACEDIM >= 2)
                c *= (rk - amrex::Math::abs(k - (kc + kk) * rk)) / static_cast<amrex::Real>(rk);
#endif
#if (AMREX_SPACEDIM == 3)
                c *= (rl - amrex::Math::abs(l - (lc + ll) * rl)) / static_cast<amrex::Real>(rl);
#endif
                tmp += c;
            }
        }
    }

    // 2) Interpolation from coarse staggered to fine nodal

    nj = 2;
#if defined(WARPX_DIM_1D_Z)
    nk = 1;
    nl = 1;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    nk = 2;
    nl = 1;
#else
    nk = 2;
    nl = 2;
#endif

    const int jn = (sj_cp == 1) ? j : j - rj / 2;
    const int kn = (sk_cp == 1) ? k : k - rk / 2;
    const int ln = (sl_cp == 1) ? l : l - rl / 2;

    jc = amrex::coarsen(jn, rj);
    kc = amrex::coarsen(kn, rk);
    lc = amrex::coarsen(ln, rl);

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                auto c = arr_coarse_zeropad(jc+jj,kc+kk,lc+ll);
                c *= (rj - amrex::Math::abs(jn - (jc + jj) * rj)) / static_cast<amrex::Real>(rj);
#if (AMREX_SPACEDIM >= 2)
                c *= (rk - amrex::Math::abs(kn - (kc + kk) * rk)) / static_cast<amrex::Real>(rk);
#endif
#if (AMREX_SPACEDIM == 3)
                c *= (rl - amrex::Math::abs(ln - (lc + ll) * rl)) / static_cast<amrex::Real>(rl);
#endif
                coarse += c;
            }
        }
    }

    // 3) Interpolation from fine staggered to fine nodal

    nj = (sj_fp == 0) ? 2 : 1;
#if defined(WARPX_DIM_1D_Z)
    nk = 1;
    nl = 1;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    nk = (sk_fp == 0) ? 2 : 1;
    nl = 1;
#else
    nk = (sk_fp == 0) ? 2 : 1;
    nl = (sl_fp == 0) ? 2 : 1;
#endif

    const int jm = (sj_fp == 0) ? j-1 : j;
#if defined(WARPX_DIM_1D_Z)
    const int km = k;
    const int lm = l;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const int km = (sk_fp == 0) ? k-1 : k;
    const int lm = l;
#else
    const int km = (sk_fp == 0) ? k-1 : k;
    const int lm = (sl_fp == 0) ? l-1 : l;
#endif

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                fine += arr_fine_zeropad(jm+jj,km+kk,lm+ll);
            }
        }
    }
    fine = fine/static_cast<amrex::Real>(nj*nk*nl);

    // Final result
    arr_aux(j,k,l) = tmp + (fine - coarse);
}
/**
 * \brief Interpolation function called within WarpX::UpdateAuxilaryDataStagToNodal
 * to interpolate data from the coarse and fine grids to the fine aux grid,
 * with momentum-conserving field gathering, hence between grids with different staggering,
 * and assuming that the aux grid is collocated.
 *
 * \param[in] j index along x of the output array
 * \param[in] k index along y (in 3D) or z (in 2D) of the output array
 * \param[in] l index along z (in 3D, l=0 in 2D) of the output array
 * \param[in,out] arr_aux output array where interpolated values are stored
 * \param[in] arr_fine input fine-patch array storing the values to interpolate
 * \param[in] arr_fine_stag IndexType of the fine-patch arrays
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp (int j, int k, int l,
                   amrex::Array4<amrex::Real      > const& arr_aux,
                   amrex::Array4<amrex::Real const> const& arr_fine,
                   const amrex::IntVect& arr_fine_stag)
{
    using namespace amrex;

    // Pad input arrays with zeros beyond ghost cells
    // for out-of-bound accesses due to large-stencil operations
    const auto arr_fine_zeropad = [arr_fine] (const int jj, const int kk, const int ll) noexcept
    {
        return arr_fine.contains(jj,kk,ll) ? arr_fine(jj,kk,ll) : 0.0_rt;
    };

    // NOTE Indices (j,k,l) in the following refer to:
    // - (z,-,-) in 1D
    // - (x,z,-) in 2D
    // - (r,z,-) in RZ
    // - (x,y,z) in 3D

    // Staggering of fine array (0: cell-centered; 1: nodal)
    const int sj_fp = arr_fine_stag[0];
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const int sk_fp = arr_fine_stag[1];
#elif defined(WARPX_DIM_3D)
    const int sk_fp = arr_fine_stag[1];
    const int sl_fp = arr_fine_stag[2];
#endif

    // Number of points used for interpolation from coarse grid to fine grid
    int nj;
    int nk;
    int nl;

    amrex::Real fine = 0.0_rt;

    // 3) Interpolation from fine staggered to fine nodal

    nj = (sj_fp == 0) ? 2 : 1;
#if defined(WARPX_DIM_1D_Z)
    nk = 1;
    nl = 1;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    nk = (sk_fp == 0) ? 2 : 1;
    nl = 1;
#else
    nk = (sk_fp == 0) ? 2 : 1;
    nl = (sl_fp == 0) ? 2 : 1;
#endif

    const int jm = (sj_fp == 0) ? j-1 : j;
#if defined(WARPX_DIM_1D_Z)
    const int km = k;
    const int lm = l;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const int km = (sk_fp == 0) ? k-1 : k;
    const int lm = l;
#else
    const int km = (sk_fp == 0) ? k-1 : k;
    const int lm = (sl_fp == 0) ? l-1 : l;
#endif

    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                fine += arr_fine_zeropad(jm+jj,km+kk,lm+ll);
            }
        }
    }
    fine = fine/static_cast<amrex::Real>(nj*nk*nl);

    // Final result
    arr_aux(j,k,l) = fine;
}

/**
 * \brief Arbitrary-order interpolation function used to center a given MultiFab between two grids
 * with different staggerings. The arbitrary-order interpolation is based on the Fornberg coefficients.
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
#if defined(WARPX_DIM_1D_Z)
    amrex::ignore_unused(nox, noy, stencil_coeffs_x, stencil_coeffs_y);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    amrex::ignore_unused(noy, stencil_coeffs_y);
#endif

    // If dst_nodal = true , we are centering from a staggered grid to a nodal grid
    // If dst_nodal = false, we are centering from a nodal grid to a staggered grid
    const bool dst_nodal = (dst_stag == amrex::IntVect::TheNodeVector());

    // See 1D examples below to understand the meaning of this integer shift
    const int shift = (dst_nodal) ? 0 : 1;

    // Staggering (s = 0 if cell-centered, s = 1 if nodal)
    const int sj = (dst_nodal) ? src_stag[0] : dst_stag[0];
#if (AMREX_SPACEDIM >= 2)
    const int sk = (dst_nodal) ? src_stag[1] : dst_stag[1];
#endif
#if defined(WARPX_DIM_3D)
    const int sl = (dst_nodal) ? src_stag[2] : dst_stag[2];
#endif

    // Interpolate along j,k,l only if source MultiFab is staggered along j,k,l
    const bool interp_j = (sj == 0);
#if (AMREX_SPACEDIM >= 2)
    const bool interp_k = (sk == 0);
#endif
#if defined(WARPX_DIM_3D)
    const bool interp_l = (sl == 0);
#endif

#if   defined(WARPX_DIM_1D_Z)
    const int noj = noz;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const int noj = nox;
    const int nok = noz;
#elif defined(WARPX_DIM_3D)
    const int noj = nox;
    const int nok = noy;
    const int nol = noz;
#endif

    // Additional normalization factor
    const amrex::Real wj = (interp_j) ? 0.5_rt : 1.0_rt;
#if   defined(WARPX_DIM_1D_Z)
    constexpr amrex::Real wk = 1.0_rt;
    constexpr amrex::Real wl = 1.0_rt;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const amrex::Real wk = (interp_k) ? 0.5_rt : 1.0_rt;
    constexpr amrex::Real wl = 1.0_rt;
#elif defined(WARPX_DIM_3D)
    const amrex::Real wk = (interp_k) ? 0.5_rt : 1.0_rt;
    const amrex::Real wl = (interp_l) ? 0.5_rt : 1.0_rt;
#endif

    // Min and max for interpolation loop along j
    const int jmin = (interp_j) ? j - noj/2 + shift     : j;
    const int jmax = (interp_j) ? j + noj/2 + shift - 1 : j;

    // Min and max for interpolation loop along k
#if defined(WARPX_DIM_1D_Z)
    // k = 0 always
    const int kmin = k;
    const int kmax = k;
#else
    const int kmin = (interp_k) ? k - nok/2 + shift     : k;
    const int kmax = (interp_k) ? k + nok/2 + shift - 1 : k;
#endif

    // Min and max for interpolation loop along l
#if (AMREX_SPACEDIM <= 2)
    // l = 0 always
    const int lmin = l;
    const int lmax = l;
#elif defined(WARPX_DIM_3D)
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

#if defined(WARPX_DIM_1D_Z)
    amrex::Real const* scj = stencil_coeffs_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    amrex::Real const* scj = stencil_coeffs_x;
    amrex::Real const* sck = stencil_coeffs_z;
#elif defined(WARPX_DIM_3D)
    amrex::Real const* scj = stencil_coeffs_x;
    amrex::Real const* sck = stencil_coeffs_y;
    amrex::Real const* scl = stencil_coeffs_z;
#endif

    for (int ll = 0; ll <= nl; ll++)
    {
#if defined(WARPX_DIM_3D)
    const amrex::Real cl = (interp_l)? scl[ll] : 1.0_rt;
#else
    const amrex::Real cl = 1.0_rt;
#endif
        for (int kk = 0; kk <= nk; kk++)
        {
#if (AMREX_SPACEDIM >= 2)
            const amrex::Real ck = (interp_k)? sck[kk] : 1.0_rt;
#else
            const amrex::Real ck = 1.0_rt;
#endif
            for (int jj = 0; jj <= nj; jj++)
            {
                const amrex::Real cj =  (interp_j)? scj[jj] : 1.0_rt;

                res += cj * ck * cl * src_arr_zeropad(jmin+jj,kmin+kk,lmin+ll);
            }
        }
    }

    dst_arr(j,k,l) = wj * wk * wl * res;
}

#endif
