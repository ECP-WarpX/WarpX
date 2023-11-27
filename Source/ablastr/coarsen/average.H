/* Copyright 2022 Edoardo Zoni, Remi Lehe, Prabhat Kumar, Axel Huebl
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_COARSEN_AVERAGE_H_
#define ABLASTR_COARSEN_AVERAGE_H_

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_BaseFwd.H>

#include <cstdlib>


/** Mesh Coarsening by Averaging
 *
 * These methods are mostly used for mesh-refinement.
 */
namespace ablastr::coarsen::average
{
    /**
     * \brief Interpolates the floating point data contained in the source Array4
     *        \c arr_src, extracted from a fine MultiFab, with weights defined in
     *        such a way that the total charge is preserved.
     *
     * Works for any staggering of the fine source and coarse d
     *
     * \param[in] arr_src       source fine floating point data to be interpolated
     * \param[in] stag_src_fine staggering of the source fine MultiFab
     * \param[in] stag_des_crse staggering of the destination coarsened MultiFab
     * \param[in] crse_ratio    coarsening ratio along each spatial direction
     * \param[in] i             index along x of the coarsened Array4 to be filled
     * \param[in] j             index along y of the coarsened Array4 to be filled
     * \param[in] k             index along z of the coarsened Array4 to be filled
     * \param[in] comp          index along the fourth component of the Array4 \c arr_src
     *                          containing the data to be interpolated
     *
     * \return Interpolated field at cell (i,j,k) of a coarsened Array4
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    Interp (
        amrex::Array4<amrex::Real const> const &arr_src,
        amrex::GpuArray<int, 3> const &stag_src_fine,
        amrex::GpuArray<int, 3> const &stag_des_crse,
        amrex::GpuArray<int, 3> const &crse_ratio,
        int const i, int const j, int const k, int const comp)
    {
        // Number of points and starting indices of source array (fine)
        amrex::GpuArray<int, 3> np, idx_min, use_half;
        amrex::Real crx_cry_crz_inv;

        CalculateCoarseningData(idx_min, np, use_half, crx_cry_crz_inv, stag_des_crse, stag_src_fine, crse_ratio);

        return InterpWithCoarseningData(arr_src, idx_min, np, use_half, crx_cry_crz_inv, crse_ratio, i, j, k, comp);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void
    CalculateCoarseningData (
        amrex::GpuArray<int, 3> &src_fine_index_min,
        amrex::GpuArray<int, 3> &num_src_points,
        amrex::GpuArray<int, 3> &use_half_weight, // Use int instead bool?
        amrex::Real &crx_cry_crz_inv,
        amrex::GpuArray<int, 3> const &stag_src_fine,
        amrex::GpuArray<int, 3> const &stag_des_crse,
        amrex::GpuArray<int, 3> const &crse_ratio)
    {
        for (int l = 0; l < 3; ++l) {
            int two_times_index_min = -crse_ratio[l]*stag_des_crse[l] + stag_src_fine[l] - 1;
            if (two_times_index_min % 2 == 0) {
                src_fine_index_min[l] = two_times_index_min/2;
                num_src_points[l] = crse_ratio[l]+1;
                use_half_weight[l] = 1; // True
            } else {
                src_fine_index_min[l] = (two_times_index_min+1)/2;
                num_src_points[l] = crse_ratio[l];
                use_half_weight[l] = 0; // False
            }
        }
        crx_cry_crz_inv = 1.0_rt / static_cast<amrex::Real>(crse_ratio[0]*crse_ratio[1]*crse_ratio[2]);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    InterpWithCoarseningData (
        amrex::Array4<amrex::Real const> const &arr_src,
        amrex::GpuArray<int, 3> const &src_fine_index_min,
        amrex::GpuArray<int, 3> const &num_src_points,
        amrex::GpuArray<int, 3> const &use_half,
        amrex::Real crx_cry_crz_inv,
        amrex::GpuArray<int, 3> const &coarsen_ratio,
        int i, int j, int k, int comp)
    {
        using namespace amrex::literals;

        // Auxiliary integer variables
        int const numx = num_src_points[0];
        int const numy = num_src_points[1];
        int const numz = num_src_points[2];
        int const iimin = i * coarsen_ratio[0] + src_fine_index_min[0];
        int const jjmin = j * coarsen_ratio[1] + src_fine_index_min[1];
        int const kkmin = k * coarsen_ratio[2] + src_fine_index_min[2];

        // Add neutral elements (=0) beyond guard cells in source array (fine)
        auto const arr_src_safe = [arr_src]
                AMREX_GPU_DEVICE(int const ix, int const iy, int const iz, int const n) noexcept {
            return arr_src.contains(ix, iy, iz) ? arr_src(ix, iy, iz, n) : 0.0_rt;
        };

        // Interpolate over points computed above. Weights are computed in order
        // to guarantee total charge conservation for any staggering.
        // Python script Source/Utils/check_interp_points_and_weights.py can be
        // used to check interpolation points and weights in 1D.
        amrex::Real c = 0.0_rt;
        for (int kref = 0; kref < numz; ++kref) {
            int kk = kkmin + kref;
            // kk = amrex::Clamp(kk, arr_src.begin.z, arr_src.end.z - 1);
            amrex::Real kfactor = (use_half[2] && (kref == 0 || kref == numz-1)) ? 0.5_rt : 1.0_rt;

            for (int jref = 0; jref < numy; ++jref) {
                int jj = jjmin + jref;
                // jj = amrex::Clamp(jj, arr_src.begin.y, arr_src.end.y - 1);
                amrex::Real jfactor = (use_half[1] && (jref == 0 || jref == numy-1)) ? 0.5_rt : 1.0_rt;

                for (int iref = 0; iref < numx; ++iref) {
                    int ii = iimin + iref;
                    // ii = amrex::Clamp(ii, arr_src.begin.x, arr_src.end.x - 1);
                    amrex::Real ifactor = (use_half[0] && (iref == 0 || iref == numx-1)) ? 0.5_rt : 1.0_rt;

                    c += ifactor * jfactor * kfactor * arr_src_safe(ii, jj, kk, comp);
                }
            }
        }

        c *= crx_cry_crz_inv;

        return c;
    }

    /**
     * \brief Loops over the boxes of the coarsened MultiFab \c mf_dst and fills
     *        them by interpolating the data contained in the fine MultiFab \c mf_src.
     *
     * \param[in,out] mf_dst     coarsened MultiFab containing the floating point data
     *                           to be filled by interpolating the source fine MultiFab
     * \param[in]     mf_src     fine MultiFab containing the floating point data to be interpolated
     * \param[in]     ncomp      number of components to loop over for the coarsened
     *                           Array4 extracted from the coarsened MultiFab \c mf_dst
     * \param[in]     ngrow      number of guard cells to fill along each spatial direction
     * \param[in]     crse_ratio coarsening ratio between the fine MultiFab \c mf_src
     *                           and the coarsened MultiFab \c mf_dst along each spatial direction
     */
    void
    Loop (
        amrex::MultiFab & mf_dst,
        amrex::MultiFab const & mf_src,
        int ncomp,
        amrex::IntVect ngrow,
        amrex::IntVect crse_ratio
    );

    /**
     * \brief Stores in the coarsened MultiFab \c mf_dst the values obtained by
     *        interpolating the data contained in the fine MultiFab \c mf_src.
     *
     * \param[in,out] mf_dst     coarsened MultiFab containing the floating point data
     *                           to be filled by interpolating the fine MultiFab \c mf_src
     * \param[in]     mf_src     fine MultiFab containing the floating point data to be interpolated
     * \param[in]     crse_ratio coarsening ratio between the fine MultiFab \c mf_src
     *                           and the coarsened MultiFab \c mf_dst along each spatial direction
     */
    void
    Coarsen (
        amrex::MultiFab & mf_dst,
        amrex::MultiFab const & mf_src,
        amrex::IntVect crse_ratio
    );

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    InterpWithWeights (
        amrex::Array4<amrex::Real const> const &arr_src,
        amrex::Vector<amrex::Real> const &weights_x,
        amrex::Vector<amrex::Real> const &weights_y,
        amrex::Vector<amrex::Real> const &weights_z,
        amrex::GpuArray<int, 3> const &src_fine_ind_min,
        amrex::GpuArray<int, 3> const &coarsen_ratio,
        int i, int j, int k, int comp)
    {
        const int iimin = i * coarsen_ratio[0] + src_fine_ind_min[0];
        const int jjmin = j * coarsen_ratio[1] + src_fine_ind_min[1];
        const int kkmin = k * coarsen_ratio[2] + src_fine_ind_min[2];

        auto const arr_src_safe = [arr_src]
                AMREX_GPU_DEVICE(int const ix, int const iy, int const iz, int const n) noexcept {
            return arr_src.contains(ix, iy, iz) ? arr_src(ix, iy, iz, n) : 0.0_rt;
        };

        amrex::Real c = 0.0_rt;
        for (int kref = 0; k < weights_z.size(); ++kref)
        {
            int kk = kkmin + kref;
            // kk = amrex::Clamp(kk, arr_src.begin.z, arr_src.end.z - 1);
            amrex::Real wz = weights_z[kref];

            for (int jref = 0; j < weights_y.size(); ++jref)
            {
                int jj = jjmin + jref;
                // jj = amrex::Clamp(jj, arr_src.begin.y, arr_src.end.y - 1);
                amrex::Real wy = weights_y[jref];

                for (int iref = 0; i < weights_x.size(); ++iref)
                {
                    int ii = iimin + iref;
                    // ii = amrex::Clamp(ii, arr_src.begin.x, arr_src.end.x - 1);
                    amrex::Real wx = weights_x[iref];

                    c += wx*wy*wz*arr_src_safe(ii, jj, kk, comp);
                }
            }
        }
        return c;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void
    CoarseningPointsAndWeights (
        amrex::GpuArray<int, 3> &src_fine_index_min,
        amrex::Vector<amrex::Real> &weights_x,
        amrex::Vector<amrex::Real> &weights_y,
        amrex::Vector<amrex::Real> &weights_z,
        amrex::GpuArray<int, 3> const &coarsen_ratio,
        amrex::GpuArray<int, 3> const &stag_src_fine,
        amrex::GpuArray<int, 3> const &stag_des_crse)
    {
        using namespace amrex::literals;

        int num_points[3];
        bool useHalf[3];

        for (int l = 0; l < 3; ++l) {
            int const twoImin = -coarsen_ratio[l]*stag_des_crse[l] + stag_src_fine[l] - 1;
            if (twoImin % 2 == 0) {
                src_fine_index_min[l] = twoImin/2;
                num_points[l] = coarsen_ratio[l]+1;
                useHalf[l] = true;
            } else {
                src_fine_index_min[l] = (twoImin+1)/2;
                num_points[l] = coarsen_ratio[l];
                useHalf[l] = false;
            }
        }

        weights_x.resize(num_points[0], 1.0_rt);
        weights_y.resize(num_points[1], 1.0_rt);
        weights_z.resize(num_points[2], 1.0_rt);

        if (useHalf[0]) {
            weights_x[0] *= 0.5_rt;
            weights_x[num_points[0]-1] *= 0.5_rt;
        }
        if (useHalf[1]) {
            weights_y[0] *= 0.5_rt;
            weights_y[num_points[1]-1] *= 0.5_rt;
        }
        if (useHalf[2]) {
            weights_z[0] *= 0.5_rt;
            weights_z[num_points[2]-1] *= 0.5_rt;
        }

    }

    template <class U, typename std::enable_if<!std::is_void<U>::value,int>::type = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U& Array4ClampedEval (amrex::Array4<U> const& arr,
                          int i, int j, int k, int n)
    {
        int iClamped = amrex::Clamp(i, arr.begin.x, arr.end.x - 1);
        int jClamped = amrex::Clamp(j, arr.begin.y, arr.end.y - 1);
        int kClamped = amrex::Clamp(k, arr.begin.x, arr.end.z - 1);
        return arr(iClamped,jClamped,kClamped,n);
    }

} // namespace ablastr::coarsen::average

#endif // ABLASTR_COARSEN_AVERAGE_H_
