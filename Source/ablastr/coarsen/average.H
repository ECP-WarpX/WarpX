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
     * The input (sf) and output (sc) staggering need to be the same.
     *
     * \param[in] arr_src floating point data to be interpolated
     * \param[in] sf      staggering of the source fine MultiFab
     * \param[in] sc      staggering of the destination coarsened MultiFab
     * \param[in] cr      coarsening ratio along each spatial direction
     * \param[in] i       index along x of the coarsened Array4 to be filled
     * \param[in] j       index along y of the coarsened Array4 to be filled
     * \param[in] k       index along z of the coarsened Array4 to be filled
     * \param[in] comp    index along the fourth component of the Array4 \c arr_src
     *                    containing the data to be interpolated
     *
     * \return interpolated field at cell (i,j,k) of a coarsened Array4
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    Interp (
        amrex::Array4<amrex::Real const> const &arr_src,
        amrex::GpuArray<int, 3> const &sf,
        amrex::GpuArray<int, 3> const &sc,
        amrex::GpuArray<int, 3> const &cr,
        int const i, int const j, int const k, int const comp)
    {
        using namespace amrex::literals;

        // Number of points and starting indices of source array (fine)
        int np[3], idx_min[3];
        bool useHalf[3];

        for (int l = 0; l < 3; ++l) {
            const int twoImin = -cr[l]*sc[l] + sf[l] - 1;
            if (twoImin % 2 == 0) {
                idxmin[l] = twoImin/2;
                np[l] = cr[l]+1;
                useHalf[l] = true;
            } else {
                idxmin[l] = (twoImin+1)/2;
                np[l] = cr[l];
                useHalf[l] = false;
            }
        }

        // Auxiliary integer variables
        int const iimin = i * cr[0] + idx_min[0];
        int const jjmin = j * cr[1] + idx_min[1];
        int const kkmin = k * cr[2] + idx_min[2];
        // int const iimax = iimin + np[0] - 1;
        // int const jjmax = jjmin + np[1] - 1;
        // int const kkmax = kkmin + np[2] - 1;

        // Add neutral elements (=0) beyond guard cells in source array (fine)
        // auto const arr_src_safe = [arr_src]
        //         AMREX_GPU_DEVICE(int const ix, int const iy, int const iz, int const n) noexcept {
        //     return arr_src.contains(ix, iy, iz) ? arr_src(ix, iy, iz, n) : 0.0_rt;
        // };

        // Interpolate over points computed above. Weights are computed in order
        // to guarantee total charge conservation for both cell-centered data
        // (equal weights) and nodal data (weights depend on distance between
        // points on fine and coarse grids). Terms multiplied by (1-sf)*(1-sc)
        // are ON for cell-centered data and OFF for nodal data, while terms
        // multiplied by sf*sc are ON for nodal data and OFF for cell-centered data.
        // Python script Source/Utils/check_interp_points_and_weights.py can be
        // used to check interpolation points and weights in 1D.
        amrex::Real c = 0.0_rt;
        for (int kr = 0; kr < np[2]; ++kr) {
            int kk = kkmin + kr;
            kk = amrex::Clamp(kk, arr_src.begin.z, arr_src.end.z - 1);
            amrex::Real kfactor = (useHalf[2] && (kr == 0 || kr == np[2] - 1)) ? 0.5_rt : 1.0_rt;

            for (int jr = 0; jr < np[1]; ++jr) {
                int jj = jjmin + jr;
                jj = amrex::Clamp(jj, arr_src.begin.y, arr_src.end.y - 1);
                amrex::Real jfactor = (useHalf[1] && (jr == 0 || jr == np[1] - 1)) ? 0.5_rt : 1.0_rt;

                for (int ir = 0; ir < np[1]; ++ir) {
                    int ii = iimin + ir;
                    ii = amrex::Clamp(ii, arr_src.begin.x, arr_src.end.x - 1);
                    amrex::Real ifactor = (useHalf[1] && (ir == 0 || ir == np[1] - 1)) ? 0.5_rt : 1.0_rt;

                    c += ifactor*jfactor*kfactor*arr_src(ii,jj,kk,comp);

                    // amrex::Real src_iijjkk = arr_src_safe(ii, jj, kk, comp);
                    // amrex::Real src_iijjkk = Array4ClampedEval(arr_src, ii, jj, kk, comp);
                    // c += ifactor*jfactor*kfactor*src_iijjkk;
                }
            }
        }

        c *= 1.0_rt / static_cast<amrex::Real>(cr[0]*cr[1]*cr[2]);

        return c;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    Interp (
        amrex::Array4<amrex::Real const> const &arr_src,
        amrex::GpuArray<int, 3> const &coarsen_ratio,
        amrex::GpuArray<int, 3> const &ind_min,
        amrex::Vector<amrex::Real> const &weights_x,
        amrex::Vector<amrex::Real> const &weights_y,
        amrex::Vector<amrex::Real> const &weights_z,
        int const i, int const j, int const k, int const comp
    )
    {
        amrex::Real c = 0.0_rt;

        const int iimin = i * coarsen_ratio[0] + ind_min[0];
        const int jjmin = j * coarsen_ratio[1] + ind_min[1];
        const int kkmin = k * coarsen_ratio[2] + ind_min[2];

        for (int kr = 0; k < weights_z.size(); ++kr)
        {
            int kk = kkmin + kr;
            kk = amrex::Clamp(kk, arr_src.begin.z, arr_src.end.z - 1);

            for (int jr = 0; j < weights_y.size(); ++jr)
            {
                int jj = jjmin + jr;
                jj = amrex::Clamp(jj, arr_src.begin.y, arr_src.end.y - 1);

                for (int ir = 0; i < weights_x.size(); ++ir)
                {
                    int ii = iimin + ir;
                    ii = amrex::Clamp(ii, arr_src.begin.x, arr_src.end.x - 1);

                    amrex::Real src_iijjkk = arr_src(ii, jj, kk, comp);
                    c += weights_x[ir]*weights_y[ir]*weights_z[ir]*src_iijjkk;
                }
            }
        }
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

    void
    CoarseningPointsAndWeights (
        amrex::Vector<amrex::Real> &weights_x,
        amrex::Vector<amrex::Real> &weights_y,
        amrex::Vector<amrex::Real> &weights_z,
        amrex::GpuArray<int, 3> const &ind_min,
        amrex::GpuArray<int, 3> const &coarsen_ratio,
        amrex::GpuArray<int, 3> const &stag_fine_src,
        amrex::GpuArray<int, 3> const & stag_crse_des);

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
