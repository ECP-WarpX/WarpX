/* Copyright 2022 Edoardo Zoni, Remi Lehe, Prabhat Kumar
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "average.H"

#include "ablastr/utils/TextMsg.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace ablastr::coarsen::average
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    Interp (
        amrex::Array4<amrex::Real const> const &arr_src,
        amrex::GpuArray<int, 3> const &sf,
        amrex::GpuArray<int, 3> const &sc,
        amrex::GpuArray<int, 3> const &cr,
        int const i, int const j, int const k, int const comp)
    {
        // Number of points and starting indices of source array (fine)
        amrex::GpuArray<int, 3> np, idx_min;
        amrex::GpuArray<bool, 3> use_half;
        amrex::Real crxyz_inv;

        CalculateCoarseningData(idx_min, np, use_half, crxyz_inv, sf, sc, cr);

        return InterpWithCoarseningData(arr_src, idx_min, np, use_half, crxyz_inv, cr, i, j, k, comp);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void
    CalculateCoarseningData (
        amrex::GpuArray<int, 3> &src_fine_index_min,
        amrex::GpuArray<int, 3> &num_src_points,
        amrex::GpuArray<bool, 3> &use_half_weight, // Use int instead bool?
        amrex::Real &crx_cry_crz_inv,
        amrex::GpuArray<int, 3> const &stag_src_fine,
        amrex::GpuArray<int, 3> const &stag_des_crse,
        amrex::GpuArray<int, 3> const &crse_ratio)
    {
        using namespace amrex::literals;
        for (int l = 0; l < 3; ++l) {
            int two_times_index_min = -crse_ratio[l]*stag_des_crse[l] + stag_src_fine[l] - 1;
            if ((two_times_index_min % 2) == 0) {
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
        amrex::GpuArray<bool, 3> const &use_half,
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
            amrex::Real kfactor = (use_half[2] && (kref == 0 || kref == numz-1)) ? 0.5_rt : 1.0_rt;

            for (int jref = 0; jref < numy; ++jref) {
                int jj = jjmin + jref;
                amrex::Real jfactor = (use_half[1] && (jref == 0 || jref == numy-1)) ? 0.5_rt : 1.0_rt;

                for (int iref = 0; iref < numx; ++iref) {
                    int ii = iimin + iref;
                    amrex::Real ifactor = (use_half[0] && (iref == 0 || iref == numx-1)) ? 0.5_rt : 1.0_rt;

                    c += ifactor * jfactor * kfactor * arr_src_safe(ii, jj, kk, comp);
                }
            }
        }

        c *= crx_cry_crz_inv;

        return c;
    }

    void
    Loop (
        amrex::MultiFab & mf_dst,
        amrex::MultiFab const & mf_src,
        int const ncomp,
        amrex::IntVect const ngrow,
        amrex::IntVect const crse_ratio
    )
    {
        // Staggering of source fine MultiFab and destination coarse MultiFab
        amrex::IntVect const stag_src = mf_src.boxArray().ixType().toIntVect();
        amrex::IntVect const stag_dst = mf_dst.boxArray().ixType().toIntVect();

        // Auxiliary integer arrays (always 3D)
        auto sf = amrex::GpuArray<int,3>{0,0,0}; // staggering of source fine MultiFab
        auto sc = amrex::GpuArray<int,3>{0,0,0}; // staggering of destination coarse MultiFab
        auto cr = amrex::GpuArray<int,3>{1,1,1}; // coarsening ratio
        for (int i=0; i<AMREX_SPACEDIM; ++i)
        {
             sf[i] = stag_src[i];
             sc[i] = stag_dst[i];
             cr[i] = crse_ratio[i];
        }

        amrex::GpuArray<int,3> idx_min, np;
        amrex::GpuArray<bool,3> use_half;
        amrex::Real crxyz_inv;

        CalculateCoarseningData(idx_min, np, use_half, crxyz_inv, sf, sc, cr);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        // Loop over boxes (or tiles if not on GPU)
        for (amrex::MFIter mfi(mf_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Tiles defined at the coarse level
            amrex::Box const & bx = mfi.growntilebox(ngrow);
            amrex::Array4<amrex::Real> const &arr_dst = mf_dst.array(mfi);
            amrex::Array4<amrex::Real const> const &arr_src = mf_src.const_array(mfi);
            amrex::ParallelFor(bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                    arr_dst(i, j, k, n) = InterpWithCoarseningData(
                        arr_src, idx_min, np, use_half, crxyz_inv, cr, i, j, k, n);
                });
        }
    }

    void
    Coarsen (
        amrex::MultiFab & mf_dst,
        amrex::MultiFab const & mf_src,
        amrex::IntVect const crse_ratio
    )
    {
        BL_PROFILE("ablastr::coarsen::Coarsen()");

        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(mf_src.ixType() == mf_dst.ixType(),
                                           "source MultiFab and destination MultiFab have different IndexType");

        // Number of guard cells to fill on coarse patch and number of components
        const amrex::IntVect ngrow = (mf_src.nGrowVect() + crse_ratio-1) / crse_ratio; // round up int division
        const int ncomp = mf_src.nComp();

        Loop(mf_dst, mf_src, ncomp, ngrow, crse_ratio);
    }

} // namespace ablastr::coarsen::average
