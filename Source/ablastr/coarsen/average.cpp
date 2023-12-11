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
        amrex::IntVect const ixtype_src = mf_src.boxArray().ixType().toIntVect();
        amrex::IntVect const ixtype_dst = mf_dst.boxArray().ixType().toIntVect();

        // Auxiliary integer arrays (always 3D)
        auto stag_src = amrex::GpuArray<int,3>{0,0,0}; // staggering of source fine MultiFab
        auto stag_des = amrex::GpuArray<int,3>{0,0,0}; // staggering of destination coarse MultiFab
        auto crse_ratio = amrex::GpuArray<int,3>{1,1,1}; // coarsening ratio
        for (int i=0; i<AMREX_SPACEDIM; ++i)
        {
            stag_src[i] = ixtype_src[i];
            stag_des[i] = ixtype_dst[i];
            crse_ratio[i] = crse_ratio[i];
        }

        amrex::GpuArray<int,3> idx_min_src, np_src;
        amrex::GpuArray<bool,3> use_half_weight;
        amrex::Real crx_cry_crz_inv;

        CalculateCoarseningData(
            idx_min_src, np_src, use_half_weight, crx_cry_crz_inv, stag_src, stag_des, crse_ratio);

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
                        arr_src, idx_min_src, np_src, use_half_weight, crx_cry_crz_inv, crse_ratio, i, j, k, n);
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
