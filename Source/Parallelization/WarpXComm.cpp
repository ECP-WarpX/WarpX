/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
#   include "BoundaryConditions/PML_RZ.H"
#endif
#include "Filter/BilinearFilter.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpXComm_K.H"
#include "WarpXSumGuardCells.H"

#include <ablastr/coarsen/average.H>
#include <ablastr/utils/Communication.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MakeType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

using namespace amrex;

void
WarpX::UpdateAuxilaryData ()
{
    WARPX_PROFILE("WarpX::UpdateAuxilaryData()");

    if (Bfield_aux[0][0]->ixType() == Bfield_fp[0][0]->ixType()) {
        UpdateAuxilaryDataSameType();
    } else {
        UpdateAuxilaryDataStagToNodal();
    }
}

void
WarpX::UpdateAuxilaryDataStagToNodal ()
{
#ifndef WARPX_USE_PSATD
    if (electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE( false,
            "WarpX::UpdateAuxilaryDataStagToNodal: PSATD solver requires "
            "WarpX build with spectral solver support.");
    }
#endif

    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>> const & Bmf = WarpX::fft_do_time_averaging ?
                                                                                Bfield_avg_fp : Bfield_fp;
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>> const & Emf = WarpX::fft_do_time_averaging ?
                                                                                Efield_avg_fp : Efield_fp;

    const amrex::IntVect& Bx_stag = Bmf[0][0]->ixType().toIntVect();
    const amrex::IntVect& By_stag = Bmf[0][1]->ixType().toIntVect();
    const amrex::IntVect& Bz_stag = Bmf[0][2]->ixType().toIntVect();

    const amrex::IntVect& Ex_stag = Emf[0][0]->ixType().toIntVect();
    const amrex::IntVect& Ey_stag = Emf[0][1]->ixType().toIntVect();
    const amrex::IntVect& Ez_stag = Emf[0][2]->ixType().toIntVect();

    // Destination MultiFab (aux) always has nodal index type when this function is called
    const amrex::IntVect& dst_stag = amrex::IntVect::TheNodeVector();

    // For level 0, we only need to do the average.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*Bfield_aux[0][0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& bx_aux = Bfield_aux[0][0]->array(mfi);
        Array4<Real> const& by_aux = Bfield_aux[0][1]->array(mfi);
        Array4<Real> const& bz_aux = Bfield_aux[0][2]->array(mfi);
        Array4<Real const> const& bx_fp = Bmf[0][0]->const_array(mfi);
        Array4<Real const> const& by_fp = Bmf[0][1]->const_array(mfi);
        Array4<Real const> const& bz_fp = Bmf[0][2]->const_array(mfi);

        Array4<Real> const& ex_aux = Efield_aux[0][0]->array(mfi);
        Array4<Real> const& ey_aux = Efield_aux[0][1]->array(mfi);
        Array4<Real> const& ez_aux = Efield_aux[0][2]->array(mfi);
        Array4<Real const> const& ex_fp = Emf[0][0]->const_array(mfi);
        Array4<Real const> const& ey_fp = Emf[0][1]->const_array(mfi);
        Array4<Real const> const& ez_fp = Emf[0][2]->const_array(mfi);

        // Loop includes ghost cells (`growntilebox`)
        // (input arrays will be padded with zeros beyond ghost cells
        // for out-of-bound accesses due to large-stencil operations)
        const Box bx = mfi.growntilebox();

        // Order of finite-order centering of fields
        const int fg_nox = WarpX::field_centering_nox;
        const int fg_noy = WarpX::field_centering_noy;
        const int fg_noz = WarpX::field_centering_noz;

        // Device vectors of stencil coefficients used for finite-order centering of fields
        amrex::Real const * stencil_coeffs_x = WarpX::device_field_centering_stencil_coeffs_x.data();
        amrex::Real const * stencil_coeffs_y = WarpX::device_field_centering_stencil_coeffs_y.data();
        amrex::Real const * stencil_coeffs_z = WarpX::device_field_centering_stencil_coeffs_z.data();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
        {
            warpx_interp(j, k, l, bx_aux, bx_fp, dst_stag, Bx_stag, fg_nox, fg_noy, fg_noz,
                         stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

            warpx_interp(j, k, l, by_aux, by_fp, dst_stag, By_stag, fg_nox, fg_noy, fg_noz,
                         stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

            warpx_interp(j, k, l, bz_aux, bz_fp, dst_stag, Bz_stag, fg_nox, fg_noy, fg_noz,
                         stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

            warpx_interp(j, k, l, ex_aux, ex_fp, dst_stag, Ex_stag, fg_nox, fg_noy, fg_noz,
                         stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

            warpx_interp(j, k, l, ey_aux, ey_fp, dst_stag, Ey_stag, fg_nox, fg_noy, fg_noz,
                         stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

            warpx_interp(j, k, l, ez_aux, ez_fp, dst_stag, Ez_stag, fg_nox, fg_noy, fg_noz,
                         stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
        });
    }

    // NOTE: high-order interpolation is not implemented for mesh refinement
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        BoxArray const& nba = Bfield_aux[lev][0]->boxArray();
        BoxArray const& cnba = amrex::coarsen(nba,2);
        DistributionMapping const& dm = Bfield_aux[lev][0]->DistributionMap();
        amrex::Periodicity const& cperiod = Geom(lev-1).periodicity();

        // Bfield
        {
            Array<std::unique_ptr<MultiFab>,3> Btmp;
            if (Bfield_cax[lev][0]) {
                for (int i = 0; i < 3; ++i) {
                    Btmp[i] = std::make_unique<MultiFab>(
                        *Bfield_cax[lev][i], amrex::make_alias, 0, 1);
                }
            } else {
                const IntVect ngtmp = Bfield_aux[lev-1][0]->nGrowVect();
                for (int i = 0; i < 3; ++i) {
                    Btmp[i] = std::make_unique<MultiFab>(cnba, dm, 1, ngtmp);
                }
            }
            Btmp[0]->setVal(0.0);
            Btmp[1]->setVal(0.0);
            Btmp[2]->setVal(0.0);
            // ParallelCopy from coarse level
            for (int i = 0; i < 3; ++i) {
                const IntVect ng = Btmp[i]->nGrowVect();
                // Guard cells may not be up to date beyond ng_FieldGather
                const amrex::IntVect& ng_src = guard_cells.ng_FieldGather;
                // Copy Bfield_aux to Btmp, using up to ng_src (=ng_FieldGather) guard cells from
                // Bfield_aux and filling up to ng (=nGrow) guard cells in Btmp
                ablastr::utils::communication::ParallelCopy(*Btmp[i], *Bfield_aux[lev - 1][i], 0, 0, 1,
                                                            ng_src, ng, WarpX::do_single_precision_comms, cperiod);
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*Bfield_aux[lev][0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Array4<Real> const& bx_aux = Bfield_aux[lev][0]->array(mfi);
                Array4<Real> const& by_aux = Bfield_aux[lev][1]->array(mfi);
                Array4<Real> const& bz_aux = Bfield_aux[lev][2]->array(mfi);
                Array4<Real const> const& bx_fp = Bfield_fp[lev][0]->const_array(mfi);
                Array4<Real const> const& by_fp = Bfield_fp[lev][1]->const_array(mfi);
                Array4<Real const> const& bz_fp = Bfield_fp[lev][2]->const_array(mfi);
                Array4<Real const> const& bx_cp = Bfield_cp[lev][0]->const_array(mfi);
                Array4<Real const> const& by_cp = Bfield_cp[lev][1]->const_array(mfi);
                Array4<Real const> const& bz_cp = Bfield_cp[lev][2]->const_array(mfi);
                Array4<Real const> const& bx_c = Btmp[0]->const_array(mfi);
                Array4<Real const> const& by_c = Btmp[1]->const_array(mfi);
                Array4<Real const> const& bz_c = Btmp[2]->const_array(mfi);

                const Box& bx = mfi.growntilebox();
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp_nd_bfield_x(j,k,l, bx_aux, bx_fp, bx_cp, bx_c);
                    warpx_interp_nd_bfield_y(j,k,l, by_aux, by_fp, by_cp, by_c);
                    warpx_interp_nd_bfield_z(j,k,l, bz_aux, bz_fp, bz_cp, bz_c);
                });
            }
        }

        // Efield
        {
            Array<std::unique_ptr<MultiFab>,3> Etmp;
            if (Efield_cax[lev][0]) {
                for (int i = 0; i < 3; ++i) {
                    Etmp[i] = std::make_unique<MultiFab>(
                        *Efield_cax[lev][i], amrex::make_alias, 0, 1);
                }
            } else {
                const IntVect ngtmp = Efield_aux[lev-1][0]->nGrowVect();
                for (int i = 0; i < 3; ++i) {
                    Etmp[i] = std::make_unique<MultiFab>(
                        cnba, dm, 1, ngtmp);
                }
            }
            Etmp[0]->setVal(0.0);
            Etmp[1]->setVal(0.0);
            Etmp[2]->setVal(0.0);
            // ParallelCopy from coarse level
            for (int i = 0; i < 3; ++i) {
                const IntVect ng = Etmp[i]->nGrowVect();
                // Guard cells may not be up to date beyond ng_FieldGather
                const amrex::IntVect& ng_src = guard_cells.ng_FieldGather;
                // Copy Efield_aux to Etmp, using up to ng_src (=ng_FieldGather) guard cells from
                // Efield_aux and filling up to ng (=nGrow) guard cells in Etmp
                ablastr::utils::communication::ParallelCopy(*Etmp[i], *Efield_aux[lev - 1][i], 0, 0, 1,
                                                            ng_src, ng, WarpX::do_single_precision_comms, cperiod);
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*Efield_aux[lev][0]); mfi.isValid(); ++mfi)
            {
                Array4<Real> const& ex_aux = Efield_aux[lev][0]->array(mfi);
                Array4<Real> const& ey_aux = Efield_aux[lev][1]->array(mfi);
                Array4<Real> const& ez_aux = Efield_aux[lev][2]->array(mfi);
                Array4<Real const> const& ex_fp = Efield_fp[lev][0]->const_array(mfi);
                Array4<Real const> const& ey_fp = Efield_fp[lev][1]->const_array(mfi);
                Array4<Real const> const& ez_fp = Efield_fp[lev][2]->const_array(mfi);
                Array4<Real const> const& ex_cp = Efield_cp[lev][0]->const_array(mfi);
                Array4<Real const> const& ey_cp = Efield_cp[lev][1]->const_array(mfi);
                Array4<Real const> const& ez_cp = Efield_cp[lev][2]->const_array(mfi);
                Array4<Real const> const& ex_c = Etmp[0]->const_array(mfi);
                Array4<Real const> const& ey_c = Etmp[1]->const_array(mfi);
                Array4<Real const> const& ez_c = Etmp[2]->const_array(mfi);

                const Box& bx = mfi.fabbox();
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp_nd_efield_x(j,k,l, ex_aux, ex_fp, ex_cp, ex_c);
                    warpx_interp_nd_efield_y(j,k,l, ey_aux, ey_fp, ey_cp, ey_c);
                    warpx_interp_nd_efield_z(j,k,l, ez_aux, ez_fp, ez_cp, ez_c);
                });
            }
        }
    }
}

void
WarpX::UpdateAuxilaryDataSameType ()
{
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const amrex::Periodicity& crse_period = Geom(lev-1).periodicity();
        const IntVect& ng = Bfield_cp[lev][0]->nGrowVect();
        const DistributionMapping& dm = Bfield_cp[lev][0]->DistributionMap();

        // B field
        {
            MultiFab dBx(Bfield_cp[lev][0]->boxArray(), dm, Bfield_cp[lev][0]->nComp(), ng);
            MultiFab dBy(Bfield_cp[lev][1]->boxArray(), dm, Bfield_cp[lev][1]->nComp(), ng);
            MultiFab dBz(Bfield_cp[lev][2]->boxArray(), dm, Bfield_cp[lev][2]->nComp(), ng);
            dBx.setVal(0.0);
            dBy.setVal(0.0);
            dBz.setVal(0.0);

            // Guard cells may not be up to date beyond ng_FieldGather
            const amrex::IntVect& ng_src = guard_cells.ng_FieldGather;
            // Copy Bfield_aux to the dB MultiFabs, using up to ng_src (=ng_FieldGather) guard
            // cells from Bfield_aux and filling up to ng (=nGrow) guard cells in the dB MultiFabs

            ablastr::utils::communication::ParallelCopy(dBx, *Bfield_aux[lev - 1][0], 0, 0,
                                                        Bfield_aux[lev - 1][0]->nComp(), ng_src, ng, WarpX::do_single_precision_comms,
                                                        crse_period);
            ablastr::utils::communication::ParallelCopy(dBy, *Bfield_aux[lev - 1][1], 0, 0,
                                                        Bfield_aux[lev - 1][1]->nComp(), ng_src, ng, WarpX::do_single_precision_comms,
                                                        crse_period);
            ablastr::utils::communication::ParallelCopy(dBz, *Bfield_aux[lev - 1][2], 0, 0,
                                                        Bfield_aux[lev - 1][2]->nComp(), ng_src, ng, WarpX::do_single_precision_comms,
                                                        crse_period);

            if (Bfield_cax[lev][0])
            {
                MultiFab::Copy(*Bfield_cax[lev][0], dBx, 0, 0, Bfield_cax[lev][0]->nComp(), ng);
                MultiFab::Copy(*Bfield_cax[lev][1], dBy, 0, 0, Bfield_cax[lev][1]->nComp(), ng);
                MultiFab::Copy(*Bfield_cax[lev][2], dBz, 0, 0, Bfield_cax[lev][2]->nComp(), ng);
            }
            MultiFab::Subtract(dBx, *Bfield_cp[lev][0], 0, 0, Bfield_cp[lev][0]->nComp(), ng);
            MultiFab::Subtract(dBy, *Bfield_cp[lev][1], 0, 0, Bfield_cp[lev][1]->nComp(), ng);
            MultiFab::Subtract(dBz, *Bfield_cp[lev][2], 0, 0, Bfield_cp[lev][2]->nComp(), ng);

            const amrex::IntVect& refinement_ratio = refRatio(lev-1);

            const amrex::IntVect& Bx_stag = Bfield_aux[lev-1][0]->ixType().toIntVect();
            const amrex::IntVect& By_stag = Bfield_aux[lev-1][1]->ixType().toIntVect();
            const amrex::IntVect& Bz_stag = Bfield_aux[lev-1][2]->ixType().toIntVect();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*Bfield_aux[lev][0]); mfi.isValid(); ++mfi)
            {
                Array4<Real> const& bx_aux = Bfield_aux[lev][0]->array(mfi);
                Array4<Real> const& by_aux = Bfield_aux[lev][1]->array(mfi);
                Array4<Real> const& bz_aux = Bfield_aux[lev][2]->array(mfi);
                Array4<Real const> const& bx_fp = Bfield_fp[lev][0]->const_array(mfi);
                Array4<Real const> const& by_fp = Bfield_fp[lev][1]->const_array(mfi);
                Array4<Real const> const& bz_fp = Bfield_fp[lev][2]->const_array(mfi);
                Array4<Real const> const& bx_c = dBx.const_array(mfi);
                Array4<Real const> const& by_c = dBy.const_array(mfi);
                Array4<Real const> const& bz_c = dBz.const_array(mfi);

                amrex::ParallelFor(Box(bx_aux), Box(by_aux), Box(bz_aux),
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, bx_aux, bx_fp, bx_c, Bx_stag, refinement_ratio);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, by_aux, by_fp, by_c, By_stag, refinement_ratio);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, bz_aux, bz_fp, bz_c, Bz_stag, refinement_ratio);
                });
            }
        }

        // E field
        {
            MultiFab dEx(Efield_cp[lev][0]->boxArray(), dm, Efield_cp[lev][0]->nComp(), ng);
            MultiFab dEy(Efield_cp[lev][1]->boxArray(), dm, Efield_cp[lev][1]->nComp(), ng);
            MultiFab dEz(Efield_cp[lev][2]->boxArray(), dm, Efield_cp[lev][2]->nComp(), ng);
            dEx.setVal(0.0);
            dEy.setVal(0.0);
            dEz.setVal(0.0);

            // Guard cells may not be up to date beyond ng_FieldGather
            const amrex::IntVect& ng_src = guard_cells.ng_FieldGather;
            // Copy Efield_aux to the dE MultiFabs, using up to ng_src (=ng_FieldGather) guard
            // cells from Efield_aux and filling up to ng (=nGrow) guard cells in the dE MultiFabs
            ablastr::utils::communication::ParallelCopy(dEx, *Efield_aux[lev - 1][0], 0, 0,
                                                        Efield_aux[lev - 1][0]->nComp(), ng_src, ng,
                                                        WarpX::do_single_precision_comms,
                                                        crse_period);
            ablastr::utils::communication::ParallelCopy(dEy, *Efield_aux[lev - 1][1], 0, 0,
                                                        Efield_aux[lev - 1][1]->nComp(), ng_src, ng,
                                                        WarpX::do_single_precision_comms,
                                                        crse_period);
            ablastr::utils::communication::ParallelCopy(dEz, *Efield_aux[lev - 1][2], 0, 0,
                                                        Efield_aux[lev - 1][2]->nComp(), ng_src, ng,
                                                        WarpX::do_single_precision_comms,
                                                        crse_period);

            if (Efield_cax[lev][0])
            {
                MultiFab::Copy(*Efield_cax[lev][0], dEx, 0, 0, Efield_cax[lev][0]->nComp(), ng);
                MultiFab::Copy(*Efield_cax[lev][1], dEy, 0, 0, Efield_cax[lev][1]->nComp(), ng);
                MultiFab::Copy(*Efield_cax[lev][2], dEz, 0, 0, Efield_cax[lev][2]->nComp(), ng);
            }
            MultiFab::Subtract(dEx, *Efield_cp[lev][0], 0, 0, Efield_cp[lev][0]->nComp(), ng);
            MultiFab::Subtract(dEy, *Efield_cp[lev][1], 0, 0, Efield_cp[lev][1]->nComp(), ng);
            MultiFab::Subtract(dEz, *Efield_cp[lev][2], 0, 0, Efield_cp[lev][2]->nComp(), ng);

            const amrex::IntVect& refinement_ratio = refRatio(lev-1);

            const amrex::IntVect& Ex_stag = Efield_aux[lev-1][0]->ixType().toIntVect();
            const amrex::IntVect& Ey_stag = Efield_aux[lev-1][1]->ixType().toIntVect();
            const amrex::IntVect& Ez_stag = Efield_aux[lev-1][2]->ixType().toIntVect();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*Efield_aux[lev][0]); mfi.isValid(); ++mfi)
            {
                Array4<Real> const& ex_aux = Efield_aux[lev][0]->array(mfi);
                Array4<Real> const& ey_aux = Efield_aux[lev][1]->array(mfi);
                Array4<Real> const& ez_aux = Efield_aux[lev][2]->array(mfi);
                Array4<Real const> const& ex_fp = Efield_fp[lev][0]->const_array(mfi);
                Array4<Real const> const& ey_fp = Efield_fp[lev][1]->const_array(mfi);
                Array4<Real const> const& ez_fp = Efield_fp[lev][2]->const_array(mfi);
                Array4<Real const> const& ex_c = dEx.const_array(mfi);
                Array4<Real const> const& ey_c = dEy.const_array(mfi);
                Array4<Real const> const& ez_c = dEz.const_array(mfi);

                amrex::ParallelFor(Box(ex_aux), Box(ey_aux), Box(ez_aux),
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, ex_aux, ex_fp, ex_c, Ex_stag, refinement_ratio);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, ey_aux, ey_fp, ey_c, Ey_stag, refinement_ratio);
                },
                [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
                {
                    warpx_interp(j, k, l, ez_aux, ez_fp, ez_c, Ez_stag, refinement_ratio);
                });
            }
        }
    }
}

void WarpX::UpdateCurrentNodalToStag (amrex::MultiFab& dst, amrex::MultiFab const& src)
{
    // If source and destination MultiFabs have the same index type, a simple copy is enough
    // (for example, this happens with the current along y in 2D, which is always fully nodal)
    if (dst.ixType() == src.ixType())
    {
        amrex::MultiFab::Copy(dst, src, 0, 0, dst.nComp(), dst.nGrowVect());
        return;
    }

    amrex::IntVect const& dst_stag = dst.ixType().toIntVect();

    // Source MultiFab always has nodal index type when this function is called
    amrex::IntVect const& src_stag = amrex::IntVect::TheNodeVector();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Loop over full box including ghost cells
        // (input arrays will be padded with zeros beyond ghost cells
        // for out-of-bound accesses due to large-stencil operations)
        const Box bx = mfi.growntilebox();

        amrex::Array4<amrex::Real const> const& src_arr = src.const_array(mfi);
        amrex::Array4<amrex::Real>       const& dst_arr = dst.array(mfi);

        // Order of finite-order centering of currents
        const int cc_nox = WarpX::current_centering_nox;
        const int cc_noy = WarpX::current_centering_noy;
        const int cc_noz = WarpX::current_centering_noz;

        // Device vectors of stencil coefficients used for finite-order centering of currents
        amrex::Real const * stencil_coeffs_x = WarpX::device_current_centering_stencil_coeffs_x.data();
        amrex::Real const * stencil_coeffs_y = WarpX::device_current_centering_stencil_coeffs_y.data();
        amrex::Real const * stencil_coeffs_z = WarpX::device_current_centering_stencil_coeffs_z.data();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
        {
            warpx_interp(j, k, l, dst_arr, src_arr, dst_stag, src_stag, cc_nox, cc_noy, cc_noz,
                         stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
        });
    }
}

void
WarpX::FillBoundaryB (IntVect ng, std::optional<bool> nodal_sync)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryB(lev, ng, nodal_sync);
    }
}

void
WarpX::FillBoundaryE (IntVect ng, std::optional<bool> nodal_sync)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryE(lev, ng, nodal_sync);
    }
}

void
WarpX::FillBoundaryF (IntVect ng, std::optional<bool> nodal_sync)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryF(lev, ng, nodal_sync);
    }
}

void
WarpX::FillBoundaryG (IntVect ng, std::optional<bool> nodal_sync)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryG(lev, ng, nodal_sync);
    }
}

void
WarpX::FillBoundaryB_avg (IntVect ng)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryB_avg(lev, ng);
    }
}

void
WarpX::FillBoundaryE_avg (IntVect ng)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryE_avg(lev, ng);
    }
}


void
WarpX::FillBoundaryE (int lev, IntVect ng, std::optional<bool> nodal_sync)
{
    FillBoundaryE(lev, PatchType::fine, ng, nodal_sync);
    if (lev > 0) FillBoundaryE(lev, PatchType::coarse, ng, nodal_sync);
}

void
WarpX::FillBoundaryE (const int lev, const PatchType patch_type, const amrex::IntVect ng, std::optional<bool> nodal_sync)
{
    std::array<amrex::MultiFab*,3> mf;
    amrex::Periodicity period;

    if (patch_type == PatchType::fine)
    {
        mf     = {Efield_fp[lev][0].get(), Efield_fp[lev][1].get(), Efield_fp[lev][2].get()};
        period = Geom(lev).periodicity();
    }
    else // coarse patch
    {
        mf     = {Efield_cp[lev][0].get(), Efield_cp[lev][1].get(), Efield_cp[lev][2].get()};
        period = Geom(lev-1).periodicity();
    }

    // Exchange data between valid domain and PML
    // Fill guard cells in PML
    if (do_pml)
    {
        if (pml[lev] && pml[lev]->ok())
        {
            const std::array<amrex::MultiFab*,3> mf_pml =
                (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();

            pml[lev]->Exchange(mf_pml, mf, patch_type, do_pml_in_domain);
            pml[lev]->FillBoundaryE(patch_type);
        }

#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
        if (pml_rz[lev])
        {
            pml_rz[lev]->FillBoundaryE(patch_type);
        }
#endif
    }

    // Fill guard cells in valid domain
    for (int i = 0; i < 3; ++i)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            ng <= mf[i]->nGrowVect(),
            "Error: in FillBoundaryE, requested more guard cells than allocated");

        const amrex::IntVect nghost = (safe_guard_cells) ? mf[i]->nGrowVect() : ng;
        ablastr::utils::communication::FillBoundary(*mf[i], nghost, WarpX::do_single_precision_comms, period, nodal_sync);
    }
}

void
WarpX::FillBoundaryB (int lev, IntVect ng, std::optional<bool> nodal_sync)
{
    FillBoundaryB(lev, PatchType::fine, ng, nodal_sync);
    if (lev > 0) FillBoundaryB(lev, PatchType::coarse, ng, nodal_sync);
}

void
WarpX::FillBoundaryB (const int lev, const PatchType patch_type, const amrex::IntVect ng, std::optional<bool> nodal_sync)
{
    std::array<amrex::MultiFab*,3> mf;
    amrex::Periodicity period;

    if (patch_type == PatchType::fine)
    {
        mf     = {Bfield_fp[lev][0].get(), Bfield_fp[lev][1].get(), Bfield_fp[lev][2].get()};
        period = Geom(lev).periodicity();
    }
    else // coarse patch
    {
        mf     = {Bfield_cp[lev][0].get(), Bfield_cp[lev][1].get(), Bfield_cp[lev][2].get()};
        period = Geom(lev-1).periodicity();
    }

    // Exchange data between valid domain and PML
    // Fill guard cells in PML
    if (do_pml)
    {
        if (pml[lev] && pml[lev]->ok())
        {
            const std::array<amrex::MultiFab*,3> mf_pml =
                (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();

            pml[lev]->Exchange(mf_pml, mf, patch_type, do_pml_in_domain);
            pml[lev]->FillBoundaryB(patch_type);
        }

#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
        if (pml_rz[lev])
        {
            pml_rz[lev]->FillBoundaryB(patch_type);
        }
#endif
    }

    // Fill guard cells in valid domain
    for (int i = 0; i < 3; ++i)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            ng <= mf[i]->nGrowVect(),
            "Error: in FillBoundaryB, requested more guard cells than allocated");

        const amrex::IntVect nghost = (safe_guard_cells) ? mf[i]->nGrowVect() : ng;
        ablastr::utils::communication::FillBoundary(*mf[i], nghost, WarpX::do_single_precision_comms, period, nodal_sync);
    }
}

void
WarpX::FillBoundaryE_avg(int lev, IntVect ng)
{
    FillBoundaryE_avg(lev, PatchType::fine, ng);
    if (lev > 0) FillBoundaryE_avg(lev, PatchType::coarse, ng);
}

void
WarpX::FillBoundaryE_avg (int lev, PatchType patch_type, IntVect ng)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev]->ok())
         {
            WARPX_ABORT_WITH_MESSAGE("Averaged Galilean PSATD with PML is not yet implemented");
         }

        const amrex::Periodicity& period = Geom(lev).periodicity();
        if ( safe_guard_cells ){
            const Vector<MultiFab*> mf{Efield_avg_fp[lev][0].get(),Efield_avg_fp[lev][1].get(),Efield_avg_fp[lev][2].get()};
            ablastr::utils::communication::FillBoundary(mf, WarpX::do_single_precision_comms, period);
        } else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Efield_avg_fp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryE_avg, requested more guard cells than allocated");
            ablastr::utils::communication::FillBoundary(*Efield_avg_fp[lev][0], ng, WarpX::do_single_precision_comms, period);
            ablastr::utils::communication::FillBoundary(*Efield_avg_fp[lev][1], ng, WarpX::do_single_precision_comms, period);
            ablastr::utils::communication::FillBoundary(*Efield_avg_fp[lev][2], ng, WarpX::do_single_precision_comms, period);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
         {
            WARPX_ABORT_WITH_MESSAGE("Averaged Galilean PSATD with PML is not yet implemented");
         }

        const amrex::Periodicity& cperiod = Geom(lev-1).periodicity();
        if ( safe_guard_cells ) {
            const Vector<MultiFab*> mf{Efield_avg_cp[lev][0].get(),Efield_avg_cp[lev][1].get(),Efield_avg_cp[lev][2].get()};
            ablastr::utils::communication::FillBoundary(mf, WarpX::do_single_precision_comms, cperiod);

        } else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Efield_avg_cp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryE, requested more guard cells than allocated");
            ablastr::utils::communication::FillBoundary(*Efield_avg_cp[lev][0], ng, WarpX::do_single_precision_comms, cperiod);
            ablastr::utils::communication::FillBoundary(*Efield_avg_cp[lev][1], ng, WarpX::do_single_precision_comms, cperiod);
            ablastr::utils::communication::FillBoundary(*Efield_avg_cp[lev][2], ng, WarpX::do_single_precision_comms, cperiod);
        }
    }
}


void
WarpX::FillBoundaryB_avg (int lev, IntVect ng)
{
    FillBoundaryB_avg(lev, PatchType::fine, ng);
    if (lev > 0) FillBoundaryB_avg(lev, PatchType::coarse, ng);
}

void
WarpX::FillBoundaryB_avg (int lev, PatchType patch_type, IntVect ng)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev]->ok())
          {
            WARPX_ABORT_WITH_MESSAGE("Averaged Galilean PSATD with PML is not yet implemented");
          }
        const amrex::Periodicity& period = Geom(lev).periodicity();
        if ( safe_guard_cells ) {
            const Vector<MultiFab*> mf{Bfield_avg_fp[lev][0].get(),Bfield_avg_fp[lev][1].get(),Bfield_avg_fp[lev][2].get()};
            ablastr::utils::communication::FillBoundary(mf, WarpX::do_single_precision_comms, period);
        } else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Bfield_fp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryB, requested more guard cells than allocated");
            ablastr::utils::communication::FillBoundary(*Bfield_avg_fp[lev][0], ng, WarpX::do_single_precision_comms, period);
            ablastr::utils::communication::FillBoundary(*Bfield_avg_fp[lev][1], ng, WarpX::do_single_precision_comms, period);
            ablastr::utils::communication::FillBoundary(*Bfield_avg_fp[lev][2], ng, WarpX::do_single_precision_comms, period);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
          {
            WARPX_ABORT_WITH_MESSAGE("Averaged Galilean PSATD with PML is not yet implemented");
          }

        const amrex::Periodicity& cperiod = Geom(lev-1).periodicity();
        if ( safe_guard_cells ){
            const Vector<MultiFab*> mf{Bfield_avg_cp[lev][0].get(),Bfield_avg_cp[lev][1].get(),Bfield_avg_cp[lev][2].get()};
            ablastr::utils::communication::FillBoundary(mf, WarpX::do_single_precision_comms, cperiod);
        } else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Bfield_avg_cp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryB_avg, requested more guard cells than allocated");
            ablastr::utils::communication::FillBoundary(*Bfield_avg_cp[lev][0], ng, WarpX::do_single_precision_comms, cperiod);
            ablastr::utils::communication::FillBoundary(*Bfield_avg_cp[lev][1], ng, WarpX::do_single_precision_comms, cperiod);
            ablastr::utils::communication::FillBoundary(*Bfield_avg_cp[lev][2], ng, WarpX::do_single_precision_comms, cperiod);
        }
    }
}

void
WarpX::FillBoundaryF (int lev, IntVect ng, std::optional<bool> nodal_sync)
{
    FillBoundaryF(lev, PatchType::fine, ng, nodal_sync);
    if (lev > 0) FillBoundaryF(lev, PatchType::coarse, ng, nodal_sync);
}

void
WarpX::FillBoundaryF (int lev, PatchType patch_type, IntVect ng, std::optional<bool> nodal_sync)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev] && pml[lev]->ok())
        {
            if (F_fp[lev]) pml[lev]->ExchangeF(patch_type, F_fp[lev].get(), do_pml_in_domain);
            pml[lev]->FillBoundaryF(patch_type);
        }

        if (F_fp[lev])
        {
            const amrex::Periodicity& period = Geom(lev).periodicity();
            const amrex::IntVect& nghost = (safe_guard_cells) ? F_fp[lev]->nGrowVect() : ng;
            ablastr::utils::communication::FillBoundary(*F_fp[lev], nghost, WarpX::do_single_precision_comms, period, nodal_sync);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev] && pml[lev]->ok())
        {
            if (F_cp[lev]) pml[lev]->ExchangeF(patch_type, F_cp[lev].get(), do_pml_in_domain);
            pml[lev]->FillBoundaryF(patch_type);
        }

        if (F_cp[lev])
        {
            const amrex::Periodicity& period = Geom(lev-1).periodicity();
            const amrex::IntVect& nghost = (safe_guard_cells) ? F_cp[lev]->nGrowVect() : ng;
            ablastr::utils::communication::FillBoundary(*F_cp[lev], nghost, WarpX::do_single_precision_comms, period, nodal_sync);
        }
    }
}

void WarpX::FillBoundaryG (int lev, IntVect ng, std::optional<bool> nodal_sync)
{
    FillBoundaryG(lev, PatchType::fine, ng, nodal_sync);

    if (lev > 0)
    {
        FillBoundaryG(lev, PatchType::coarse, ng, nodal_sync);
    }
}

void WarpX::FillBoundaryG (int lev, PatchType patch_type, IntVect ng, std::optional<bool> nodal_sync)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev] && pml[lev]->ok())
        {
            if (G_fp[lev]) pml[lev]->ExchangeG(patch_type, G_fp[lev].get(), do_pml_in_domain);
            pml[lev]->FillBoundaryG(patch_type);
        }

        if (G_fp[lev])
        {
            const amrex::Periodicity& period = Geom(lev).periodicity();
            const amrex::IntVect& nghost = (safe_guard_cells) ? G_fp[lev]->nGrowVect() : ng;
            ablastr::utils::communication::FillBoundary(*G_fp[lev], nghost, WarpX::do_single_precision_comms, period, nodal_sync);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev] && pml[lev]->ok())
        {
            if (G_cp[lev]) pml[lev]->ExchangeG(patch_type, G_cp[lev].get(), do_pml_in_domain);
            pml[lev]->FillBoundaryG(patch_type);
        }

        if (G_cp[lev])
        {
            const amrex::Periodicity& period = Geom(lev-1).periodicity();
            const amrex::IntVect& nghost = (safe_guard_cells) ? G_cp[lev]->nGrowVect() : ng;
            ablastr::utils::communication::FillBoundary(*G_cp[lev], nghost, WarpX::do_single_precision_comms, period, nodal_sync);
        }
    }
}

void
WarpX::FillBoundaryAux (IntVect ng)
{
    for (int lev = 0; lev <= finest_level-1; ++lev)
    {
        FillBoundaryAux(lev, ng);
    }
}

void
WarpX::FillBoundaryAux (int lev, IntVect ng)
{
    const amrex::Periodicity& period = Geom(lev).periodicity();
    ablastr::utils::communication::FillBoundary(*Efield_aux[lev][0], ng, WarpX::do_single_precision_comms, period);
    ablastr::utils::communication::FillBoundary(*Efield_aux[lev][1], ng, WarpX::do_single_precision_comms, period);
    ablastr::utils::communication::FillBoundary(*Efield_aux[lev][2], ng, WarpX::do_single_precision_comms, period);
    ablastr::utils::communication::FillBoundary(*Bfield_aux[lev][0], ng, WarpX::do_single_precision_comms, period);
    ablastr::utils::communication::FillBoundary(*Bfield_aux[lev][1], ng, WarpX::do_single_precision_comms, period);
    ablastr::utils::communication::FillBoundary(*Bfield_aux[lev][2], ng, WarpX::do_single_precision_comms, period);
}

void
WarpX::SyncCurrent (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp)
{
    WARPX_PROFILE("WarpX::SyncCurrent()");

    // If warpx.do_current_centering = 1, center currents from nodal grid to staggered grid
    if (WarpX::do_current_centering)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(finest_level <= 1,
                                         "warpx.do_current_centering=1 not supported with more than one fine levels");
        for (int lev = 0; lev <= finest_level; lev++)
        {
            WarpX::UpdateCurrentNodalToStag(*J_fp[lev][0], *current_fp_nodal[lev][0]);
            WarpX::UpdateCurrentNodalToStag(*J_fp[lev][1], *current_fp_nodal[lev][1]);
            WarpX::UpdateCurrentNodalToStag(*J_fp[lev][2], *current_fp_nodal[lev][2]);
        }
    }

    // If there is a single level, we apply the filter on the fp data and
    // then call SumBoundary that adds data from different boxes. This needs
    // to be done because a particle near a box boundary may deposit current
    // at a given (i,j,k) that is on the edge or in the ghost region of that
    // box while at the same time that (i,j,k) is also in the valid region
    // of another box. After SumBoundary, the result is as if there is only
    // a single box on a single process. Also note we need to call
    // SumBoundary even if there is only a single process, because a process
    // may have multiple boxes. Furthermore, even if there is only a single
    // box on a single process, SumBoundary should also be called if there
    // are periodic boundaries. So we always call SumBounary even if it
    // might be a no-op in some cases, because the function does not perform
    // any communication if not necessary.
    //
    // When there are multiple levels, we need to send data from fine levels
    // to coarse levels. In the implementation below, we loop over levels
    // from the finest to the coarsest. On each level, filtering and
    // SumBoundary are done as the last two things. So the communication
    // data on the sender side are always unfiltered and unsummed. The
    // receivers are responsible for filtering that is dependent on the grid
    // resolution. On the finest level, we coarsen the unsummed fp data onto
    // cp grids, which are the coarsened version of the fp grids with the
    // same DistributionMapping. Then on the level below, We use ParallelAdd
    // to add the finer level's cp data onto the current level's fp
    // data. After that, we apply filter and SumBoundary to the finer
    // level's cp MultiFab. At this time, the finer level's fp and cp data
    // have all been properly filtered and summed. For the current level, if
    // there are levels below this, we need to process this level's cp data
    // just like we have done for the finer level. The iteration continues
    // until we reache level 0. There are however two additional
    // complications.
    //
    // The first complication is that simply calling ParallelAdd to add the
    // finer level's cp data to the current level's fp data does not
    // work. Suppose there are multiple boxes on the current level (or just
    // a single box with periodic bounaries). A given (i,j,k) can be present
    // in more than one box for nodal data in AMReX.
    // At the time of calling ParallelAdd, the current
    // level's fp data have not been summed. Because of how ParallelAdd
    // works, all boxes with that (i,j,k) will receive the data. So there is
    // a double counting issue of those data points existing on multiple boxes. Note
    // that at this time, the current level's fp data have not been summed
    // and we will call SumBoundary on the fp data. That would overcount the
    // finer level's cp data. So we fix this issue by creating a temporary
    // MultiFab to receive the finer level's cp data. We also create a mask
    // that can mark only one instance of the data as the owner if there are
    // overlapping points among boxes. Using the mask, we can add the
    // temporary MultiFab's data to the fp MultiFab only if the source owns
    // the data.
    //
    // The other complication is there might be a current buffer depending a
    // runtime parameter. The current buffer data, if they exist, need to be
    // communicated to the coarser level's fp MultiFab just like the cp
    // data. A simple approach would be to call another ParallelAdd in
    // additional the ParallelAdd for the cp data. But we like to minimize
    // parallel communication. So we add the cp data to the current buffer
    // MultiFab and use the latter as the source of ParallelAdd
    // communication to the coarser level. Note that we can do it this way
    // but not the other way around of adding the current buffer data to the
    // cp MultiFab because we still need to use the original cp data whereas
    // the buffer data are no longer needed once they have been sent to the
    // coarser level. So there are two cases. If there is no current buffer,
    // the cp MultiFab is the source of communication. If there is a current
    // buffer, the buffer MultiFab is the source instead. In the
    // implementation below, we use an alias MultiFab to manage this.

    std::unique_ptr<MultiFab> mf_comm; // for communication between levels
    for (int idim = 0; idim < 3; ++idim)
    {
        for (int lev = finest_level; lev >= 0; --lev)
        {
            const int ncomp = J_fp[lev][idim]->nComp();
            auto const& period = Geom(lev).periodicity();

            if (lev < finest_level)
            {
                // On a coarse level, the data in mf_comm comes from the
                // coarse patch of the fine level. They are unfiltered and uncommunicated.
                // We need to add it to the fine patch of the current level.
                MultiFab fine_lev_cp(J_fp[lev][idim]->boxArray(),
                                     J_fp[lev][idim]->DistributionMap(),
                                     ncomp, 0);
                fine_lev_cp.setVal(0.0);
                fine_lev_cp.ParallelAdd(*mf_comm, 0, 0, ncomp, mf_comm->nGrowVect(),
                                        IntVect(0), period);
                // We now need to create a mask to fix the double counting.
                auto owner_mask = amrex::OwnerMask(fine_lev_cp, period);
                auto const& mma = owner_mask->const_arrays();
                auto const& sma = fine_lev_cp.const_arrays();
                auto const& dma = J_fp[lev][idim]->arrays();
                amrex::ParallelFor(fine_lev_cp, IntVect(0), ncomp,
                [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
                {
                    if (mma[bno](i,j,k) && sma[bno](i,j,k,n) != 0.0_rt) {
                        dma[bno](i,j,k,n) += sma[bno](i,j,k,n);
                    }
                });
                // Now it's safe to apply filter and sumboundary on J_cp
                if (use_filter)
                {
                    ApplyFilterJ(J_cp, lev+1, idim);
                }
                SumBoundaryJ(J_cp, lev+1, idim, period);
            }

            if (lev > 0)
            {
                // On a fine level, we need to coarsen the current onto the
                // coarse level. This needs to be done before filtering because
                // filtering depends on the level. This is also done before any
                // same-level communication because it's easier this way to
                // avoid double counting.
                J_cp[lev][idim]->setVal(0.0);
                ablastr::coarsen::average::Coarsen(*J_cp[lev][idim],
                                                   *J_fp[lev][idim],
                                                   refRatio(lev-1));
                if (current_buf[lev][idim])
                {
                    IntVect const& ng = J_cp[lev][idim]->nGrowVect();
                    AMREX_ASSERT(ng.allLE(current_buf[lev][idim]->nGrowVect()));
                    MultiFab::Add(*current_buf[lev][idim], *J_cp[lev][idim],
                                  0, 0, ncomp, ng);
                    mf_comm = std::make_unique<MultiFab>
                        (*current_buf[lev][idim], amrex::make_alias, 0, ncomp);
                }
                else
                {
                    mf_comm = std::make_unique<MultiFab>
                        (*J_cp[lev][idim], amrex::make_alias, 0, ncomp);
                }
            }

            if (use_filter)
            {
                ApplyFilterJ(J_fp, lev, idim);
            }
            SumBoundaryJ(J_fp, lev, idim, period);
        }
    }
}

void
WarpX::SyncRho () {
    SyncRho(rho_fp, rho_cp);
}

void
WarpX::SyncRho (
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_fp,
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_cp)
{
    WARPX_PROFILE("WarpX::SyncRho()");

    if (!charge_fp[0]) return;
    const int ncomp = charge_fp[0]->nComp();

    // See comments in WarpX::SyncCurrent for an explanation of the algorithm.

    std::unique_ptr<MultiFab> mf_comm; // for communication between levels
    for (int lev = finest_level; lev >= 0; --lev)
    {
        if (lev < finest_level)
        {
            auto const& period = Geom(lev).periodicity();

            // On a coarse level, the data in mf_comm comes from the
            // coarse patch of the fine level. They are unfiltered and uncommunicated.
            // We need to add it to the fine patch of the current level.
            MultiFab fine_lev_cp(charge_fp[lev]->boxArray(),
                                 charge_fp[lev]->DistributionMap(),
                                 ncomp, 0);
            fine_lev_cp.setVal(0.0);
            fine_lev_cp.ParallelAdd(*mf_comm, 0, 0, ncomp, mf_comm->nGrowVect(),
                                    IntVect(0), period);
            // We now need to create a mask to fix the double counting.
            auto owner_mask = amrex::OwnerMask(fine_lev_cp, period);
            auto const& mma = owner_mask->const_arrays();
            auto const& sma = fine_lev_cp.const_arrays();
            auto const& dma = charge_fp[lev]->arrays();
            amrex::ParallelFor(fine_lev_cp, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
            {
                if (mma[bno](i,j,k) && sma[bno](i,j,k,n) != 0.0_rt) {
                    dma[bno](i,j,k,n) += sma[bno](i,j,k,n);
                }
            });
            // Now it's safe to apply filter and sumboundary on charge_cp
            ApplyFilterandSumBoundaryRho(lev+1, lev, *charge_cp[lev+1], 0, ncomp);
        }

        if (lev > 0)
        {
            // On a fine level, we need to coarsen the data onto the coarse
            // level. This needs to be done before filtering because
            // filtering depends on the level. This is also done before any
            // same-level communication because it's easier this way to
            // avoid double counting.
            charge_cp[lev]->setVal(0.0);
            ablastr::coarsen::average::Coarsen(*charge_cp[lev],
                                               *charge_fp[lev],
                                               refRatio(lev-1));
            if (charge_buf[lev])
            {
                IntVect const& ng = charge_cp[lev]->nGrowVect();
                AMREX_ASSERT(ng.allLE(charge_buf[lev]->nGrowVect()));
                MultiFab::Add(*charge_buf[lev], *charge_cp[lev], 0, 0, ncomp, ng);
                mf_comm = std::make_unique<MultiFab>
                    (*charge_buf[lev], amrex::make_alias, 0, ncomp);
            }
            else
            {
                mf_comm = std::make_unique<MultiFab>
                    (*charge_cp[lev], amrex::make_alias, 0, ncomp);
            }
        }

        ApplyFilterandSumBoundaryRho(lev, lev, *charge_fp[lev], 0, ncomp);
    }
}

/** \brief Fills the values of the current on the coarse patch by
 *  averaging the values of the current of the fine patch (on the same level).
 */
void WarpX::RestrictCurrentFromFineToCoarsePatch (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp,
    const int lev)
{
    J_cp[lev][0]->setVal(0.0);
    J_cp[lev][1]->setVal(0.0);
    J_cp[lev][2]->setVal(0.0);

    const IntVect& refinement_ratio = refRatio(lev-1);

    std::array<const MultiFab*,3> fine { J_fp[lev][0].get(),
                                         J_fp[lev][1].get(),
                                         J_fp[lev][2].get() };
    std::array<      MultiFab*,3> crse { J_cp[lev][0].get(),
                                         J_cp[lev][1].get(),
                                         J_cp[lev][2].get() };
    ablastr::coarsen::average::Coarsen(*crse[0], *fine[0], refinement_ratio );
    ablastr::coarsen::average::Coarsen(*crse[1], *fine[1], refinement_ratio );
    ablastr::coarsen::average::Coarsen(*crse[2], *fine[2], refinement_ratio );
}

void WarpX::ApplyFilterJ (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& current,
    const int lev,
    const int idim)
{
    amrex::MultiFab& J = *current[lev][idim];

    const int ncomp = J.nComp();
    const amrex::IntVect ngrow = J.nGrowVect();
    amrex::MultiFab Jf(J.boxArray(), J.DistributionMap(), ncomp, ngrow);
    bilinear_filter.ApplyStencil(Jf, J, lev);

    const int srccomp = 0;
    const int dstcomp = 0;
    amrex::MultiFab::Copy(J, Jf, srccomp, dstcomp, ncomp, ngrow);
}

void WarpX::ApplyFilterJ (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& current,
    const int lev)
{
    for (int idim=0; idim<3; ++idim)
    {
        ApplyFilterJ(current, lev, idim);
    }
}

void WarpX::SumBoundaryJ (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& current,
    const int lev,
    const int idim,
    const amrex::Periodicity& period)
{
    amrex::MultiFab& J = *current[lev][idim];

    const amrex::IntVect ng = J.nGrowVect();
    amrex::IntVect ng_depos_J = get_ng_depos_J();

    if (WarpX::do_current_centering)
    {
#if   defined(WARPX_DIM_1D_Z)
        ng_depos_J[0] += WarpX::current_centering_noz / 2;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        ng_depos_J[0] += WarpX::current_centering_nox / 2;
        ng_depos_J[1] += WarpX::current_centering_noz / 2;
#elif defined(WARPX_DIM_3D)
        ng_depos_J[0] += WarpX::current_centering_nox / 2;
        ng_depos_J[1] += WarpX::current_centering_noy / 2;
        ng_depos_J[2] += WarpX::current_centering_noz / 2;
#endif
    }

    if (use_filter)
    {
        ng_depos_J += bilinear_filter.stencil_length_each_dir - amrex::IntVect(1);
    }

    ng_depos_J.min(ng);

    const amrex::IntVect src_ngrow = ng_depos_J;
    const int icomp = 0;
    const int ncomp = J.nComp();
    WarpXSumGuardCells(J, period, src_ngrow, icomp, ncomp);
}

void WarpX::SumBoundaryJ (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& current,
    const int lev,
    const amrex::Periodicity& period)
{
    for (int idim=0; idim<3; ++idim)
    {
        SumBoundaryJ(current, lev, idim, period);
    }
}

/* /brief Update the currents of `lev` by adding the currents from particles
*         that are in the mesh refinement patches at `lev+1`
*
* More precisely, apply filter and sum boundaries for the current of:
* - the fine patch of `lev`
* - the coarse patch of `lev+1` (same resolution)
* - the buffer regions of the coarse patch of `lev+1` (i.e. for particules
* that are within the mesh refinement patch, but do not deposit on the
* mesh refinement patch because they are too close to the boundary)
*
* Then update the fine patch of `lev` by adding the currents for the coarse
* patch (and buffer region) of `lev+1`
*/
void WarpX::AddCurrentFromFineLevelandSumBoundary (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp,
    const int lev)
{
    const amrex::Periodicity& period = Geom(lev).periodicity();

    if (use_filter)
    {
        ApplyFilterJ(J_fp, lev);
    }
    SumBoundaryJ(J_fp, lev, period);

    if (lev < finest_level)
    {
        // When there are current buffers, unlike coarse patch,
        // we don't care about the final state of them.

        for (int idim=0; idim<3; ++idim)
        {
            MultiFab mf(J_fp[lev][idim]->boxArray(),
                        J_fp[lev][idim]->DistributionMap(), J_fp[lev][idim]->nComp(), 0);
            mf.setVal(0.0);

            const IntVect ng = J_cp[lev+1][idim]->nGrowVect();

            if (use_filter && current_buf[lev+1][idim])
            {
                ApplyFilterJ(J_cp, lev+1, idim);
                ApplyFilterJ(current_buf, lev+1, idim);

                MultiFab::Add(
                    *current_buf[lev+1][idim], *J_cp[lev+1][idim],
                    0, 0, current_buf[lev+1][idim]->nComp(), ng);

                ablastr::utils::communication::ParallelAdd(
                    mf, *current_buf[lev+1][idim], 0, 0,
                    current_buf[lev+1][idim]->nComp(),
                    ng, amrex::IntVect(0),
                    do_single_precision_comms, period);
            }
            else if (use_filter) // but no buffer
            {
                ApplyFilterJ(J_cp, lev+1, idim);

                ablastr::utils::communication::ParallelAdd(
                    mf, *J_cp[lev+1][idim], 0, 0,
                    J_cp[lev+1][idim]->nComp(),
                    ng, amrex::IntVect(0),
                    do_single_precision_comms, period);
            }
            else if (current_buf[lev+1][idim]) // but no filter
            {
                MultiFab::Add(
                    *current_buf[lev+1][idim], *J_cp[lev+1][idim],
                    0, 0, current_buf[lev+1][idim]->nComp(), ng);

                ablastr::utils::communication::ParallelAdd(
                    mf, *current_buf[lev+1][idim], 0, 0,
                    current_buf[lev+1][idim]->nComp(),
                    ng, amrex::IntVect(0),
                    do_single_precision_comms, period);
            }
            else // no filter, no buffer
            {
                ablastr::utils::communication::ParallelAdd(
                    mf, *J_cp[lev+1][idim], 0, 0,
                    J_cp[lev+1][idim]->nComp(),
                    ng, amrex::IntVect(0),
                    do_single_precision_comms, period);
            }
            SumBoundaryJ(J_cp, lev+1, idim, period);
            MultiFab::Add(*J_fp[lev][idim], mf, 0, 0, J_fp[lev+1][idim]->nComp(), 0);
        }
    }
}

void WarpX::RestrictRhoFromFineToCoarsePatch (
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_fp,
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_cp,
    const int lev)
{
    if (charge_fp[lev]) {
        charge_cp[lev]->setVal(0.0);
        const IntVect& refinement_ratio = refRatio(lev-1);
        ablastr::coarsen::average::Coarsen(*charge_cp[lev], *charge_fp[lev], refinement_ratio );
    }
}

void WarpX::ApplyFilterandSumBoundaryRho (
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_fp,
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_cp,
    const int lev,
    PatchType patch_type,
    const int icomp,
    const int ncomp)
{
    const int glev = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::unique_ptr<amrex::MultiFab>& rho = (patch_type == PatchType::fine) ?
                                                  charge_fp[lev] : charge_cp[lev];
    if (rho == nullptr) return;
    ApplyFilterandSumBoundaryRho(lev, glev, *rho, icomp, ncomp);
}

void WarpX::ApplyFilterandSumBoundaryRho (int /*lev*/, int glev, amrex::MultiFab& rho, int icomp, int ncomp)
{
    const amrex::Periodicity& period = Geom(glev).periodicity();
    IntVect ng = rho.nGrowVect();
    IntVect ng_depos_rho = get_ng_depos_rho();
    if (use_filter) {
        ng += bilinear_filter.stencil_length_each_dir-1;
        ng_depos_rho += bilinear_filter.stencil_length_each_dir-1;
        ng_depos_rho.min(ng);
        MultiFab rf(rho.boxArray(), rho.DistributionMap(), ncomp, ng);
        bilinear_filter.ApplyStencil(rf, rho, glev, icomp, 0, ncomp);
        WarpXSumGuardCells(rho, rf, period, ng_depos_rho, icomp, ncomp );
    } else {
        ng_depos_rho.min(ng);
        WarpXSumGuardCells(rho, period, ng_depos_rho, icomp, ncomp);
    }
}

/* /brief Update the charge density of `lev` by adding the charge density from particles
*         that are in the mesh refinement patches at `lev+1`
*
* More precisely, apply filter and sum boundaries for the charge density of:
* - the fine patch of `lev`
* - the coarse patch of `lev+1` (same resolution)
* - the buffer regions of the coarse patch of `lev+1` (i.e. for particules
* that are within the mesh refinement patch, but do not deposit on the
* mesh refinement patch because they are too close to the boundary)
*
* Then update the fine patch of `lev` by adding the charge density for the coarse
* patch (and buffer region) of `lev+1`
*/
void WarpX::AddRhoFromFineLevelandSumBoundary (
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_fp,
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_cp,
    const int lev,
    const int icomp,
    const int ncomp)
{
    if (!charge_fp[lev]) return;

    ApplyFilterandSumBoundaryRho(charge_fp, charge_cp, lev, PatchType::fine, icomp, ncomp);

    if (lev < finest_level){

        const amrex::Periodicity& period = Geom(lev).periodicity();
        MultiFab mf(charge_fp[lev]->boxArray(),
                    charge_fp[lev]->DistributionMap(),
                    ncomp, 0);
        mf.setVal(0.0);
        IntVect ng = charge_cp[lev+1]->nGrowVect();
        IntVect ng_depos_rho = get_ng_depos_rho();
        if (use_filter && charge_buf[lev+1])
        {
            // coarse patch of fine level
            ng += bilinear_filter.stencil_length_each_dir-1;
            ng_depos_rho += bilinear_filter.stencil_length_each_dir-1;
            ng_depos_rho.min(ng);
            MultiFab rhofc(charge_cp[lev+1]->boxArray(),
                           charge_cp[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rhofc, *charge_cp[lev+1], lev+1, icomp, 0, ncomp);

            // buffer patch of fine level
            MultiFab rhofb(charge_buf[lev+1]->boxArray(),
                           charge_buf[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rhofb, *charge_buf[lev+1], lev+1, icomp, 0, ncomp);

            MultiFab::Add(rhofb, rhofc, 0, 0, ncomp, ng);

            ablastr::utils::communication::ParallelAdd(mf, rhofb, 0, 0, ncomp, ng, IntVect::TheZeroVector(),
                                                       WarpX::do_single_precision_comms, period);
            WarpXSumGuardCells( *charge_cp[lev+1], rhofc, period, ng_depos_rho, icomp, ncomp );
        }
        else if (use_filter) // but no buffer
        {
            ng += bilinear_filter.stencil_length_each_dir-1;
            ng_depos_rho += bilinear_filter.stencil_length_each_dir-1;
            ng_depos_rho.min(ng);
            MultiFab rf(charge_cp[lev+1]->boxArray(), charge_cp[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rf, *charge_cp[lev+1], lev+1, icomp, 0, ncomp);

            ablastr::utils::communication::ParallelAdd(mf, rf, 0, 0, ncomp, ng, IntVect::TheZeroVector(),
                                                       WarpX::do_single_precision_comms, period);
            WarpXSumGuardCells( *charge_cp[lev+1], rf, period, ng_depos_rho, icomp, ncomp );
        }
        else if (charge_buf[lev+1]) // but no filter
        {
            ng_depos_rho.min(ng);
            MultiFab::Add(*charge_buf[lev+1],
                          *charge_cp[lev+1], icomp, icomp, ncomp,
                           charge_cp[lev+1]->nGrowVect());

            ablastr::utils::communication::ParallelAdd(mf, *charge_buf[lev + 1], icomp, 0,
                                                       ncomp,
                                                       charge_buf[lev + 1]->nGrowVect(),
                                                       IntVect::TheZeroVector(), WarpX::do_single_precision_comms,
                                                       period);
            WarpXSumGuardCells(*(charge_cp[lev+1]), period, ng_depos_rho, icomp, ncomp);
        }
        else // no filter, no buffer
        {
            ng_depos_rho.min(ng);
            ablastr::utils::communication::ParallelAdd(mf, *charge_cp[lev + 1], icomp, 0, ncomp,
                                                       charge_cp[lev + 1]->nGrowVect(),
                                                       IntVect::TheZeroVector(), WarpX::do_single_precision_comms,
                                                       period);
            WarpXSumGuardCells(*(charge_cp[lev+1]), period, ng_depos_rho, icomp, ncomp);
        }
        MultiFab::Add(*charge_fp[lev], mf, 0, icomp, ncomp, 0);
    }
}

void WarpX::NodalSyncJ (
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
    const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp,
    const int lev,
    PatchType patch_type)
{
    if (!override_sync_intervals.contains(istep[0])) return;

    if (patch_type == PatchType::fine)
    {
        const amrex::Periodicity& period = Geom(lev).periodicity();
        ablastr::utils::communication::OverrideSync(*J_fp[lev][0], WarpX::do_single_precision_comms, period);
        ablastr::utils::communication::OverrideSync(*J_fp[lev][1], WarpX::do_single_precision_comms, period);
        ablastr::utils::communication::OverrideSync(*J_fp[lev][2], WarpX::do_single_precision_comms, period);
    }
    else if (patch_type == PatchType::coarse)
    {
        const amrex::Periodicity& cperiod = Geom(lev-1).periodicity();
        ablastr::utils::communication::OverrideSync(*J_cp[lev][0], WarpX::do_single_precision_comms, cperiod);
        ablastr::utils::communication::OverrideSync(*J_cp[lev][1], WarpX::do_single_precision_comms, cperiod);
        ablastr::utils::communication::OverrideSync(*J_cp[lev][2], WarpX::do_single_precision_comms, cperiod);
    }
}

void WarpX::NodalSyncRho (
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_fp,
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& charge_cp,
    const int lev,
    PatchType patch_type,
    const int icomp,
    const int ncomp)
{
    if (!override_sync_intervals.contains(istep[0])) return;

    if (patch_type == PatchType::fine && charge_fp[lev])
    {
        const amrex::Periodicity& period = Geom(lev).periodicity();
        MultiFab rhof(*charge_fp[lev], amrex::make_alias, icomp, ncomp);
        ablastr::utils::communication::OverrideSync(rhof, WarpX::do_single_precision_comms, period);
    }
    else if (patch_type == PatchType::coarse && charge_cp[lev])
    {
        const amrex::Periodicity& cperiod = Geom(lev-1).periodicity();
        MultiFab rhoc(*charge_cp[lev], amrex::make_alias, icomp, ncomp);
        ablastr::utils::communication::OverrideSync(rhoc, WarpX::do_single_precision_comms, cperiod);
    }
}

void WarpX::NodalSyncPML ()
{
    for (int lev = 0; lev <= finest_level; lev++) {
        NodalSyncPML(lev);
    }
}

void WarpX::NodalSyncPML (int lev)
{
    NodalSyncPML(lev, PatchType::fine);
    if (lev > 0) NodalSyncPML(lev, PatchType::coarse);
}

void WarpX::NodalSyncPML (int lev, PatchType patch_type)
{
    if (pml[lev] && pml[lev]->ok())
    {
        const std::array<amrex::MultiFab*,3>& pml_E = (patch_type == PatchType::fine) ?
                                                      pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
        const std::array<amrex::MultiFab*,3>& pml_B = (patch_type == PatchType::fine) ?
                                                      pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        amrex::MultiFab* const pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();
        amrex::MultiFab* const pml_G = (patch_type == PatchType::fine) ? pml[lev]->GetG_fp() : pml[lev]->GetG_cp();

        // Always synchronize nodal points
        const amrex::Periodicity& period = Geom(lev).periodicity();
        ablastr::utils::communication::OverrideSync(*pml_E[0], WarpX::do_single_precision_comms, period);
        ablastr::utils::communication::OverrideSync(*pml_E[1], WarpX::do_single_precision_comms, period);
        ablastr::utils::communication::OverrideSync(*pml_E[2], WarpX::do_single_precision_comms, period);
        ablastr::utils::communication::OverrideSync(*pml_B[0], WarpX::do_single_precision_comms, period);
        ablastr::utils::communication::OverrideSync(*pml_B[1], WarpX::do_single_precision_comms, period);
        ablastr::utils::communication::OverrideSync(*pml_B[2], WarpX::do_single_precision_comms, period);
        if (pml_F) {
            ablastr::utils::communication::OverrideSync(*pml_F, WarpX::do_single_precision_comms, period);
        }
        if (pml_G) {
            ablastr::utils::communication::OverrideSync(*pml_G, WarpX::do_single_precision_comms, period);
        }
    }

#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
    if (pml_rz[lev])
    {
        // This is not actually needed with RZ PSATD since the
        // arrays are always cell centered. Keep for now since
        // it may be useful if the PML is used with FDTD
        const std::array<amrex::MultiFab*,2> pml_rz_E = pml_rz[lev]->GetE_fp();
        const std::array<amrex::MultiFab*,2> pml_rz_B = pml_rz[lev]->GetB_fp();

        // Always synchronize nodal points
        const amrex::Periodicity& period = Geom(lev).periodicity();
        pml_rz_E[0]->OverrideSync(period);
        pml_rz_E[1]->OverrideSync(period);
        pml_rz_B[0]->OverrideSync(period);
        pml_rz_B[1]->OverrideSync(period);
    }
#endif
}
