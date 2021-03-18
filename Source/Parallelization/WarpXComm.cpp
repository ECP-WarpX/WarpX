/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXComm_K.H"
#include "WarpX.H"
#include "WarpXSumGuardCells.H"
#include "Utils/CoarsenMR.H"
#ifdef WARPX_USE_PSATD
#include "FieldSolver/SpectralSolver/SpectralKSpace.H"
#endif

#include <algorithm>
#include <cstdlib>
#include <memory>

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
    if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
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
    for (MFIter mfi(*Bfield_aux[0][0]); mfi.isValid(); ++mfi)
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

        // Loop over full box including ghost cells
        // (input arrays will be padded with zeros beyond ghost cells
        // for out-of-bound accesses due to large-stencil operations)
        Box bx = mfi.fabbox();

        if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {

#ifdef WARPX_USE_PSATD

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
                warpx_interp<true>(j, k, l, bx_aux, bx_fp, dst_stag, Bx_stag, fg_nox, fg_noy, fg_noz,
                             stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

                warpx_interp<true>(j, k, l, by_aux, by_fp, dst_stag, By_stag, fg_nox, fg_noy, fg_noz,
                             stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

                warpx_interp<true>(j, k, l, bz_aux, bz_fp, dst_stag, Bz_stag, fg_nox, fg_noy, fg_noz,
                             stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

                warpx_interp<true>(j, k, l, ex_aux, ex_fp, dst_stag, Ex_stag, fg_nox, fg_noy, fg_noz,
                             stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

                warpx_interp<true>(j, k, l, ey_aux, ey_fp, dst_stag, Ey_stag, fg_nox, fg_noy, fg_noz,
                             stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);

                warpx_interp<true>(j, k, l, ez_aux, ez_fp, dst_stag, Ez_stag, fg_nox, fg_noy, fg_noz,
                             stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
            });
#endif
        } else { // FDTD
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
            {
                warpx_interp(j, k, l, bx_aux, bx_fp, dst_stag, Bx_stag);
                warpx_interp(j, k, l, by_aux, by_fp, dst_stag, By_stag);
                warpx_interp(j, k, l, bz_aux, bz_fp, dst_stag, Bz_stag);
                warpx_interp(j, k, l, ex_aux, ex_fp, dst_stag, Ex_stag);
                warpx_interp(j, k, l, ey_aux, ey_fp, dst_stag, Ey_stag);
                warpx_interp(j, k, l, ez_aux, ez_fp, dst_stag, Ez_stag);
            });
        }
    }

    // NOTE: high-order interpolation is not implemented for mesh refinement
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        BoxArray const& nba = Bfield_aux[lev][0]->boxArray();
        BoxArray const& cnba = amrex::coarsen(nba,2);
        DistributionMapping const& dm = Bfield_aux[lev][0]->DistributionMap();
        auto const& cperiod = Geom(lev-1).periodicity();

        // Bfield
        {
            Array<std::unique_ptr<MultiFab>,3> Btmp;
            if (Bfield_cax[lev][0]) {
                for (int i = 0; i < 3; ++i) {
                    Btmp[i] = std::make_unique<MultiFab>(
                        *Bfield_cax[lev][i], amrex::make_alias, 0, 1);
                }
            } else {
                IntVect ngtmp = Bfield_aux[lev-1][0]->nGrowVect();
                for (int i = 0; i < 3; ++i) {
                    Btmp[i] = std::make_unique<MultiFab>(cnba, dm, 1, ngtmp);
                }
            }
            // ParallelCopy from coarse level
            for (int i = 0; i < 3; ++i) {
                IntVect ng = Btmp[i]->nGrowVect();
                Btmp[i]->ParallelCopy(*Bfield_aux[lev-1][i], 0, 0, 1, ng, ng, cperiod);
            }

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
                Array4<Real const> const& bx_cp = Bfield_cp[lev][0]->const_array(mfi);
                Array4<Real const> const& by_cp = Bfield_cp[lev][1]->const_array(mfi);
                Array4<Real const> const& bz_cp = Bfield_cp[lev][2]->const_array(mfi);
                Array4<Real const> const& bx_c = Btmp[0]->const_array(mfi);
                Array4<Real const> const& by_c = Btmp[1]->const_array(mfi);
                Array4<Real const> const& bz_c = Btmp[2]->const_array(mfi);

                const Box& bx = mfi.fabbox();
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
                IntVect ngtmp = Efield_aux[lev-1][0]->nGrowVect();
                for (int i = 0; i < 3; ++i) {
                    Etmp[i] = std::make_unique<MultiFab>(
                        cnba, dm, 1, ngtmp);
                }
            }
            // ParallelCopy from coarse level
            for (int i = 0; i < 3; ++i) {
                IntVect ng = Etmp[i]->nGrowVect();
                Etmp[i]->ParallelCopy(*Efield_aux[lev-1][i], 0, 0, 1, ng, ng, cperiod);
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
        const auto& crse_period = Geom(lev-1).periodicity();
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
            dBx.ParallelCopy(*Bfield_aux[lev-1][0], 0, 0, Bfield_aux[lev-1][0]->nComp(), ng, ng, crse_period);
            dBy.ParallelCopy(*Bfield_aux[lev-1][1], 0, 0, Bfield_aux[lev-1][1]->nComp(), ng, ng, crse_period);
            dBz.ParallelCopy(*Bfield_aux[lev-1][2], 0, 0, Bfield_aux[lev-1][2]->nComp(), ng, ng, crse_period);
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
            dEx.ParallelCopy(*Efield_aux[lev-1][0], 0, 0, Efield_aux[lev-1][0]->nComp(), ng, ng, crse_period);
            dEy.ParallelCopy(*Efield_aux[lev-1][1], 0, 0, Efield_aux[lev-1][1]->nComp(), ng, ng, crse_period);
            dEz.ParallelCopy(*Efield_aux[lev-1][2], 0, 0, Efield_aux[lev-1][2]->nComp(), ng, ng, crse_period);
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

    for (MFIter mfi(dst); mfi.isValid(); ++mfi)
    {
        // Loop over full box including ghost cells
        // (input arrays will be padded with zeros beyond ghost cells
        // for out-of-bound accesses due to large-stencil operations)
        Box bx = mfi.fabbox();

        amrex::Array4<amrex::Real const> const& src_arr = src.const_array(mfi);
        amrex::Array4<amrex::Real>       const& dst_arr = dst.array(mfi);

        if (maxwell_solver_id == MaxwellSolverAlgo::PSATD)
        {
#ifdef WARPX_USE_PSATD

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
                warpx_interp<true>(j, k, l, dst_arr, src_arr, dst_stag, src_stag, cc_nox, cc_noy, cc_noz,
                                   stencil_coeffs_x, stencil_coeffs_y, stencil_coeffs_z);
            });
#endif
        }

        else // FDTD
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int j, int k, int l) noexcept
            {
                warpx_interp<false>(j, k, l, dst_arr, src_arr, dst_stag, src_stag);
            });
        }
    }
}

void
WarpX::FillBoundaryB (IntVect ng)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryB(lev, ng);
    }
}

void
WarpX::FillBoundaryE (IntVect ng)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryE(lev, ng);
    }
}

void
WarpX::FillBoundaryF (IntVect ng)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryF(lev, ng);
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
WarpX::FillBoundaryE(int lev, IntVect ng)
{
    FillBoundaryE(lev, PatchType::fine, ng);
    if (lev > 0) FillBoundaryE(lev, PatchType::coarse, ng);
}

void
WarpX::FillBoundaryE (int lev, PatchType patch_type, IntVect ng)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev]->ok())
        {
            pml[lev]->ExchangeE(patch_type,
                                { Efield_fp[lev][0].get(),
                                  Efield_fp[lev][1].get(),
                                  Efield_fp[lev][2].get() },
                                do_pml_in_domain);
            pml[lev]->FillBoundaryE(patch_type);
        }

        const auto& period = Geom(lev).periodicity();
        if ( safe_guard_cells ){
            Vector<MultiFab*> mf{Efield_fp[lev][0].get(),Efield_fp[lev][1].get(),Efield_fp[lev][2].get()};
            amrex::FillBoundary(mf, period);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Efield_fp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryE, requested more guard cells than allocated");
            Efield_fp[lev][0]->FillBoundary(ng, period);
            Efield_fp[lev][1]->FillBoundary(ng, period);
            Efield_fp[lev][2]->FillBoundary(ng, period);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
        {
            pml[lev]->ExchangeE(patch_type,
                                { Efield_cp[lev][0].get(),
                                  Efield_cp[lev][1].get(),
                                  Efield_cp[lev][2].get() },
                                do_pml_in_domain);
            pml[lev]->FillBoundaryE(patch_type);
        }
        const auto& cperiod = Geom(lev-1).periodicity();
        if ( safe_guard_cells ) {
            Vector<MultiFab*> mf{Efield_cp[lev][0].get(),Efield_cp[lev][1].get(),Efield_cp[lev][2].get()};
            amrex::FillBoundary(mf, cperiod);

        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Efield_cp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryE, requested more guard cells than allocated");
            Efield_cp[lev][0]->FillBoundary(ng, cperiod);
            Efield_cp[lev][1]->FillBoundary(ng, cperiod);
            Efield_cp[lev][2]->FillBoundary(ng, cperiod);
        }
    }
}

void
WarpX::FillBoundaryB (int lev, IntVect ng)
{
    FillBoundaryB(lev, PatchType::fine, ng);
    if (lev > 0) FillBoundaryB(lev, PatchType::coarse, ng);
}

void
WarpX::FillBoundaryB (int lev, PatchType patch_type, IntVect ng)
{
    if (patch_type == PatchType::fine)
    {
        if (do_pml && pml[lev]->ok())
        {
            pml[lev]->ExchangeB(patch_type,
                            { Bfield_fp[lev][0].get(),
                              Bfield_fp[lev][1].get(),
                              Bfield_fp[lev][2].get() },
                              do_pml_in_domain);
        pml[lev]->FillBoundaryB(patch_type);
        }
        const auto& period = Geom(lev).periodicity();
        if ( safe_guard_cells ) {
            Vector<MultiFab*> mf{Bfield_fp[lev][0].get(),Bfield_fp[lev][1].get(),Bfield_fp[lev][2].get()};
            amrex::FillBoundary(mf, period);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Bfield_fp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryB, requested more guard cells than allocated");
            Bfield_fp[lev][0]->FillBoundary(ng, period);
            Bfield_fp[lev][1]->FillBoundary(ng, period);
            Bfield_fp[lev][2]->FillBoundary(ng, period);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
        {
        pml[lev]->ExchangeB(patch_type,
                      { Bfield_cp[lev][0].get(),
                        Bfield_cp[lev][1].get(),
                        Bfield_cp[lev][2].get() },
                        do_pml_in_domain);
        pml[lev]->FillBoundaryB(patch_type);
        }
        const auto& cperiod = Geom(lev-1).periodicity();
        if ( safe_guard_cells ){
            Vector<MultiFab*> mf{Bfield_cp[lev][0].get(),Bfield_cp[lev][1].get(),Bfield_cp[lev][2].get()};
            amrex::FillBoundary(mf, cperiod);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Bfield_cp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryB, requested more guard cells than allocated");
            Bfield_cp[lev][0]->FillBoundary(ng, cperiod);
            Bfield_cp[lev][1]->FillBoundary(ng, cperiod);
            Bfield_cp[lev][2]->FillBoundary(ng, cperiod);
        }
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
            amrex::Abort("Averaged Galilean PSATD with PML is not yet implemented");
         }

        const auto& period = Geom(lev).periodicity();
        if ( safe_guard_cells ){
            Vector<MultiFab*> mf{Efield_avg_fp[lev][0].get(),Efield_avg_fp[lev][1].get(),Efield_avg_fp[lev][2].get()};
            amrex::FillBoundary(mf, period);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Efield_avg_fp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryE_avg, requested more guard cells than allocated");
            Efield_avg_fp[lev][0]->FillBoundary(ng, period);
            Efield_avg_fp[lev][1]->FillBoundary(ng, period);
            Efield_avg_fp[lev][2]->FillBoundary(ng, period);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
         {
            amrex::Abort("Averaged Galilean PSATD with PML is not yet implemented");
         }

        const auto& cperiod = Geom(lev-1).periodicity();
        if ( safe_guard_cells ) {
            Vector<MultiFab*> mf{Efield_avg_cp[lev][0].get(),Efield_avg_cp[lev][1].get(),Efield_avg_cp[lev][2].get()};
            amrex::FillBoundary(mf, cperiod);

        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Efield_avg_cp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryE, requested more guard cells than allocated");
            Efield_avg_cp[lev][0]->FillBoundary(ng, cperiod);
            Efield_avg_cp[lev][1]->FillBoundary(ng, cperiod);
            Efield_avg_cp[lev][2]->FillBoundary(ng, cperiod);
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
            amrex::Abort("Averaged Galilean PSATD with PML is not yet implemented");
          }
        const auto& period = Geom(lev).periodicity();
        if ( safe_guard_cells ) {
            Vector<MultiFab*> mf{Bfield_avg_fp[lev][0].get(),Bfield_avg_fp[lev][1].get(),Bfield_avg_fp[lev][2].get()};
            amrex::FillBoundary(mf, period);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Bfield_fp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryB, requested more guard cells than allocated");
            Bfield_avg_fp[lev][0]->FillBoundary(ng, period);
            Bfield_avg_fp[lev][1]->FillBoundary(ng, period);
            Bfield_avg_fp[lev][2]->FillBoundary(ng, period);
        }
    }
    else if (patch_type == PatchType::coarse)
    {
        if (do_pml && pml[lev]->ok())
          {
            amrex::Abort("Averaged Galilean PSATD with PML is not yet implemented");
          }

        const auto& cperiod = Geom(lev-1).periodicity();
        if ( safe_guard_cells ){
            Vector<MultiFab*> mf{Bfield_avg_cp[lev][0].get(),Bfield_avg_cp[lev][1].get(),Bfield_avg_cp[lev][2].get()};
            amrex::FillBoundary(mf, cperiod);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= Bfield_avg_cp[lev][0]->nGrowVect(),
                "Error: in FillBoundaryB_avg, requested more guard cells than allocated");
            Bfield_avg_cp[lev][0]->FillBoundary(ng, cperiod);
            Bfield_avg_cp[lev][1]->FillBoundary(ng, cperiod);
            Bfield_avg_cp[lev][2]->FillBoundary(ng, cperiod);
        }
    }
}


void
WarpX::FillBoundaryF (int lev, IntVect ng)
{
    FillBoundaryF(lev, PatchType::fine, ng);
    if (lev > 0) FillBoundaryF(lev, PatchType::coarse, ng);
}

void
WarpX::FillBoundaryF (int lev, PatchType patch_type, IntVect ng)
{
    if (patch_type == PatchType::fine && F_fp[lev])
    {
        if (do_pml && pml[lev]->ok())
        {
            pml[lev]->ExchangeF(patch_type, F_fp[lev].get(),
                                do_pml_in_domain);
            pml[lev]->FillBoundaryF(patch_type);
        }

        const auto& period = Geom(lev).periodicity();
        if ( safe_guard_cells ) {
            F_fp[lev]->FillBoundary(period);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= F_fp[lev]->nGrowVect(),
                "Error: in FillBoundaryF, requested more guard cells than allocated");
            F_fp[lev]->FillBoundary(ng, period);
        }
    }
    else if (patch_type == PatchType::coarse && F_cp[lev])
    {
        if (do_pml && pml[lev]->ok())
        {
        pml[lev]->ExchangeF(patch_type, F_cp[lev].get(),
                            do_pml_in_domain);
        pml[lev]->FillBoundaryF(patch_type);
        }

        const auto& cperiod = Geom(lev-1).periodicity();
        if ( safe_guard_cells ) {
            F_cp[lev]->FillBoundary(cperiod);
        } else {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ng <= F_cp[lev]->nGrowVect(),
                "Error: in FillBoundaryF, requested more guard cells than allocated");
            F_cp[lev]->FillBoundary(ng, cperiod);
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
    const auto& period = Geom(lev).periodicity();
    Efield_aux[lev][0]->FillBoundary(ng, period);
    Efield_aux[lev][1]->FillBoundary(ng, period);
    Efield_aux[lev][2]->FillBoundary(ng, period);
    Bfield_aux[lev][0]->FillBoundary(ng, period);
    Bfield_aux[lev][1]->FillBoundary(ng, period);
    Bfield_aux[lev][2]->FillBoundary(ng, period);
}

void
WarpX::SyncCurrent ()
{
    WARPX_PROFILE("WarpX::SyncCurrent()");

    // If warpx.do_current_centering = 1, center currents from nodal grid to staggered grid
    if (WarpX::do_current_centering)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            WarpX::UpdateCurrentNodalToStag(*current_fp[lev][0], *current_fp_nodal[lev][0]);
            WarpX::UpdateCurrentNodalToStag(*current_fp[lev][1], *current_fp_nodal[lev][1]);
            WarpX::UpdateCurrentNodalToStag(*current_fp[lev][2], *current_fp_nodal[lev][2]);
        }
    }

    // Restrict fine patch current onto the coarse patch, before
    // summing the guard cells of the fine patch
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        current_cp[lev][0]->setVal(0.0);
        current_cp[lev][1]->setVal(0.0);
        current_cp[lev][2]->setVal(0.0);

        const IntVect& refinement_ratio = refRatio(lev-1);

        std::array<const MultiFab*,3> fine { current_fp[lev][0].get(),
                                             current_fp[lev][1].get(),
                                             current_fp[lev][2].get() };
        std::array<      MultiFab*,3> crse { current_cp[lev][0].get(),
                                             current_cp[lev][1].get(),
                                             current_cp[lev][2].get() };
        CoarsenMR::Coarsen( *crse[0], *fine[0], refinement_ratio );
        CoarsenMR::Coarsen( *crse[1], *fine[1], refinement_ratio );
        CoarsenMR::Coarsen( *crse[2], *fine[2], refinement_ratio );
    }

    // For each level
    // - apply filter to the coarse patch/buffer of `lev+1` and fine patch of `lev` (same resolution)
    // - add the coarse patch/buffer of `lev+1` into the fine patch of `lev`
    // - sum guard cells of the coarse patch of `lev+1` and fine patch of `lev`
    for (int lev=0; lev <= finest_level; ++lev) {
        AddCurrentFromFineLevelandSumBoundary(lev);
    }
}

void
WarpX::SyncRho ()
{
    WARPX_PROFILE("WarpX::SyncRho()");

    if (!rho_fp[0]) return;
    const int ncomp = rho_fp[0]->nComp();

    // Restrict fine patch onto the coarse patch,
    // before summing the guard cells of the fine patch
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        rho_cp[lev]->setVal(0.0);
        const IntVect& refinement_ratio = refRatio(lev-1);
        CoarsenMR::Coarsen( *rho_cp[lev], *rho_fp[lev], refinement_ratio );
    }

    // For each level
    // - apply filter to the coarse patch/buffer of `lev+1` and fine patch of `lev` (same resolution)
    // - add the coarse patch/buffer of `lev+1` into the fine patch of `lev`
    // - sum guard cells of the coarse patch of `lev+1` and fine patch of `lev`
    for (int lev=0; lev <= finest_level; ++lev) {
        AddRhoFromFineLevelandSumBoundary(lev, 0, ncomp);
    }
}

/** \brief Fills the values of the current on the coarse patch by
 *  averaging the values of the current of the fine patch (on the same level).
 */
void
WarpX::RestrictCurrentFromFineToCoarsePatch (int lev)
{
    current_cp[lev][0]->setVal(0.0);
    current_cp[lev][1]->setVal(0.0);
    current_cp[lev][2]->setVal(0.0);

    const IntVect& refinement_ratio = refRatio(lev-1);

    std::array<const MultiFab*,3> fine { current_fp[lev][0].get(),
                                         current_fp[lev][1].get(),
                                         current_fp[lev][2].get() };
    std::array<      MultiFab*,3> crse { current_cp[lev][0].get(),
                                         current_cp[lev][1].get(),
                                         current_cp[lev][2].get() };
    CoarsenMR::Coarsen( *crse[0], *fine[0], refinement_ratio );
    CoarsenMR::Coarsen( *crse[1], *fine[1], refinement_ratio );
    CoarsenMR::Coarsen( *crse[2], *fine[2], refinement_ratio );
}

void
WarpX::ApplyFilterandSumBoundaryJ (int lev, PatchType patch_type)
{
    const int glev = (patch_type == PatchType::fine) ? lev : lev-1;
    const auto& period = Geom(glev).periodicity();
    auto& j = (patch_type == PatchType::fine) ? current_fp[lev] : current_cp[lev];
    for (int idim = 0; idim < 3; ++idim) {
        if (use_filter) {
            IntVect ng = j[idim]->nGrowVect();
            ng += bilinear_filter.stencil_length_each_dir-1;
            MultiFab jf(j[idim]->boxArray(), j[idim]->DistributionMap(), j[idim]->nComp(), ng);
            bilinear_filter.ApplyStencil(jf, *j[idim]);
            WarpXSumGuardCells(*(j[idim]), jf, period, 0, (j[idim])->nComp());
        } else {
            WarpXSumGuardCells(*(j[idim]), period, 0, (j[idim])->nComp());
        }
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
void
WarpX::AddCurrentFromFineLevelandSumBoundary (int lev)
{
    ApplyFilterandSumBoundaryJ(lev, PatchType::fine);

    if (lev < finest_level) {
        // When there are current buffers, unlike coarse patch,
        // we don't care about the final state of them.

        const auto& period = Geom(lev).periodicity();
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab mf(current_fp[lev][idim]->boxArray(),
                        current_fp[lev][idim]->DistributionMap(), current_fp[lev][idim]->nComp(), 0);
            mf.setVal(0.0);
            if (use_filter && current_buf[lev+1][idim])
            {
                // coarse patch of fine level
                IntVect ng = current_cp[lev+1][idim]->nGrowVect();
                ng += bilinear_filter.stencil_length_each_dir-1;
                MultiFab jfc(current_cp[lev+1][idim]->boxArray(),
                             current_cp[lev+1][idim]->DistributionMap(), current_cp[lev+1][idim]->nComp(), ng);
                bilinear_filter.ApplyStencil(jfc, *current_cp[lev+1][idim]);

                // buffer patch of fine level
                MultiFab jfb(current_buf[lev+1][idim]->boxArray(),
                             current_buf[lev+1][idim]->DistributionMap(), current_buf[lev+1][idim]->nComp(), ng);
                bilinear_filter.ApplyStencil(jfb, *current_buf[lev+1][idim]);

                MultiFab::Add(jfb, jfc, 0, 0, current_buf[lev+1][idim]->nComp(), ng);
                mf.ParallelAdd(jfb, 0, 0, current_buf[lev+1][idim]->nComp(), ng, IntVect::TheZeroVector(), period);

                WarpXSumGuardCells(*current_cp[lev+1][idim], jfc, period, 0, current_cp[lev+1][idim]->nComp());
            }
            else if (use_filter) // but no buffer
            {
                // coarse patch of fine level
                IntVect ng = current_cp[lev+1][idim]->nGrowVect();
                ng += bilinear_filter.stencil_length_each_dir-1;
                MultiFab jf(current_cp[lev+1][idim]->boxArray(),
                            current_cp[lev+1][idim]->DistributionMap(), current_cp[lev+1][idim]->nComp(), ng);
                bilinear_filter.ApplyStencil(jf, *current_cp[lev+1][idim]);
                mf.ParallelAdd(jf, 0, 0, current_cp[lev+1][idim]->nComp(), ng, IntVect::TheZeroVector(), period);
                WarpXSumGuardCells(*current_cp[lev+1][idim], jf, period, 0, current_cp[lev+1][idim]->nComp());
            }
            else if (current_buf[lev+1][idim]) // but no filter
            {
                MultiFab::Add(*current_buf[lev+1][idim],
                               *current_cp [lev+1][idim], 0, 0, current_buf[lev+1][idim]->nComp(),
                               current_cp[lev+1][idim]->nGrowVect());
                mf.ParallelAdd(*current_buf[lev+1][idim], 0, 0, current_buf[lev+1][idim]->nComp(),
                               current_buf[lev+1][idim]->nGrowVect(), IntVect::TheZeroVector(),
                               period);
                WarpXSumGuardCells(*(current_cp[lev+1][idim]), period, 0, current_cp[lev+1][idim]->nComp());
            }
            else // no filter, no buffer
            {
                mf.ParallelAdd(*current_cp[lev+1][idim], 0, 0, current_cp[lev+1][idim]->nComp(),
                               current_cp[lev+1][idim]->nGrowVect(), IntVect::TheZeroVector(),
                               period);
                WarpXSumGuardCells(*(current_cp[lev+1][idim]), period, 0, current_cp[lev+1][idim]->nComp());
            }
            MultiFab::Add(*current_fp[lev][idim], mf, 0, 0, current_fp[lev+1][idim]->nComp(), 0);
        }
        NodalSyncJ(lev+1, PatchType::coarse);
    }
    NodalSyncJ(lev, PatchType::fine);
}

void
WarpX::RestrictRhoFromFineToCoarsePatch (int lev)
{
    if (rho_fp[lev]) {
        rho_cp[lev]->setVal(0.0);
        const IntVect& refinement_ratio = refRatio(lev-1);
        CoarsenMR::Coarsen( *rho_cp[lev], *rho_fp[lev], refinement_ratio );
    }
}

void
WarpX::ApplyFilterandSumBoundaryRho (int lev, PatchType patch_type, int icomp, int ncomp)
{
    const int glev = (patch_type == PatchType::fine) ? lev : lev-1;
    auto& r = (patch_type == PatchType::fine) ? rho_fp[lev] : rho_cp[lev];
    if (r == nullptr) return;
    ApplyFilterandSumBoundaryRho(lev, glev, *r, icomp, ncomp);
}

void
WarpX::ApplyFilterandSumBoundaryRho (int /*lev*/, int glev, amrex::MultiFab& rho, int icomp, int ncomp)
{
    const auto& period = Geom(glev).periodicity();
    if (use_filter) {
        IntVect ng = rho.nGrowVect();
        ng += bilinear_filter.stencil_length_each_dir-1;
        MultiFab rf(rho.boxArray(), rho.DistributionMap(), ncomp, ng);
        bilinear_filter.ApplyStencil(rf, rho, icomp, 0, ncomp);
        WarpXSumGuardCells(rho, rf, period, icomp, ncomp );
    } else {
        WarpXSumGuardCells(rho, period, icomp, ncomp);
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
void
WarpX::AddRhoFromFineLevelandSumBoundary(int lev, int icomp, int ncomp)
{
    if (!rho_fp[lev]) return;

    ApplyFilterandSumBoundaryRho(lev, PatchType::fine, icomp, ncomp);

    if (lev < finest_level){

        const auto& period = Geom(lev).periodicity();
        MultiFab mf(rho_fp[lev]->boxArray(),
                    rho_fp[lev]->DistributionMap(),
                    ncomp, 0);
        mf.setVal(0.0);
        if (use_filter && charge_buf[lev+1])
        {
            // coarse patch of fine level
            IntVect ng = rho_cp[lev+1]->nGrowVect();
            ng += bilinear_filter.stencil_length_each_dir-1;
            MultiFab rhofc(rho_cp[lev+1]->boxArray(),
                         rho_cp[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rhofc, *rho_cp[lev+1], icomp, 0, ncomp);

            // buffer patch of fine level
            MultiFab rhofb(charge_buf[lev+1]->boxArray(),
                           charge_buf[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rhofb, *charge_buf[lev+1], icomp, 0, ncomp);

            MultiFab::Add(rhofb, rhofc, 0, 0, ncomp, ng);
            mf.ParallelAdd(rhofb, 0, 0, ncomp, ng, IntVect::TheZeroVector(), period);
            WarpXSumGuardCells( *rho_cp[lev+1], rhofc, period, icomp, ncomp );
        }
        else if (use_filter) // but no buffer
        {
            IntVect ng = rho_cp[lev+1]->nGrowVect();
            ng += bilinear_filter.stencil_length_each_dir-1;
            MultiFab rf(rho_cp[lev+1]->boxArray(), rho_cp[lev+1]->DistributionMap(), ncomp, ng);
            bilinear_filter.ApplyStencil(rf, *rho_cp[lev+1], icomp, 0, ncomp);
            mf.ParallelAdd(rf, 0, 0, ncomp, ng, IntVect::TheZeroVector(), period);
            WarpXSumGuardCells( *rho_cp[lev+1], rf, period, icomp, ncomp );
        }
        else if (charge_buf[lev+1]) // but no filter
        {
            MultiFab::Add(*charge_buf[lev+1],
                           *rho_cp[lev+1], icomp, icomp, ncomp,
                           rho_cp[lev+1]->nGrowVect());
            mf.ParallelAdd(*charge_buf[lev+1], icomp, 0,
                           ncomp,
                           charge_buf[lev+1]->nGrowVect(), IntVect::TheZeroVector(),
                           period);
            WarpXSumGuardCells(*(rho_cp[lev+1]), period, icomp, ncomp);
        }
        else // no filter, no buffer
        {
            mf.ParallelAdd(*rho_cp[lev+1], icomp, 0, ncomp,
                           rho_cp[lev+1]->nGrowVect(), IntVect::TheZeroVector(),
                           period);
            WarpXSumGuardCells(*(rho_cp[lev+1]), period, icomp, ncomp);
        }
        MultiFab::Add(*rho_fp[lev], mf, 0, icomp, ncomp, 0);
        NodalSyncRho(lev+1, PatchType::coarse, icomp, ncomp);
    }

    NodalSyncRho(lev, PatchType::fine, icomp, ncomp);
}

void
WarpX::NodalSyncJ (int lev, PatchType patch_type)
{
    if (!override_sync_intervals.contains(istep[0])) return;

    if (patch_type == PatchType::fine)
    {
        const auto& period = Geom(lev).periodicity();
        current_fp[lev][0]->OverrideSync(period);
        current_fp[lev][1]->OverrideSync(period);
        current_fp[lev][2]->OverrideSync(period);
    }
    else if (patch_type == PatchType::coarse)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        current_cp[lev][0]->OverrideSync(cperiod);
        current_cp[lev][1]->OverrideSync(cperiod);
        current_cp[lev][2]->OverrideSync(cperiod);
    }
}

void
WarpX::NodalSyncRho (int lev, PatchType patch_type, int icomp, int ncomp)
{
    if (!override_sync_intervals.contains(istep[0])) return;

    if (patch_type == PatchType::fine && rho_fp[lev])
    {
        const auto& period = Geom(lev).periodicity();
        MultiFab rhof(*rho_fp[lev], amrex::make_alias, icomp, ncomp);
        rhof.OverrideSync(period);
    }
    else if (patch_type == PatchType::coarse && rho_cp[lev])
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        MultiFab rhoc(*rho_cp[lev], amrex::make_alias, icomp, ncomp);
        rhoc.OverrideSync(cperiod);
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
    if (pml[lev]->ok())
    {
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();

        // Always synchronize nodal points
        const auto& period = Geom(lev).periodicity();
        pml_E[0]->OverrideSync(period);
        pml_E[1]->OverrideSync(period);
        pml_E[2]->OverrideSync(period);
        pml_B[0]->OverrideSync(period);
        pml_B[1]->OverrideSync(period);
        pml_B[2]->OverrideSync(period);
        if (pml_F) {
            pml_F->OverrideSync(period);
        }
    }
}

void WarpX::NodalSyncE ()
{
    if (!override_sync_intervals.contains(istep[0]) && !do_pml) return;

    for (int lev = 0; lev <= finest_level; lev++) {
        NodalSyncE(lev);
    }
}

void WarpX::NodalSyncE (int lev)
{
    NodalSyncE(lev, PatchType::fine);
    if (lev > 0) NodalSyncE(lev, PatchType::coarse);
}

void WarpX::NodalSyncE (int lev, PatchType patch_type)
{
    if (patch_type == PatchType::fine)
    {
        const auto& period = Geom(lev).periodicity();
        Efield_fp[lev][0]->OverrideSync(period);
        Efield_fp[lev][1]->OverrideSync(period);
        Efield_fp[lev][2]->OverrideSync(period);
    }
    else if (patch_type == PatchType::coarse)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        Efield_cp[lev][0]->OverrideSync(cperiod);
        Efield_cp[lev][1]->OverrideSync(cperiod);
        Efield_cp[lev][2]->OverrideSync(cperiod);
    }
}

void WarpX::NodalSyncB ()
{
    if (!override_sync_intervals.contains(istep[0]) && !do_pml) return;

    for (int lev = 0; lev <= finest_level; lev++) {
        NodalSyncB(lev);
    }
}

void WarpX::NodalSyncB (int lev)
{
    NodalSyncB(lev, PatchType::fine);
    if (lev > 0) NodalSyncB(lev, PatchType::coarse);
}

void WarpX::NodalSyncB (int lev, PatchType patch_type)
{
    if (patch_type == PatchType::fine)
    {
        const auto& period = Geom(lev).periodicity();
        Bfield_fp[lev][0]->OverrideSync(period);
        Bfield_fp[lev][1]->OverrideSync(period);
        Bfield_fp[lev][2]->OverrideSync(period);
    }
    else if (patch_type == PatchType::coarse)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        Bfield_cp[lev][0]->OverrideSync(cperiod);
        Bfield_cp[lev][1]->OverrideSync(cperiod);
        Bfield_cp[lev][2]->OverrideSync(cperiod);
    }
}
