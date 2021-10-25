/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * Maxence Thevenet, Remi Lehe, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PML.H"

#include "BoundaryConditions/PML.H"
#include "BoundaryConditions/PMLComponent.H"
#ifdef WARPX_USE_PSATD
#   include "FieldSolver/SpectralSolver/SpectralFieldData.H"
#endif
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"
#include "Parallelization/WarpXCommUtil.H"

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FBI.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealVect.H>
#include <AMReX_SPACE.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>

using namespace amrex;

namespace
{
    static void FillLo (int idim, Sigma& sigma, Sigma& sigma_cumsum,
                        Sigma& sigma_star, Sigma& sigma_star_cumsum,
                        const Box& overlap, const Box& grid, Real fac)
    {
        int glo = grid.smallEnd(idim);
        int olo = overlap.smallEnd(idim);
        int ohi = overlap.bigEnd(idim);
        int slo = sigma.m_lo;
        int sslo = sigma_star.m_lo;

        const int N = ohi+1-olo+1;
        Real* p_sigma = sigma.data();
        Real* p_sigma_cumsum = sigma_cumsum.data();
        Real* p_sigma_star = sigma_star.data();
        Real* p_sigma_star_cumsum = sigma_star_cumsum.data();
        amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            i += olo;
            Real offset = static_cast<Real>(glo-i);
            p_sigma[i-slo] = fac*(offset*offset);
            // sigma_cumsum is the analytical integral of sigma function at same points than sigma
            p_sigma_cumsum[i-slo] = (fac*(offset*offset*offset)/3._rt)/PhysConst::c;
            if (i <= ohi+1) {
                offset = static_cast<Real>(glo-i) - 0.5_rt;
                p_sigma_star[i-sslo] = fac*(offset*offset);
                // sigma_star_cumsum is the analytical integral of sigma function at same points than sigma_star
                p_sigma_star_cumsum[i-sslo] = (fac*(offset*offset*offset)/3._rt)/PhysConst::c;
            }
        });
    }

    static void FillHi (int idim, Sigma& sigma, Sigma& sigma_cumsum,
                        Sigma& sigma_star, Sigma& sigma_star_cumsum,
                        const Box& overlap, const Box& grid, Real fac)
    {
        int ghi = grid.bigEnd(idim);
        int olo = overlap.smallEnd(idim);
        int ohi = overlap.bigEnd(idim);
        int slo = sigma.m_lo;
        int sslo = sigma_star.m_lo;

        const int N = ohi+1-olo+1;
        Real* p_sigma = sigma.data();
        Real* p_sigma_cumsum = sigma_cumsum.data();
        Real* p_sigma_star = sigma_star.data();
        Real* p_sigma_star_cumsum = sigma_star_cumsum.data();
        amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            i += olo;
            Real offset = static_cast<Real>(i-ghi-1);
            p_sigma[i-slo] = fac*(offset*offset);
            p_sigma_cumsum[i-slo] = (fac*(offset*offset*offset)/3._rt)/PhysConst::c;
            if (i <= ohi+1) {
                offset = static_cast<Real>(i-ghi) - 0.5_rt;
                p_sigma_star[i-sslo] = fac*(offset*offset);
                p_sigma_star_cumsum[i-sslo] = (fac*(offset*offset*offset)/3._rt)/PhysConst::c;
            }
        });
    }

    static void FillZero (int idim, Sigma& sigma, Sigma& sigma_cumsum,
                          Sigma& sigma_star, Sigma& sigma_star_cumsum,
                          const Box& overlap)
    {
        int olo = overlap.smallEnd(idim);
        int ohi = overlap.bigEnd(idim);
        int slo = sigma.m_lo;
        int sslo = sigma_star.m_lo;

        const int N = ohi+1-olo+1;
        Real* p_sigma = sigma.data();
        Real* p_sigma_cumsum = sigma_cumsum.data();
        Real* p_sigma_star = sigma_star.data();
        Real* p_sigma_star_cumsum = sigma_star_cumsum.data();
        amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            i += olo;
            p_sigma[i-slo] = Real(0.0);
            p_sigma_cumsum[i-slo] = Real(0.0);
            if (i <= ohi+1) {
                p_sigma_star[i-sslo] = Real(0.0);
                p_sigma_star_cumsum[i-sslo] = Real(0.0);
            }
        });
    }
}

SigmaBox::SigmaBox (const Box& box, const BoxArray& grids, const Real* dx, int ncell, int delta)
{
    BL_ASSERT(box.cellCentered());

    const IntVect& sz = box.size();
    const int*     lo = box.loVect();
    const int*     hi = box.hiVect();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        sigma                [idim].resize(sz[idim]+1);
        sigma_cumsum         [idim].resize(sz[idim]+1);
        sigma_star           [idim].resize(sz[idim]+1);
        sigma_star_cumsum    [idim].resize(sz[idim]+1);
        sigma_fac            [idim].resize(sz[idim]+1);
        sigma_cumsum_fac     [idim].resize(sz[idim]+1);
        sigma_star_fac       [idim].resize(sz[idim]+1);
        sigma_star_cumsum_fac[idim].resize(sz[idim]+1);

        sigma                [idim].m_lo = lo[idim];
        sigma                [idim].m_hi = hi[idim]+1;
        sigma_cumsum         [idim].m_lo = lo[idim];
        sigma_cumsum         [idim].m_hi = hi[idim]+1;
        sigma_star           [idim].m_lo = lo[idim];
        sigma_star           [idim].m_hi = hi[idim]+1;
        sigma_star_cumsum    [idim].m_lo = lo[idim];
        sigma_star_cumsum    [idim].m_hi = hi[idim]+1;
        sigma_fac            [idim].m_lo = lo[idim];
        sigma_fac            [idim].m_hi = hi[idim]+1;
        sigma_cumsum_fac     [idim].m_lo = lo[idim];
        sigma_cumsum_fac     [idim].m_hi = hi[idim]+1;
        sigma_star_fac       [idim].m_lo = lo[idim];
        sigma_star_fac       [idim].m_hi = hi[idim]+1;
        sigma_star_cumsum_fac[idim].m_lo = lo[idim];
        sigma_star_cumsum_fac[idim].m_hi = hi[idim]+1;
    }

    Array<Real,AMREX_SPACEDIM> fac;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fac[idim] = 4.0_rt*PhysConst::c/(dx[idim]*static_cast<Real>(delta*delta));
    }

    const std::vector<std::pair<int,Box> >& isects = grids.intersections(box, false, ncell);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        int jdim = (idim+1) % AMREX_SPACEDIM;
#if (AMREX_SPACEDIM == 3)
        int kdim = (idim+2) % AMREX_SPACEDIM;
#endif

        Vector<int> direct_faces, side_faces, direct_side_edges, side_side_edges, corners;
        for (const auto& kv : isects)
        {
            const Box& grid_box = grids[kv.first];

            if (amrex::grow(grid_box, idim, ncell).intersects(box))
            {
                direct_faces.push_back(kv.first);
            }
            else if (amrex::grow(grid_box, jdim, ncell).intersects(box))
            {
                side_faces.push_back(kv.first);
            }
#if (AMREX_SPACEDIM == 3)
            else if (amrex::grow(grid_box, kdim, ncell).intersects(box))
            {
                side_faces.push_back(kv.first);
            }
            else if (amrex::grow(amrex::grow(grid_box,idim,ncell),
                                 jdim,ncell).intersects(box))
            {
                direct_side_edges.push_back(kv.first);
            }
            else if (amrex::grow(amrex::grow(grid_box,idim,ncell),
                                 kdim,ncell).intersects(box))
            {
                direct_side_edges.push_back(kv.first);
            }
            else if (amrex::grow(amrex::grow(grid_box,jdim,ncell),
                                 kdim,ncell).intersects(box))
            {
                side_side_edges.push_back(kv.first);
            }
#endif
            else
            {
                corners.push_back(kv.first);
            }
        }

        for (auto gid : corners)
        {
            const Box& grid_box = grids[gid];

            Box lobox = amrex::adjCellLo(grid_box, idim, ncell);
            lobox.grow(jdim,ncell);
#if (AMREX_SPACEDIM == 3)
            lobox.grow(kdim,ncell);
#endif
            Box looverlap = lobox & box;
            if (looverlap.ok()) {
                FillLo(idim, sigma[idim], sigma_cumsum[idim],
                       sigma_star[idim], sigma_star_cumsum[idim],
                       looverlap, grid_box, fac[idim]);
            }

            Box hibox = amrex::adjCellHi(grid_box, idim, ncell);
            hibox.grow(jdim,ncell);
#if (AMREX_SPACEDIM == 3)
            hibox.grow(kdim,ncell);
#endif
            Box hioverlap = hibox & box;
            if (hioverlap.ok()) {
                FillHi(idim, sigma[idim], sigma_cumsum[idim],
                       sigma_star[idim],  sigma_star_cumsum[idim],
                       hioverlap, grid_box, fac[idim]);
            }

            if (!looverlap.ok() && !hioverlap.ok()) {
                amrex::Abort("SigmaBox::SigmaBox(): corners, how did this happen?\n");
            }
        }

#if (AMREX_SPACEDIM == 3)
        for (auto gid : side_side_edges)
        {
            const Box& grid_box = grids[gid];
            const Box& overlap = amrex::grow(amrex::grow(grid_box,jdim,ncell),kdim,ncell) & box;
            if (overlap.ok()) {
                FillZero(idim, sigma[idim], sigma_cumsum[idim],
                        sigma_star[idim], sigma_star_cumsum[idim], overlap);
            }
            else {
                amrex::Abort("SigmaBox::SigmaBox(): side_side_edges, how did this happen?\n");
            }
        }

        for (auto gid : direct_side_edges)
        {
            const Box& grid_box = grids[gid];

            Box lobox = amrex::adjCellLo(grid_box, idim, ncell);
            Box looverlap = lobox.grow(jdim,ncell).grow(kdim,ncell) & box;
            if (looverlap.ok()) {
                FillLo(idim, sigma[idim], sigma_cumsum[idim],
                      sigma_star[idim],  sigma_star_cumsum[idim],
                      looverlap, grid_box, fac[idim]);
            }

            Box hibox = amrex::adjCellHi(grid_box, idim, ncell);
            Box hioverlap = hibox.grow(jdim,ncell).grow(kdim,ncell) & box;
            if (hioverlap.ok()) {
                FillHi(idim, sigma[idim], sigma_cumsum[idim],
                      sigma_star[idim],  sigma_star_cumsum[idim],
                      hioverlap, grid_box, fac[idim]);
            }

            if (!looverlap.ok() && !hioverlap.ok()) {
                amrex::Abort("SigmaBox::SigmaBox(): direct_side_edges, how did this happen?\n");
            }
        }
#endif

        for (auto gid : side_faces)
        {
            const Box& grid_box = grids[gid];
#if (AMREX_SPACEDIM == 2)
            const Box& overlap = amrex::grow(grid_box,jdim,ncell) & box;
#else
            const Box& overlap = amrex::grow(amrex::grow(grid_box,jdim,ncell),kdim,ncell) & box;
#endif
            if (overlap.ok()) {
                FillZero(idim, sigma[idim], sigma_cumsum[idim],
                        sigma_star[idim], sigma_star_cumsum[idim], overlap);
            } else {
                amrex::Abort("SigmaBox::SigmaBox(): side_faces, how did this happen?\n");
            }
        }

        for (auto gid : direct_faces)
        {
            const Box& grid_box = grids[gid];

            const Box& lobox = amrex::adjCellLo(grid_box, idim, ncell);
            Box looverlap = lobox & box;
            if (looverlap.ok()) {
                FillLo(idim, sigma[idim], sigma_cumsum[idim],
                      sigma_star[idim],  sigma_star_cumsum[idim],
                      looverlap, grid_box, fac[idim]);
            }

            const Box& hibox = amrex::adjCellHi(grid_box, idim, ncell);
            Box hioverlap = hibox & box;
            if (hioverlap.ok()) {
                FillHi(idim, sigma[idim], sigma_cumsum[idim],
                      sigma_star[idim],  sigma_star_cumsum[idim],
                      hioverlap, grid_box, fac[idim]);
            }

            if (!looverlap.ok() && !hioverlap.ok()) {
                amrex::Abort("SigmaBox::SigmaBox(): direct faces, how did this happen?\n");
            }
        }

        if (direct_faces.size() > 1) {
            amrex::Abort("SigmaBox::SigmaBox(): direct_faces.size() > 1, Box gaps not wide enough?\n");
        }
    }

    amrex::Gpu::synchronize();
}


void
SigmaBox::ComputePMLFactorsB (const Real* a_dx, Real dt)
{
    GpuArray<Real*,AMREX_SPACEDIM> p_sigma_star_fac;
    GpuArray<Real*,AMREX_SPACEDIM> p_sigma_star_cumsum_fac;
    GpuArray<Real const*,AMREX_SPACEDIM> p_sigma_star;
    GpuArray<Real const*,AMREX_SPACEDIM> p_sigma_star_cumsum;
    GpuArray<int, AMREX_SPACEDIM> N;
    GpuArray<Real, AMREX_SPACEDIM> dx;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        p_sigma_star_fac[idim] = sigma_star_fac[idim].data();
        p_sigma_star_cumsum_fac[idim] = sigma_star_cumsum_fac[idim].data();
        p_sigma_star[idim] = sigma_star[idim].data();
        p_sigma_star_cumsum[idim] = sigma_star_cumsum[idim].data();
        N[idim] = sigma_star[idim].size();
        dx[idim] = a_dx[idim];
    }
    amrex::ParallelFor(amrex::max(AMREX_D_DECL(N[0],N[1],N[2])),
    [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (i < N[idim]) {
                p_sigma_star_fac[idim][i] = std::exp(-p_sigma_star[idim][i]*dt);
                p_sigma_star_cumsum_fac[idim][i] = std::exp(-p_sigma_star_cumsum[idim][i]*dx[idim]);
            }
        }
    });
}

void
SigmaBox::ComputePMLFactorsE (const Real* a_dx, Real dt)
{
    GpuArray<Real*,AMREX_SPACEDIM> p_sigma_fac;
    GpuArray<Real*,AMREX_SPACEDIM> p_sigma_cumsum_fac;
    GpuArray<Real const*,AMREX_SPACEDIM> p_sigma;
    GpuArray<Real const*,AMREX_SPACEDIM> p_sigma_cumsum;
    GpuArray<int, AMREX_SPACEDIM> N;
    GpuArray<Real, AMREX_SPACEDIM> dx;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        p_sigma_fac[idim] = sigma_fac[idim].data();
        p_sigma_cumsum_fac[idim] = sigma_cumsum_fac[idim].data();
        p_sigma[idim] = sigma[idim].data();
        p_sigma_cumsum[idim] = sigma_cumsum[idim].data();
        N[idim] = sigma[idim].size();
        dx[idim] = a_dx[idim];
    }
    amrex::ParallelFor(amrex::max(AMREX_D_DECL(N[0],N[1],N[2])),
    [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (i < N[idim]) {
                p_sigma_fac[idim][i] = std::exp(-p_sigma[idim][i]*dt);
                p_sigma_cumsum_fac[idim][i] = std::exp(-p_sigma_cumsum[idim][i]*dx[idim]);
            }
        }
    });
}

MultiSigmaBox::MultiSigmaBox (const BoxArray& ba, const DistributionMapping& dm,
                              const BoxArray& grid_ba, const Real* dx, int ncell, int delta)
    : FabArray<SigmaBox>(ba,dm,1,0,MFInfo(),
                         SigmaBoxFactory(grid_ba,dx,ncell,delta))
{}

void
MultiSigmaBox::ComputePMLFactorsB (const Real* dx, Real dt)
{
    if (dt == dt_B) return;

    dt_B = dt;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        (*this)[mfi].ComputePMLFactorsB(dx, dt);
    }
}

void
MultiSigmaBox::ComputePMLFactorsE (const Real* dx, Real dt)
{
    if (dt == dt_E) return;

    dt_E = dt;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        (*this)[mfi].ComputePMLFactorsE(dx, dt);
    }
}

PML::PML (const int lev, const BoxArray& grid_ba, const DistributionMapping& /*grid_dm*/,
          const Geometry* geom, const Geometry* cgeom,
          int ncell, int delta, amrex::IntVect ref_ratio,
          Real dt, int nox_fft, int noy_fft, int noz_fft, bool do_nodal,
          int do_moving_window, int /*pml_has_particles*/, int do_pml_in_domain,
          const bool J_linear_in_time,
          const bool do_pml_dive_cleaning, const bool do_pml_divb_cleaning,
          const amrex::IntVect do_pml_Lo, const amrex::IntVect do_pml_Hi)
    : m_dive_cleaning(do_pml_dive_cleaning),
      m_divb_cleaning(do_pml_divb_cleaning),
      m_geom(geom),
      m_cgeom(cgeom)
{
    // When `do_pml_in_domain` is true, the PML overlap with the last `ncell` of the physical domain
    // (instead of extending `ncell` outside of the physical domain)
    // In order to implement this, a reduced domain is created here (decreased by ncells in all direction)
    // and passed to `MakeBoxArray`, which surrounds it by PML boxes
    // (thus creating the PML boxes at the right position, where they overlap with the original domain)
    // minimalBox provides the bounding box around grid_ba for level, lev.
    // Note that this is okay to build pml inside domain for a single patch, or joint patches
    // with same [min,max]. But it does not support multiple disjoint refinement patches.
    Box domain0 = grid_ba.minimalBox();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (do_pml_Lo[idim]){
            domain0.growLo(idim, -ncell);
        }
        if (do_pml_Hi[idim]){
            domain0.growHi(idim, -ncell);
        }
    }
    const BoxArray grid_ba_reduced = BoxArray(grid_ba.boxList().intersect(domain0));

    const BoxArray& ba = (do_pml_in_domain)?
          MakeBoxArray(*geom, grid_ba_reduced, ncell, do_pml_in_domain, do_pml_Lo, do_pml_Hi) :
          MakeBoxArray(*geom, grid_ba, ncell, do_pml_in_domain, do_pml_Lo, do_pml_Hi);
    if (ba.empty()) {
        m_ok = false;
        return;
    } else {
        m_ok = true;
    }

    DistributionMapping dm{ba};

    // Define the number of guard cells in each direction, for E, B, and F
    IntVect nge = IntVect(AMREX_D_DECL(2, 2, 2));
    IntVect ngb = IntVect(AMREX_D_DECL(2, 2, 2));
    int ngf_int = 0;
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::CKC) ngf_int = std::max( ngf_int, 1 );
    IntVect ngf = IntVect(AMREX_D_DECL(ngf_int, ngf_int, ngf_int));

    if (do_moving_window) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev <= 1,
            "The number of grow cells for the moving window currently assumes 2 levels max.");
        int rr = ref_ratio[WarpX::moving_window_dir];
        nge[WarpX::moving_window_dir] = std::max(nge[WarpX::moving_window_dir], rr);
        ngb[WarpX::moving_window_dir] = std::max(ngb[WarpX::moving_window_dir], rr);
        ngf[WarpX::moving_window_dir] = std::max(ngf[WarpX::moving_window_dir], rr);
    }

    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        // Increase the number of guard cells, in order to fit the extent
        // of the stencil for the spectral solver
        int ngFFt_x = do_nodal ? nox_fft : nox_fft/2;
        int ngFFt_y = do_nodal ? noy_fft : noy_fft/2;
        int ngFFt_z = do_nodal ? noz_fft : noz_fft/2;

        ParmParse pp_psatd("psatd");
        queryWithParser(pp_psatd, "nx_guard", ngFFt_x);
        queryWithParser(pp_psatd, "ny_guard", ngFFt_y);
        queryWithParser(pp_psatd, "nz_guard", ngFFt_z);

#if (AMREX_SPACEDIM == 3)
        IntVect ngFFT = IntVect(ngFFt_x, ngFFt_y, ngFFt_z);
#elif (AMREX_SPACEDIM == 2)
        IntVect ngFFT = IntVect(ngFFt_x, ngFFt_z);
#endif

        // Set the number of guard cells to the maximum of each field
        // (all fields should have the same number of guard cells)
        ngFFT = ngFFT.max(nge);
        ngFFT = ngFFT.max(ngb);
        ngFFT = ngFFT.max(ngf);
        nge = ngFFT;
        ngb = ngFFT;
        ngf = ngFFT;
    }

    // Allocate diagonal components (xx,yy,zz) only with divergence cleaning
    const int ncompe = (m_dive_cleaning) ? 3 : 2;
    const int ncompb = (m_divb_cleaning) ? 3 : 2;

    pml_E_fp[0] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getEfield_fp(0,0).ixType().toIntVect() ), dm, ncompe, nge );
    pml_E_fp[1] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getEfield_fp(0,1).ixType().toIntVect() ), dm, ncompe, nge );
    pml_E_fp[2] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getEfield_fp(0,2).ixType().toIntVect() ), dm, ncompe, nge );

    pml_B_fp[0] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getBfield_fp(0,0).ixType().toIntVect() ), dm, ncompb, ngb );
    pml_B_fp[1] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getBfield_fp(0,1).ixType().toIntVect() ), dm, ncompb, ngb );
    pml_B_fp[2] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getBfield_fp(0,2).ixType().toIntVect() ), dm, ncompb, ngb );

    pml_E_fp[0]->setVal(0.0);
    pml_E_fp[1]->setVal(0.0);
    pml_E_fp[2]->setVal(0.0);
    pml_B_fp[0]->setVal(0.0);
    pml_B_fp[1]->setVal(0.0);
    pml_B_fp[2]->setVal(0.0);

    pml_j_fp[0] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getcurrent_fp(0,0).ixType().toIntVect() ), dm, 1, ngb );
    pml_j_fp[1] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getcurrent_fp(0,1).ixType().toIntVect() ), dm, 1, ngb );
    pml_j_fp[2] = std::make_unique<MultiFab>(amrex::convert( ba,
        WarpX::GetInstance().getcurrent_fp(0,2).ixType().toIntVect() ), dm, 1, ngb );

    pml_j_fp[0]->setVal(0.0);
    pml_j_fp[1]->setVal(0.0);
    pml_j_fp[2]->setVal(0.0);

    if (m_dive_cleaning)
    {
        const amrex::IntVect& F_nodal_flag = amrex::IntVect::TheNodeVector();
        pml_F_fp = std::make_unique<MultiFab>(amrex::convert(ba, F_nodal_flag), dm, 3, ngf);
        pml_F_fp->setVal(0.0);
    }

    if (m_divb_cleaning)
    {
        // TODO Shall we define a separate guard cells parameter ngG?
        const amrex::IntVect& G_nodal_flag = (do_nodal) ? amrex::IntVect::TheNodeVector()
                                                        : amrex::IntVect::TheCellVector();
        pml_G_fp = std::make_unique<MultiFab>(amrex::convert(ba, G_nodal_flag), dm, 3, ngf);
        pml_G_fp->setVal(0.0);
    }

    if (do_pml_in_domain){
        sigba_fp = std::make_unique<MultiSigmaBox>(ba, dm, grid_ba_reduced, geom->CellSize(), ncell, delta);
    }
    else {
        sigba_fp = std::make_unique<MultiSigmaBox>(ba, dm, grid_ba, geom->CellSize(), ncell, delta);
    }

    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
#ifndef WARPX_USE_PSATD
        amrex::ignore_unused(lev, dt, J_linear_in_time);
#   if(AMREX_SPACEDIM!=3)
        amrex::ignore_unused(noy_fft);
#   endif
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                                         "PML: PSATD solver selected but not built.");
#else
        // Flags passed to the spectral solver constructor
        const bool in_pml = true;
        const bool periodic_single_box = false;
        const bool update_with_rho = false;
        const bool fft_do_time_averaging = false;
        const RealVect dx{AMREX_D_DECL(geom->CellSize(0), geom->CellSize(1), geom->CellSize(2))};
        // Get the cell-centered box, with guard cells
        BoxArray realspace_ba = ba; // Copy box
        amrex::Vector<amrex::Real> const v_galilean_zero = {0., 0., 0.};
        amrex::Vector<amrex::Real> const v_comoving_zero = {0., 0., 0.};
        realspace_ba.enclosedCells().grow(nge); // cell-centered + guard cells
        spectral_solver_fp = std::make_unique<SpectralSolver>(lev, realspace_ba, dm,
            nox_fft, noy_fft, noz_fft, do_nodal, WarpX::fill_guards, v_galilean_zero,
            v_comoving_zero, dx, dt, in_pml, periodic_single_box, update_with_rho,
            fft_do_time_averaging, J_linear_in_time, m_dive_cleaning, m_divb_cleaning);
#endif
    }

    if (cgeom)
    {
        if (WarpX::maxwell_solver_id != MaxwellSolverAlgo::PSATD) {
            nge = IntVect(AMREX_D_DECL(1, 1, 1));
            ngb = IntVect(AMREX_D_DECL(1, 1, 1));
        }

        BoxArray grid_cba = grid_ba;
        grid_cba.coarsen(ref_ratio);

        // assuming that the bounding box around grid_cba is a single patch, and not disjoint patches, similar to fine patch.
        amrex::Box domain1 = grid_cba.minimalBox();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (do_pml_Lo[idim]){
                // ncell is divided by refinement ratio to ensure that the
                // physical width of the PML region is equal in fine and coarse patch
                domain1.growLo(idim, -ncell/ref_ratio[idim]);
            }
            if (do_pml_Hi[idim]){
                // ncell is divided by refinement ratio to ensure that the
                // physical width of the PML region is equal in fine and coarse patch
                domain1.growHi(idim, -ncell/ref_ratio[idim]);
            }
        }
        const BoxArray grid_cba_reduced = BoxArray(grid_cba.boxList().intersect(domain1));

        // Assuming that refinement ratio is equal in all dimensions
        const BoxArray& cba = (do_pml_in_domain) ?
            MakeBoxArray(*cgeom, grid_cba_reduced, ncell/ref_ratio[0], do_pml_in_domain, do_pml_Lo, do_pml_Hi) :
            MakeBoxArray(*cgeom, grid_cba, ncell, do_pml_in_domain, do_pml_Lo, do_pml_Hi);
        DistributionMapping cdm{cba};

        pml_E_cp[0] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getEfield_cp(1,0).ixType().toIntVect() ), cdm, ncompe, nge );
        pml_E_cp[1] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getEfield_cp(1,1).ixType().toIntVect() ), cdm, ncompe, nge );
        pml_E_cp[2] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getEfield_cp(1,2).ixType().toIntVect() ), cdm, ncompe, nge );

        pml_B_cp[0] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getBfield_cp(1,0).ixType().toIntVect() ), cdm, ncompb, ngb );
        pml_B_cp[1] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getBfield_cp(1,1).ixType().toIntVect() ), cdm, ncompb, ngb );
        pml_B_cp[2] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getBfield_cp(1,2).ixType().toIntVect() ), cdm, ncompb, ngb );

        pml_E_cp[0]->setVal(0.0);
        pml_E_cp[1]->setVal(0.0);
        pml_E_cp[2]->setVal(0.0);
        pml_B_cp[0]->setVal(0.0);
        pml_B_cp[1]->setVal(0.0);
        pml_B_cp[2]->setVal(0.0);

        if (m_dive_cleaning)
        {
            const amrex::IntVect& F_nodal_flag = amrex::IntVect::TheNodeVector();
            pml_F_cp = std::make_unique<MultiFab>(amrex::convert(cba, F_nodal_flag), cdm, 3, ngf);
            pml_F_cp->setVal(0.0);
        }

        if (m_divb_cleaning)
        {
            // TODO Shall we define a separate guard cells parameter ngG?
            const amrex::IntVect& G_nodal_flag = (do_nodal) ? amrex::IntVect::TheNodeVector()
                                                            : amrex::IntVect::TheCellVector();
            pml_G_cp = std::make_unique<MultiFab>(amrex::convert(cba, G_nodal_flag), cdm, 3, ngf);
            pml_G_cp->setVal(0.0);
        }

        pml_j_cp[0] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getcurrent_cp(1,0).ixType().toIntVect() ), cdm, 1, ngb );
        pml_j_cp[1] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getcurrent_cp(1,1).ixType().toIntVect() ), cdm, 1, ngb );
        pml_j_cp[2] = std::make_unique<MultiFab>(amrex::convert( cba,
            WarpX::GetInstance().getcurrent_cp(1,2).ixType().toIntVect() ), cdm, 1, ngb );

        pml_j_cp[0]->setVal(0.0);
        pml_j_cp[1]->setVal(0.0);
        pml_j_cp[2]->setVal(0.0);

        if (do_pml_in_domain){
            // Note - assuming that the refinement ratio is equal in all dimensions
            sigba_cp = std::make_unique<MultiSigmaBox>(cba, cdm, grid_cba_reduced, cgeom->CellSize(), ncell/ref_ratio[0], delta/ref_ratio[0]);
        } else {
            sigba_cp = std::make_unique<MultiSigmaBox>(cba, cdm, grid_cba, cgeom->CellSize(), ncell, delta);
        }

        if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
#ifndef WARPX_USE_PSATD
            amrex::ignore_unused(dt);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
                                             "PML: PSATD solver selected but not built.");
#else
            // Flags passed to the spectral solver constructor
            const bool in_pml = true;
            const bool periodic_single_box = false;
            const bool update_with_rho = false;
            const bool fft_do_time_averaging = false;
            const RealVect cdx{AMREX_D_DECL(cgeom->CellSize(0), cgeom->CellSize(1), cgeom->CellSize(2))};
            // Get the cell-centered box, with guard cells
            BoxArray realspace_cba = cba; // Copy box
            amrex::Vector<amrex::Real> const v_galilean_zero = {0., 0., 0.};
            amrex::Vector<amrex::Real> const v_comoving_zero = {0., 0., 0.};
            realspace_cba.enclosedCells().grow(nge); // cell-centered + guard cells
            spectral_solver_cp = std::make_unique<SpectralSolver>(lev, realspace_cba, cdm,
                nox_fft, noy_fft, noz_fft, do_nodal, WarpX::fill_guards, v_galilean_zero,
                v_comoving_zero, cdx, dt, in_pml, periodic_single_box, update_with_rho,
                fft_do_time_averaging, J_linear_in_time, m_dive_cleaning, m_divb_cleaning);
#endif
        }
    }
}

BoxArray
PML::MakeBoxArray (const amrex::Geometry& geom, const amrex::BoxArray& grid_ba,
                   int ncell, int do_pml_in_domain,
                   const amrex::IntVect do_pml_Lo, const amrex::IntVect do_pml_Hi)
{
    Box domain = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (do_pml_Lo[idim]){
            domain.growLo(idim, ncell);
        }
        if (do_pml_Hi[idim]){
            domain.growHi(idim, ncell);
        }
    }
    BoxList bl;
    for (int i = 0, N = grid_ba.size(); i < N; ++i)
    {
        const Box& grid_bx = grid_ba[i];
        const IntVect& grid_bx_sz = grid_bx.size();

        if (do_pml_in_domain == 0) {
            // Make sure that, in the case of several distinct refinement patches,
            //  the PML cells surrounding these patches cannot overlap
            // The check is only needed along the axis where PMLs are being used.
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (do_pml_Lo[idim] || do_pml_Hi[idim]) {
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                        grid_bx.length(idim) > ncell,
                        "Consider using larger amr.blocking_factor with PMLs");
                }
            }
        }

        Box bx = grid_bx;
        bx.grow(ncell);
        bx &= domain;

        Vector<Box> bndryboxes;
#if (AMREX_SPACEDIM == 3)
        int kbegin = -1, kend = 1;
#else
        int kbegin =  0, kend = 0;
#endif
        for (int kk = kbegin; kk <= kend; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
                for (int ii = -1; ii <= 1; ++ii) {
                    if (ii != 0 || jj != 0 || kk != 0) {
                        Box b = grid_bx;
                        b.shift(grid_bx_sz * IntVect{AMREX_D_DECL(ii,jj,kk)});
                        b &= bx;
                        if (b.ok()) {
                            bndryboxes.push_back(b);
                        }
                    }
                }
            }
        }

        const BoxList& noncovered = grid_ba.complementIn(bx);
        for (const Box& b : noncovered) {
            for (const auto& bb : bndryboxes) {
                Box ib = b & bb;
                if (ib.ok()) {
                    bl.push_back(ib);
                }
            }
        }
    }

    BoxArray ba(bl);
    ba.removeOverlap(false);

    return ba;
}

void
PML::ComputePMLFactors (amrex::Real dt)
{
    if (sigba_fp) {
        sigba_fp->ComputePMLFactorsB(m_geom->CellSize(), dt);
        sigba_fp->ComputePMLFactorsE(m_geom->CellSize(), dt);
    }
    if (sigba_cp) {
        sigba_cp->ComputePMLFactorsB(m_cgeom->CellSize(), dt);
        sigba_cp->ComputePMLFactorsE(m_cgeom->CellSize(), dt);
    }
}

std::array<MultiFab*,3>
PML::GetE_fp ()
{
    return {pml_E_fp[0].get(), pml_E_fp[1].get(), pml_E_fp[2].get()};
}

std::array<MultiFab*,3>
PML::GetB_fp ()
{
    return {pml_B_fp[0].get(), pml_B_fp[1].get(), pml_B_fp[2].get()};
}

std::array<MultiFab*,3>
PML::Getj_fp ()
{
    return {pml_j_fp[0].get(), pml_j_fp[1].get(), pml_j_fp[2].get()};
}

std::array<MultiFab*,3>
PML::GetE_cp ()
{
    return {pml_E_cp[0].get(), pml_E_cp[1].get(), pml_E_cp[2].get()};
}

std::array<MultiFab*,3>
PML::GetB_cp ()
{
    return {pml_B_cp[0].get(), pml_B_cp[1].get(), pml_B_cp[2].get()};
}

std::array<MultiFab*,3>
PML::Getj_cp ()
{
    return {pml_j_cp[0].get(), pml_j_cp[1].get(), pml_j_cp[2].get()};
}

MultiFab*
PML::GetF_fp ()
{
    return pml_F_fp.get();
}

MultiFab*
PML::GetF_cp ()
{
    return pml_F_cp.get();
}

MultiFab*
PML::GetG_fp ()
{
    return pml_G_fp.get();
}

MultiFab*
PML::GetG_cp ()
{
    return pml_G_cp.get();
}

void
PML::ExchangeB (const std::array<amrex::MultiFab*,3>& B_fp,
                const std::array<amrex::MultiFab*,3>& B_cp,
                int do_pml_in_domain)
{
  ExchangeB(PatchType::fine, B_fp, do_pml_in_domain);
  ExchangeB(PatchType::coarse, B_cp, do_pml_in_domain);
}

void
PML::ExchangeB (PatchType patch_type,
                const std::array<amrex::MultiFab*,3>& Bp,
                int do_pml_in_domain)
{
    if (patch_type == PatchType::fine && pml_B_fp[0] && Bp[0])
    {
        Exchange(*pml_B_fp[0], *Bp[0], *m_geom, do_pml_in_domain);
        Exchange(*pml_B_fp[1], *Bp[1], *m_geom, do_pml_in_domain);
        Exchange(*pml_B_fp[2], *Bp[2], *m_geom, do_pml_in_domain);
    }
    else if (patch_type == PatchType::coarse && pml_B_cp[0] && Bp[0])
    {
        Exchange(*pml_B_cp[0], *Bp[0], *m_cgeom, do_pml_in_domain);
        Exchange(*pml_B_cp[1], *Bp[1], *m_cgeom, do_pml_in_domain);
        Exchange(*pml_B_cp[2], *Bp[2], *m_cgeom, do_pml_in_domain);
    }
}

void
PML::ExchangeE (const std::array<amrex::MultiFab*,3>& E_fp,
                const std::array<amrex::MultiFab*,3>& E_cp,
                int do_pml_in_domain)
{
    ExchangeE(PatchType::fine, E_fp, do_pml_in_domain);
    ExchangeE(PatchType::coarse, E_cp, do_pml_in_domain);
}

void
PML::ExchangeE (PatchType patch_type,
                const std::array<amrex::MultiFab*,3>& Ep,
                int do_pml_in_domain)
{
    if (patch_type == PatchType::fine && pml_E_fp[0] && Ep[0])
    {
        Exchange(*pml_E_fp[0], *Ep[0], *m_geom, do_pml_in_domain);
        Exchange(*pml_E_fp[1], *Ep[1], *m_geom, do_pml_in_domain);
        Exchange(*pml_E_fp[2], *Ep[2], *m_geom, do_pml_in_domain);
    }
    else if (patch_type == PatchType::coarse && pml_E_cp[0] && Ep[0])
    {
        Exchange(*pml_E_cp[0], *Ep[0], *m_cgeom, do_pml_in_domain);
        Exchange(*pml_E_cp[1], *Ep[1], *m_cgeom, do_pml_in_domain);
        Exchange(*pml_E_cp[2], *Ep[2], *m_cgeom, do_pml_in_domain);
    }
}

void
PML::CopyJtoPMLs (PatchType patch_type,
                const std::array<amrex::MultiFab*,3>& jp)
{
    if (patch_type == PatchType::fine && pml_j_fp[0] && jp[0])
    {
        CopyToPML(*pml_j_fp[0], *jp[0], *m_geom);
        CopyToPML(*pml_j_fp[1], *jp[1], *m_geom);
        CopyToPML(*pml_j_fp[2], *jp[2], *m_geom);
    }
    else if (patch_type == PatchType::coarse && pml_j_cp[0] && jp[0])
    {
        CopyToPML(*pml_j_cp[0], *jp[0], *m_cgeom);
        CopyToPML(*pml_j_cp[1], *jp[1], *m_cgeom);
        CopyToPML(*pml_j_cp[2], *jp[2], *m_cgeom);
    }
}

void
PML::CopyJtoPMLs (const std::array<amrex::MultiFab*,3>& j_fp,
                const std::array<amrex::MultiFab*,3>& j_cp)
{
    CopyJtoPMLs(PatchType::fine, j_fp);
    CopyJtoPMLs(PatchType::coarse, j_cp);
}

void
PML::ExchangeF (amrex::MultiFab* F_fp, amrex::MultiFab* F_cp, int do_pml_in_domain)
{
    ExchangeF(PatchType::fine, F_fp, do_pml_in_domain);
    ExchangeF(PatchType::coarse, F_cp, do_pml_in_domain);
}

void
PML::ExchangeF (PatchType patch_type, amrex::MultiFab* Fp, int do_pml_in_domain)
{
    if (patch_type == PatchType::fine && pml_F_fp && Fp) {
        Exchange(*pml_F_fp, *Fp, *m_geom, do_pml_in_domain);
    } else if (patch_type == PatchType::coarse && pml_F_cp && Fp) {
        Exchange(*pml_F_cp, *Fp, *m_cgeom, do_pml_in_domain);
    }
}

void PML::ExchangeG (amrex::MultiFab* G_fp, amrex::MultiFab* G_cp, int do_pml_in_domain)
{
    ExchangeG(PatchType::fine, G_fp, do_pml_in_domain);
    ExchangeG(PatchType::coarse, G_cp, do_pml_in_domain);
}

void PML::ExchangeG (PatchType patch_type, amrex::MultiFab* Gp, int do_pml_in_domain)
{
    if (patch_type == PatchType::fine && pml_G_fp && Gp)
    {
        Exchange(*pml_G_fp, *Gp, *m_geom, do_pml_in_domain);
    }
    else if (patch_type == PatchType::coarse && pml_G_cp && Gp)
    {
        Exchange(*pml_G_cp, *Gp, *m_cgeom, do_pml_in_domain);
    }
}

void
PML::Exchange (MultiFab& pml, MultiFab& reg, const Geometry& geom,
                int do_pml_in_domain)
{
    WARPX_PROFILE("PML::Exchange");

    const IntVect& ngr = reg.nGrowVect();
    const IntVect& ngp = pml.nGrowVect();
    const int ncp = pml.nComp();
    const auto& period = geom.periodicity();

    // Create temporary MultiFab to copy to and from the PML
    MultiFab tmpregmf(reg.boxArray(), reg.DistributionMap(), ncp, ngr);
    tmpregmf.setVal(0.0);

    // Create the sum of the split fields, in the PML
    MultiFab totpmlmf(pml.boxArray(), pml.DistributionMap(), 1, 0); // Allocate
    MultiFab::LinComb(totpmlmf, 1.0, pml, 0, 1.0, pml, 1, 0, 1, 0); // Sum
    if (ncp == 3) {
        MultiFab::Add(totpmlmf,pml,2,0,1,0); // Sum the third split component
    }

    // Copy from the sum of PML split field to valid cells of regular grid
    if (do_pml_in_domain){
        // Valid cells of the PML and of the regular grid overlap
        // Copy from valid cells of the PML to valid cells of the regular grid
        WarpXCommUtil::ParallelCopy(reg, totpmlmf, 0, 0, 1, IntVect(0), IntVect(0), period);
    } else {
        // Valid cells of the PML only overlap with guard cells of regular grid
        // (and outermost valid cell of the regular grid, for nodal direction)
        // Copy from valid cells of PML to ghost cells of regular grid
        // but avoid updating the outermost valid cell
        if (ngr.max() > 0) {
            MultiFab::Copy(tmpregmf, reg, 0, 0, 1, ngr);
            WarpXCommUtil::ParallelCopy(tmpregmf, totpmlmf, 0, 0, 1, IntVect(0), ngr, period);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(reg); mfi.isValid(); ++mfi)
            {
                const FArrayBox& src = tmpregmf[mfi];
                FArrayBox& dst = reg[mfi];
                const auto srcarr = src.array();
                auto dstarr = dst.array();
                const BoxList& bl = amrex::boxDiff(dst.box(), mfi.validbox());
                // boxDiff avoids the outermost valid cell
                for (const Box& bx : bl) {
                    amrex::ParallelFor(bx,
                                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                                       {
                                           dstarr(i,j,k,0) = srcarr(i,j,k,0);
                                       });
                }
            }
        }
    }

    // Copy from valid cells of the regular grid to guard cells of the PML
    // (and outermost valid cell in the nodal direction)
    // More specifically, copy from regular data to PML's first component
    // Zero out the second (and third) component
    MultiFab::Copy(tmpregmf,reg,0,0,1,0); // Fill first component of tmpregmf
    tmpregmf.setVal(0.0, 1, ncp-1, 0); // Zero out the second (and third) component
    if (do_pml_in_domain){
        // Where valid cells of tmpregmf overlap with PML valid cells,
        // copy the PML (this is order to avoid overwriting PML valid cells,
        // in the next `ParallelCopy`)
        WarpXCommUtil::ParallelCopy(tmpregmf, pml,0, 0, ncp, IntVect(0), IntVect(0), period);
    }
    WarpXCommUtil::ParallelCopy(pml, tmpregmf, 0, 0, ncp, IntVect(0), ngp, period);
}


void
PML::CopyToPML (MultiFab& pml, MultiFab& reg, const Geometry& geom)
{
  const IntVect& ngp = pml.nGrowVect();
  const auto& period = geom.periodicity();

  WarpXCommUtil::ParallelCopy(pml, reg, 0, 0, 1, IntVect(0), ngp, period);
}

void
PML::FillBoundary ()
{
    FillBoundaryE();
    FillBoundaryB();
    FillBoundaryF();
    FillBoundaryG();
}

void
PML::FillBoundaryE ()
{
    FillBoundaryE(PatchType::fine);
    FillBoundaryE(PatchType::coarse);
}

void
PML::FillBoundaryE (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_E_fp[0] && pml_E_fp[0]->nGrowVect().max() > 0)
    {
        const auto& period = m_geom->periodicity();
        Vector<MultiFab*> mf{pml_E_fp[0].get(),pml_E_fp[1].get(),pml_E_fp[2].get()};
        WarpXCommUtil::FillBoundary(mf, period);
    }
    else if (patch_type == PatchType::coarse && pml_E_cp[0] && pml_E_cp[0]->nGrowVect().max() > 0)
    {
        const auto& period = m_cgeom->periodicity();
        Vector<MultiFab*> mf{pml_E_cp[0].get(),pml_E_cp[1].get(),pml_E_cp[2].get()};
        WarpXCommUtil::FillBoundary(mf, period);
    }
}

void
PML::FillBoundaryB ()
{
    FillBoundaryB(PatchType::fine);
    FillBoundaryB(PatchType::coarse);
}

void
PML::FillBoundaryB (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_B_fp[0])
    {
        const auto& period = m_geom->periodicity();
        Vector<MultiFab*> mf{pml_B_fp[0].get(),pml_B_fp[1].get(),pml_B_fp[2].get()};
        WarpXCommUtil::FillBoundary(mf, period);
    }
    else if (patch_type == PatchType::coarse && pml_B_cp[0])
    {
        const auto& period = m_cgeom->periodicity();
        Vector<MultiFab*> mf{pml_B_cp[0].get(),pml_B_cp[1].get(),pml_B_cp[2].get()};
        WarpXCommUtil::FillBoundary(mf, period);
    }
}

void
PML::FillBoundaryF ()
{
    FillBoundaryF(PatchType::fine);
    FillBoundaryF(PatchType::coarse);
}

void
PML::FillBoundaryF (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_F_fp && pml_F_fp->nGrowVect().max() > 0)
    {
        const auto& period = m_geom->periodicity();
        WarpXCommUtil::FillBoundary(*pml_F_fp, period);
    }
    else if (patch_type == PatchType::coarse && pml_F_cp && pml_F_cp->nGrowVect().max() > 0)
    {
        const auto& period = m_cgeom->periodicity();
        WarpXCommUtil::FillBoundary(*pml_F_cp, period);
    }
}

void
PML::FillBoundaryG ()
{
    FillBoundaryG(PatchType::fine);
    FillBoundaryG(PatchType::coarse);
}

void
PML::FillBoundaryG (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_G_fp && pml_G_fp->nGrowVect().max() > 0)
    {
        const auto& period = m_geom->periodicity();
        WarpXCommUtil::FillBoundary(*pml_G_fp, period);
    }
    else if (patch_type == PatchType::coarse && pml_G_cp && pml_G_cp->nGrowVect().max() > 0)
    {
        const auto& period = m_cgeom->periodicity();
        WarpXCommUtil::FillBoundary(*pml_G_cp, period);
    }
}

void
PML::CheckPoint (const std::string& dir) const
{
    if (pml_E_fp[0])
    {
        VisMF::AsyncWrite(*pml_E_fp[0], dir+"_Ex_fp");
        VisMF::AsyncWrite(*pml_E_fp[1], dir+"_Ey_fp");
        VisMF::AsyncWrite(*pml_E_fp[2], dir+"_Ez_fp");
        VisMF::AsyncWrite(*pml_B_fp[0], dir+"_Bx_fp");
        VisMF::AsyncWrite(*pml_B_fp[1], dir+"_By_fp");
        VisMF::AsyncWrite(*pml_B_fp[2], dir+"_Bz_fp");
    }

    if (pml_E_cp[0])
    {
        VisMF::AsyncWrite(*pml_E_cp[0], dir+"_Ex_cp");
        VisMF::AsyncWrite(*pml_E_cp[1], dir+"_Ey_cp");
        VisMF::AsyncWrite(*pml_E_cp[2], dir+"_Ez_cp");
        VisMF::AsyncWrite(*pml_B_cp[0], dir+"_Bx_cp");
        VisMF::AsyncWrite(*pml_B_cp[1], dir+"_By_cp");
        VisMF::AsyncWrite(*pml_B_cp[2], dir+"_Bz_cp");
    }
}

void
PML::Restart (const std::string& dir)
{
    if (pml_E_fp[0])
    {
        VisMF::Read(*pml_E_fp[0], dir+"_Ex_fp");
        VisMF::Read(*pml_E_fp[1], dir+"_Ey_fp");
        VisMF::Read(*pml_E_fp[2], dir+"_Ez_fp");
        VisMF::Read(*pml_B_fp[0], dir+"_Bx_fp");
        VisMF::Read(*pml_B_fp[1], dir+"_By_fp");
        VisMF::Read(*pml_B_fp[2], dir+"_Bz_fp");
    }

    if (pml_E_cp[0])
    {
        VisMF::Read(*pml_E_cp[0], dir+"_Ex_cp");
        VisMF::Read(*pml_E_cp[1], dir+"_Ey_cp");
        VisMF::Read(*pml_E_cp[2], dir+"_Ez_cp");
        VisMF::Read(*pml_B_cp[0], dir+"_Bx_cp");
        VisMF::Read(*pml_B_cp[1], dir+"_By_cp");
        VisMF::Read(*pml_B_cp[2], dir+"_Bz_cp");
    }
}

#ifdef WARPX_USE_PSATD
void
PML::PushPSATD (const int lev) {

    // Update the fields on the fine and coarse patch
    PushPMLPSATDSinglePatch(lev, *spectral_solver_fp, pml_E_fp, pml_B_fp, pml_F_fp, pml_G_fp);
    if (spectral_solver_cp) {
        PushPMLPSATDSinglePatch(lev, *spectral_solver_cp, pml_E_cp, pml_B_cp, pml_F_cp, pml_G_cp);
    }
}

void
PushPMLPSATDSinglePatch (
    const int lev,
    SpectralSolver& solver,
    std::array<std::unique_ptr<amrex::MultiFab>,3>& pml_E,
    std::array<std::unique_ptr<amrex::MultiFab>,3>& pml_B,
    std::unique_ptr<amrex::MultiFab>& pml_F,
    std::unique_ptr<amrex::MultiFab>& pml_G)
{
    const SpectralFieldIndex& Idx = solver.m_spectral_index;

    // Perform forward Fourier transforms
    solver.ForwardTransform(lev, *pml_E[0], Idx.Exy, PMLComp::xy);
    solver.ForwardTransform(lev, *pml_E[0], Idx.Exz, PMLComp::xz);
    solver.ForwardTransform(lev, *pml_E[1], Idx.Eyx, PMLComp::yx);
    solver.ForwardTransform(lev, *pml_E[1], Idx.Eyz, PMLComp::yz);
    solver.ForwardTransform(lev, *pml_E[2], Idx.Ezx, PMLComp::zx);
    solver.ForwardTransform(lev, *pml_E[2], Idx.Ezy, PMLComp::zy);
    solver.ForwardTransform(lev, *pml_B[0], Idx.Bxy, PMLComp::xy);
    solver.ForwardTransform(lev, *pml_B[0], Idx.Bxz, PMLComp::xz);
    solver.ForwardTransform(lev, *pml_B[1], Idx.Byx, PMLComp::yx);
    solver.ForwardTransform(lev, *pml_B[1], Idx.Byz, PMLComp::yz);
    solver.ForwardTransform(lev, *pml_B[2], Idx.Bzx, PMLComp::zx);
    solver.ForwardTransform(lev, *pml_B[2], Idx.Bzy, PMLComp::zy);

    // WarpX::do_pml_dive_cleaning = true
    if (pml_F)
    {
        solver.ForwardTransform(lev, *pml_E[0], Idx.Exx, PMLComp::xx);
        solver.ForwardTransform(lev, *pml_E[1], Idx.Eyy, PMLComp::yy);
        solver.ForwardTransform(lev, *pml_E[2], Idx.Ezz, PMLComp::zz);
        solver.ForwardTransform(lev, *pml_F, Idx.Fx, PMLComp::x);
        solver.ForwardTransform(lev, *pml_F, Idx.Fy, PMLComp::y);
        solver.ForwardTransform(lev, *pml_F, Idx.Fz, PMLComp::z);
    }

    // WarpX::do_pml_divb_cleaning = true
    if (pml_G)
    {
        solver.ForwardTransform(lev, *pml_B[0], Idx.Bxx, PMLComp::xx);
        solver.ForwardTransform(lev, *pml_B[1], Idx.Byy, PMLComp::yy);
        solver.ForwardTransform(lev, *pml_B[2], Idx.Bzz, PMLComp::zz);
        solver.ForwardTransform(lev, *pml_G, Idx.Gx, PMLComp::x);
        solver.ForwardTransform(lev, *pml_G, Idx.Gy, PMLComp::y);
        solver.ForwardTransform(lev, *pml_G, Idx.Gz, PMLComp::z);
    }

    // Advance fields in spectral space
    solver.pushSpectralFields();

    // Perform backward Fourier transforms
    solver.BackwardTransform(lev, *pml_E[0], Idx.Exy, PMLComp::xy);
    solver.BackwardTransform(lev, *pml_E[0], Idx.Exz, PMLComp::xz);
    solver.BackwardTransform(lev, *pml_E[1], Idx.Eyx, PMLComp::yx);
    solver.BackwardTransform(lev, *pml_E[1], Idx.Eyz, PMLComp::yz);
    solver.BackwardTransform(lev, *pml_E[2], Idx.Ezx, PMLComp::zx);
    solver.BackwardTransform(lev, *pml_E[2], Idx.Ezy, PMLComp::zy);
    solver.BackwardTransform(lev, *pml_B[0], Idx.Bxy, PMLComp::xy);
    solver.BackwardTransform(lev, *pml_B[0], Idx.Bxz, PMLComp::xz);
    solver.BackwardTransform(lev, *pml_B[1], Idx.Byx, PMLComp::yx);
    solver.BackwardTransform(lev, *pml_B[1], Idx.Byz, PMLComp::yz);
    solver.BackwardTransform(lev, *pml_B[2], Idx.Bzx, PMLComp::zx);
    solver.BackwardTransform(lev, *pml_B[2], Idx.Bzy, PMLComp::zy);

    // WarpX::do_pml_dive_cleaning = true
    if (pml_F)
    {
        solver.BackwardTransform(lev, *pml_E[0], Idx.Exx, PMLComp::xx);
        solver.BackwardTransform(lev, *pml_E[1], Idx.Eyy, PMLComp::yy);
        solver.BackwardTransform(lev, *pml_E[2], Idx.Ezz, PMLComp::zz);
        solver.BackwardTransform(lev, *pml_F, Idx.Fx, PMLComp::x);
        solver.BackwardTransform(lev, *pml_F, Idx.Fy, PMLComp::y);
        solver.BackwardTransform(lev, *pml_F, Idx.Fz, PMLComp::z);
    }

    // WarpX::do_pml_divb_cleaning = true
    if (pml_G)
    {
        solver.BackwardTransform(lev, *pml_B[0], Idx.Bxx, PMLComp::xx);
        solver.BackwardTransform(lev, *pml_B[1], Idx.Byy, PMLComp::yy);
        solver.BackwardTransform(lev, *pml_B[2], Idx.Bzz, PMLComp::zz);
        solver.BackwardTransform(lev, *pml_G, Idx.Gx, PMLComp::x);
        solver.BackwardTransform(lev, *pml_G, Idx.Gy, PMLComp::y);
        solver.BackwardTransform(lev, *pml_G, Idx.Gz, PMLComp::z);
    }
}
#endif
