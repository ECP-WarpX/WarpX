/* Copyright 2021 David Grote
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PML_RZ.H"

#include "BoundaryConditions/PML_RZ.H"
#ifdef WARPX_USE_PSATD
#   include "FieldSolver/SpectralSolver/SpectralFieldDataRZ.H"
#endif
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/utils/Communication.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_RealVect.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <memory>

using namespace amrex;

PML_RZ::PML_RZ (const int lev, const amrex::BoxArray& grid_ba, const amrex::DistributionMapping& grid_dm,
                const amrex::Geometry* geom, const int ncell, const int do_pml_in_domain)
    : m_ncell(ncell),
      m_do_pml_in_domain(do_pml_in_domain),
      m_geom(geom)
{

    const amrex::MultiFab & Er_fp = WarpX::GetInstance().getEfield_fp(lev,0);
    const amrex::MultiFab & Et_fp = WarpX::GetInstance().getEfield_fp(lev,1);
    const amrex::BoxArray ba_Er = amrex::convert(grid_ba, Er_fp.ixType().toIntVect());
    const amrex::BoxArray ba_Et = amrex::convert(grid_ba, Et_fp.ixType().toIntVect());
    WarpX::AllocInitMultiFab(pml_E_fp[0], ba_Er, grid_dm, Er_fp.nComp(), Er_fp.nGrowVect(), "pml_E_fp[0]", 0.0_rt);
    WarpX::AllocInitMultiFab(pml_E_fp[1], ba_Et, grid_dm, Et_fp.nComp(), Et_fp.nGrowVect(), "pml_E_fp[1]", 0.0_rt);

    const amrex::MultiFab & Br_fp = WarpX::GetInstance().getBfield_fp(lev,0);
    const amrex::MultiFab & Bt_fp = WarpX::GetInstance().getBfield_fp(lev,1);
    const amrex::BoxArray ba_Br = amrex::convert(grid_ba, Br_fp.ixType().toIntVect());
    const amrex::BoxArray ba_Bt = amrex::convert(grid_ba, Bt_fp.ixType().toIntVect());
    WarpX::AllocInitMultiFab(pml_B_fp[0], ba_Br, grid_dm, Br_fp.nComp(), Br_fp.nGrowVect(), "pml_B_fp[0]", 0.0_rt);
    WarpX::AllocInitMultiFab(pml_B_fp[1], ba_Bt, grid_dm, Bt_fp.nComp(), Bt_fp.nGrowVect(), "pml_B_fp[1]", 0.0_rt);

}

void
PML_RZ::ApplyDamping (amrex::MultiFab* Et_fp, amrex::MultiFab* Ez_fp,
                      amrex::MultiFab* Bt_fp, amrex::MultiFab* Bz_fp,
                      amrex::Real dt)
{

    const amrex::Real dr = m_geom->CellSize(0);
    const amrex::Real cdt_over_dr = PhysConst::c*dt/dr;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(*Et_fp, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        amrex::Array4<amrex::Real> const& Et_arr = Et_fp->array(mfi);
        amrex::Array4<amrex::Real> const& Ez_arr = Ez_fp->array(mfi);
        amrex::Array4<amrex::Real> const& Bt_arr = Bt_fp->array(mfi);
        amrex::Array4<amrex::Real> const& Bz_arr = Bz_fp->array(mfi);

        amrex::Array4<amrex::Real> const& pml_Et_arr = pml_E_fp[1]->array(mfi);
        amrex::Array4<amrex::Real> const& pml_Bt_arr = pml_B_fp[1]->array(mfi);

        // Get the tileboxes from Efield and Bfield so that they include the guard cells
        // They are all the same, cell centered
        amrex::Box tilebox = amrex::convert((*Et_fp)[mfi].box(), Et_fp->ixType().toIntVect());

        // Box for the whole simulation domain
        amrex::Box const& domain = m_geom->Domain();
        const int nr_domain = domain.bigEnd(0);

        // Set tilebox to only include the upper radial cells
        const int nr_damp = m_ncell;
        int nr_damp_min;
        if (m_do_pml_in_domain) {
            nr_damp_min = nr_domain - nr_damp;
        } else {
            nr_damp_min = nr_domain;
        }
        tilebox.setSmall(0, nr_damp_min + 1);

        amrex::ParallelFor( tilebox, Et_fp->nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int icomp)
            {
                const amrex::Real rr = static_cast<amrex::Real>(i - nr_damp_min);
                const amrex::Real wr = rr/nr_damp;
                const amrex::Real damp_factor = std::exp( -4._rt * cdt_over_dr * wr*wr );

                // Substract the theta PML fields from the regular theta fields
                Et_arr(i,j,k,icomp) -= pml_Et_arr(i,j,k,icomp);
                Bt_arr(i,j,k,icomp) -= pml_Bt_arr(i,j,k,icomp);
                // Damp the theta PML fields
                pml_Et_arr(i,j,k,icomp) *= damp_factor;
                pml_Bt_arr(i,j,k,icomp) *= damp_factor;
                // Add the theta PML fields back to the regular theta fields
                Et_arr(i,j,k,icomp) += pml_Et_arr(i,j,k,icomp);
                Bt_arr(i,j,k,icomp) += pml_Bt_arr(i,j,k,icomp);
                // Damp the z fields
                Bz_arr(i,j,k,icomp) *= damp_factor;
                Ez_arr(i,j,k,icomp) *= damp_factor;
            });
    }
}

std::array<amrex::MultiFab*,2>
PML_RZ::GetE_fp ()
{
    return {pml_E_fp[0].get(), pml_E_fp[1].get()};
}

std::array<amrex::MultiFab*,2>
PML_RZ::GetB_fp ()
{
    return {pml_B_fp[0].get(), pml_B_fp[1].get()};
}

void
PML_RZ::FillBoundaryE ()
{
    FillBoundaryE(PatchType::fine);
}

void
PML_RZ::FillBoundaryE (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_E_fp[0] && pml_E_fp[0]->nGrowVect().max() > 0)
    {
        const amrex::Periodicity& period = m_geom->periodicity();
        Vector<amrex::MultiFab*> mf{pml_E_fp[0].get(),pml_E_fp[1].get()};
        ablastr::utils::communication::FillBoundary(mf, WarpX::do_single_precision_comms, period);
    }
}

void
PML_RZ::FillBoundaryB ()
{
    FillBoundaryB(PatchType::fine);
}

void
PML_RZ::FillBoundaryB (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_B_fp[0])
    {
        const amrex::Periodicity& period = m_geom->periodicity();
        Vector<amrex::MultiFab*> mf{pml_B_fp[0].get(),pml_B_fp[1].get()};
        ablastr::utils::communication::FillBoundary(mf, WarpX::do_single_precision_comms, period);
    }
}

void
PML_RZ::CheckPoint (const std::string& dir) const
{
    if (pml_E_fp[0])
    {
        VisMF::AsyncWrite(*pml_E_fp[0], dir+"_Er_fp");
        VisMF::AsyncWrite(*pml_E_fp[1], dir+"_Et_fp");
        VisMF::AsyncWrite(*pml_B_fp[0], dir+"_Br_fp");
        VisMF::AsyncWrite(*pml_B_fp[1], dir+"_Bt_fp");
    }
}

void
PML_RZ::Restart (const std::string& dir)
{
    if (pml_E_fp[0])
    {
        VisMF::Read(*pml_E_fp[0], dir+"_Er_fp");
        VisMF::Read(*pml_E_fp[1], dir+"_Et_fp");
        VisMF::Read(*pml_B_fp[0], dir+"_Br_fp");
        VisMF::Read(*pml_B_fp[1], dir+"_Bt_fp");
    }
}

#ifdef WARPX_USE_PSATD
void
PML_RZ::PushPSATD (const int lev)
{
    // Update the fields on the fine and coarse patch
    WarpX& warpx = WarpX::GetInstance();
    SpectralSolverRZ& solver = warpx.get_spectral_solver_fp(lev);
    PushPMLPSATDSinglePatchRZ(lev, solver, pml_E_fp, pml_B_fp);
}

void
PML_RZ::PushPMLPSATDSinglePatchRZ (
    const int lev,
    SpectralSolverRZ& solver,
    std::array<std::unique_ptr<amrex::MultiFab>,2>& pml_E,
    std::array<std::unique_ptr<amrex::MultiFab>,2>& pml_B)
{
    const SpectralFieldIndex& Idx = solver.m_spectral_index;

    // Perform forward Fourier transforms
    solver.ForwardTransform(lev, *pml_E[0], Idx.Er_pml, *pml_E[1], Idx.Et_pml);
    solver.ForwardTransform(lev, *pml_B[0], Idx.Br_pml, *pml_B[1], Idx.Bt_pml);

    // Advance fields in spectral space
    const bool doing_pml = true;
    solver.pushSpectralFields(doing_pml);

    // Perform backward Fourier transforms
    solver.BackwardTransform(lev, *pml_E[0], Idx.Er_pml, *pml_E[1], Idx.Et_pml);
    solver.BackwardTransform(lev, *pml_B[0], Idx.Br_pml, *pml_B[1], Idx.Bt_pml);
}
#endif
