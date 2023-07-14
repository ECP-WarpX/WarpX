#include "RhoFunctor.H"

#include "Diagnostics/ComputeDiagFunctors/ComputeDiagFunctor.H"
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
    #include "FieldSolver/SpectralSolver/SpectralFieldData.H"
    #include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
    #include "Utils/WarpXAlgorithmSelection.H"
#endif
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <memory>

RhoFunctor::RhoFunctor (const int lev,
                        const amrex::IntVect crse_ratio,
                        const int species_index,
                        bool convertRZmodes2cartesian,
                        const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio),
      m_lev(lev),
      m_species_index(species_index),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{}

void
RhoFunctor::operator() ( amrex::MultiFab& mf_dst, const int dcomp, const int /*i_buffer*/ ) const
{
    // Prepare temporary MultiFabs, in order to sum across levels
    auto& warpx = WarpX::GetInstance();
    const int nlevs_max = warpx.maxLevel() + 1;
    amrex::Vector< std::unique_ptr<amrex::MultiFab> > rho_fp;
    amrex::Vector< std::unique_ptr<amrex::MultiFab> > rho_cp;
    amrex::Vector< std::unique_ptr<amrex::MultiFab> > rho_buf;
    rho_cp.resize(nlevs_max);
    rho_fp.resize(nlevs_max);
    rho_buf.resize(nlevs_max);
    const int ng_rho = warpx.get_ng_depos_rho().max();
    for (int lev=m_lev; lev<=warpx.maxLevel(); lev++) {
        // Allocate rho_fp
        const amrex::DistributionMapping& dm = warpx.DistributionMap(lev);
        amrex::BoxArray nba = warpx.boxArray(lev);
        nba.surroundingNodes();
        rho_fp[lev] = std::make_unique<amrex::MultiFab>(nba,dm,WarpX::ncomps,ng_rho);
        rho_fp[lev]->setVal(0.0);
        // Allocate rho_cp
        if (lev>m_lev) {
            nba.coarsen(warpx.refRatio(lev-1));
            rho_cp[lev] = std::make_unique<amrex::MultiFab>(nba,dm,WarpX::ncomps,ng_rho);
        }
    }

    // Deposit charge density
    // Call this with local=true since the parallel transfers will be handled
    // by ApplyFilterandSumBoundaryRho

    // Dump total rho
    if (m_species_index == -1) {
        auto& mypc = warpx.GetPartContainer();
        for (int lev=m_lev; lev<=warpx.finestLevel(); lev++) {
            rho_fp[lev] = mypc.GetChargeDensity(lev, true);
        }
    }
    // Dump rho per species
    else {
        auto& mypc = warpx.GetPartContainer().GetParticleContainer(m_species_index);
        for (int lev=m_lev; lev<=warpx.finestLevel(); lev++) {
            rho_fp[lev] = mypc.GetChargeDensity(lev, true);
        }
    }

    // Handle parallel transfers of guard cells,
    // summation between levels, and filtering
    warpx.SyncRho(rho_fp, rho_cp, rho_buf, m_lev);

    // Select the correct level
    std::unique_ptr<amrex::MultiFab>& rho = rho_fp[m_lev];

#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
    // Apply k-space filtering when using the PSATD solver
    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD)
    {
        if (WarpX::use_kspace_filter) {
            auto & solver = warpx.get_spectral_solver_fp(m_lev);
            const SpectralFieldIndex& Idx = solver.m_spectral_index;
            solver.ForwardTransform(m_lev, *rho, Idx.rho_new);
            solver.ApplyFilter(m_lev, Idx.rho_new);
            solver.BackwardTransform(m_lev, *rho, Idx.rho_new);
        }
    }
#endif


#ifdef WARPX_DIM_RZ
    if (m_convertRZmodes2cartesian) {
        // In cylindrical geometry, sum real part of all modes of rho in
        // temporary MultiFab mf_dst_stag, and cell-center it to mf_dst
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into one single component");
        amrex::MultiFab mf_dst_stag( rho->boxArray(), warpx.DistributionMap(m_lev), 1, rho->nGrowVect() );
        // Mode 0
        amrex::MultiFab::Copy( mf_dst_stag, *rho, 0, 0, 1, rho->nGrowVect() );
        for (int ic=1 ; ic < rho->nComp() ; ic += 2) {
            // Real part of all modes > 0
            amrex::MultiFab::Add( mf_dst_stag, *rho, ic, 0, 1, rho->nGrowVect() );
        }
        ablastr::coarsen::sample::Coarsen( mf_dst, mf_dst_stag, dcomp, 0, nComp(), 0, m_crse_ratio );
    } else {
        ablastr::coarsen::sample::Coarsen( mf_dst, *rho, dcomp, 0, nComp(), 0, m_crse_ratio );
    }
#else
    // In Cartesian geometry, coarsen and interpolate from temporary MultiFab rho
    // to output diagnostic MultiFab mf_dst
    ablastr::coarsen::sample::Coarsen(mf_dst, *rho, dcomp, 0, nComp(), mf_dst.nGrowVect(), m_crse_ratio );
    amrex::ignore_unused(m_convertRZmodes2cartesian);
#endif
}
