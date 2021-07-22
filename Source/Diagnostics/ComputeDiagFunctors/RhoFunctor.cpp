#include "RhoFunctor.H"

#include "Diagnostics/ComputeDiagFunctors/ComputeDiagFunctor.H"
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
    #include "FieldSolver/SpectralSolver/SpectralFieldData.H"
    #include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
    #include "Utils/WarpXAlgorithmSelection.H"
#endif
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/CoarsenIO.H"
#include "WarpX.H"

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
    auto& warpx = WarpX::GetInstance();
    std::unique_ptr<amrex::MultiFab> rho;

    // Deposit charge density
    // Call this with local=true since the parallel transfers will be handled
    // by ApplyFilterandSumBoundaryRho

    // Dump total rho
    if (m_species_index == -1) {
        auto& mypc = warpx.GetPartContainer();
        rho = mypc.GetChargeDensity(m_lev, true);
    }
    // Dump rho per species
    else {
        auto& mypc = warpx.GetPartContainer().GetParticleContainer(m_species_index);
        rho = mypc.GetChargeDensity(m_lev, true);
    }

    // Handle the parallel transfers of guard cells and
    // apply the filtering if requested.
    warpx.ApplyFilterandSumBoundaryRho(m_lev, m_lev, *rho, 0, rho->nComp());

#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
    // Apply k-space filtering when using the PSATD solver
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD)
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
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into one single component");
        amrex::MultiFab mf_dst_stag( rho->boxArray(), warpx.DistributionMap(m_lev), 1, rho->nGrowVect() );
        // Mode 0
        amrex::MultiFab::Copy( mf_dst_stag, *rho, 0, 0, 1, rho->nGrowVect() );
        for (int ic=1 ; ic < rho->nComp() ; ic += 2) {
            // Real part of all modes > 0
            amrex::MultiFab::Add( mf_dst_stag, *rho, ic, 0, 1, rho->nGrowVect() );
        }
        CoarsenIO::Coarsen( mf_dst, mf_dst_stag, dcomp, 0, nComp(), 0, m_crse_ratio );
    } else {
        CoarsenIO::Coarsen( mf_dst, *rho, dcomp, 0, nComp(), 0, m_crse_ratio );
    }
#else
    // In Cartesian geometry, coarsen and interpolate from temporary MultiFab rho
    // to output diagnostic MultiFab mf_dst
    CoarsenIO::Coarsen( mf_dst, *rho, dcomp, 0, nComp(), mf_dst.nGrowVect(), m_crse_ratio );
    amrex::ignore_unused(m_convertRZmodes2cartesian);
#endif
}
