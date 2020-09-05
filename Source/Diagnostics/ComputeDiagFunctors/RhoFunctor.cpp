#include "WarpX.H"
#include "RhoFunctor.H"
#include "Utils/CoarsenIO.H"

#ifdef WARPX_DIM_RZ
#include "FieldSolver/SpectralSolver/SpectralFieldData.H"
#endif

RhoFunctor::RhoFunctor ( const int lev, const amrex::IntVect crse_ratio,
                         bool convertRZmodes2cartesian, const int ncomp )
    : ComputeDiagFunctor(ncomp, crse_ratio), m_lev(lev),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{}

void
RhoFunctor::operator() ( amrex::MultiFab& mf_dst, const int dcomp, const int /*i_buffer*/ ) const
{
    auto& warpx = WarpX::GetInstance();
    auto& mypc  = warpx.GetPartContainer();

    // Deposit charge density
    // Call this with local=true since the parallel transfers will be handled
    // by ApplyFilterandSumBoundaryRho
    std::unique_ptr<amrex::MultiFab> rho = mypc.GetChargeDensity( m_lev, true );

    // Call this in case filtering is turned on
    warpx.ApplyFilterandSumBoundaryRho(m_lev, m_lev, *rho, 0, rho->nComp());

#ifdef WARPX_DIM_RZ
    using Idx = SpectralAvgFieldIndex;
    if (WarpX::use_kspace_filter) {
        auto & solver = warpx.get_spectral_solver_fp(m_lev);
        solver.ForwardTransform(*rho, Idx::rho_new);
        solver.ApplyFilter(Idx::rho_new);
        solver.BackwardTransform(*rho, Idx::rho_new);
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
    CoarsenIO::Coarsen( mf_dst, *rho, dcomp, 0, nComp(), mf_dst.nGrow(0), m_crse_ratio );
#endif
}
