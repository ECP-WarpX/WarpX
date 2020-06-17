#include "WarpX.H"
#include "RhoFunctor.H"
#include "Utils/CoarsenIO.H"

RhoFunctor::RhoFunctor ( const amrex::MultiFab* const mf_src, const int lev,
                         const amrex::IntVect crse_ratio, bool convertRZmodes2cartesian,
                         const int ncomp )
    : ComputeDiagFunctor(ncomp, crse_ratio), m_lev(lev),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{
    // mf_src should not be used, let's make sure it is null
    AMREX_ALWAYS_ASSERT(mf_src == nullptr);
    // Write only in one output component
    AMREX_ALWAYS_ASSERT(ncomp == 1);
}

void
RhoFunctor::operator() ( amrex::MultiFab& mf_dst, const int dcomp ) const
{
    auto& warpx = WarpX::GetInstance();
    auto& mypc  = warpx.GetPartContainer();

    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab the guard-cell data is not needed, especially considering
    // the operations performed in CoarsenIO::Coarsen
    constexpr int ng = 1;

    // For staggered and nodal calculations, rho is computed on the nodes:
    // thus, the temporary MultiFab rho is generated with nodal index type
    const amrex::BoxArray& ba = amrex::convert( warpx.boxArray(m_lev), amrex::IntVect::TheUnitVector() );
    std::unique_ptr<amrex::MultiFab> rho_ptr;
    rho_ptr = std::unique_ptr<amrex::MultiFab>
              (new amrex::MultiFab( ba, warpx.DistributionMap(m_lev), 2*warpx.n_rz_azimuthal_modes-1, ng) );
    rho_ptr = mypc.GetChargeDensity( m_lev );
    amrex::MultiFab& rho = *rho_ptr;

#ifdef WARPX_DIM_RZ
    if (m_convertRZmodes2cartesian) {
        // In cylindrical geometry, sum real part of all modes of rho in
        // temporary multifab mf_dst_stag, and cell-center it to mf_dst
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into 1 single component");
        amrex::MultiFab mf_dst_stag( rho.boxArray(), warpx.DistributionMap(m_lev), 1, rho.nGrowVect() );
        // Mode 0
        amrex::MultiFab::Copy( mf_dst_stag, rho, 0, 0, 1, rho.nGrowVect() );
        for (int ic=1 ; ic < rho.nComp() ; ic += 2) {
            // Real part of all modes > 0
            amrex::MultiFab::Add(mf_dst_stag, rho, ic, 0, 1, rho.nGrowVect() );
        }
        CoarsenIO::Coarsen( mf_dst, mf_dst_stag, dcomp, 0, nComp(), 0, m_crse_ratio );
    } else {
        CoarsenIO::Coarsen( mf_dst, rho, dcomp, 0, nComp(), 0, m_crse_ratio );
    }
#else
    // In cartesian geometry, coarsen and interpolate from simulation MultiFab rho
    // to output diagnostic MultiFab mf_dst
    CoarsenIO::Coarsen( mf_dst, rho, dcomp, 0, nComp(), 0, m_crse_ratio );
#endif
}
