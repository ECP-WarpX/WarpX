#include "CellCenterFunctor.H"
#include "Utils/CoarsenIO.H"

CellCenterFunctor::CellCenterFunctor(amrex::MultiFab const * mf_src, int lev,
                                     amrex::IntVect crse_ratio,
                                     bool convertRZmodes2cartesian, int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{}

void
CellCenterFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp, int i_buffer) const
{
    amrex::Print() << " cc operator : " << dcomp << " ibuffer " << i_buffer << "\n";
#ifdef WARPX_DIM_RZ
    if (m_convertRZmodes2cartesian) {
        // In cylindrical geometry, sum real part of all modes of m_mf_src in
        // temporary multifab mf_dst_stag, and cell-center it to mf_dst.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into 1 single component");
        auto& warpx = WarpX::GetInstance();
        amrex::MultiFab mf_dst_stag(m_mf_src->boxArray(), warpx.DistributionMap(m_lev), 1, m_mf_src->nGrowVect());
        // Mode 0
        amrex::MultiFab::Copy(mf_dst_stag, *m_mf_src, 0, 0, 1, m_mf_src->nGrowVect());
        for (int ic=1 ; ic < m_mf_src->nComp() ; ic += 2) {
            // All modes > 0
            amrex::MultiFab::Add(mf_dst_stag, *m_mf_src, ic, 0, 1, m_mf_src->nGrowVect());
        }
        CoarsenIO::Coarsen( mf_dst, mf_dst_stag, dcomp, 0, nComp(), 0,  m_crse_ratio);
    } else {
        CoarsenIO::Coarsen( mf_dst, *m_mf_src, dcomp, 0, nComp(), 0, m_crse_ratio);
    }
#else
    // In cartesian geometry, coarsen and interpolate from simulation MultiFab, m_mf_src,
    // to output diagnostic MultiFab, mf_dst.
    amrex::Print() << " before coarsen " << dcomp << " nComp() " << nComp() << "\n";
    CoarsenIO::Coarsen( mf_dst, *m_mf_src, dcomp, 0, nComp(), mf_dst.nGrow(0), m_crse_ratio);
    amrex::Print() << " after coarsen " << dcomp << " nComp() " << nComp() << "\n";
#endif
}
