#include "DivEFunctor.H"
#include "Utils/CoarsenIO.H"

using namespace amrex;

DivEFunctor::DivEFunctor(const std::array<const amrex::MultiFab* const, 3> arr_mf_src, const int lev,
                         const amrex::IntVect crse_ratio,
                         bool convertRZmodes2cartesian, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_arr_mf_src(arr_mf_src), m_lev(lev),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{}

void
DivEFunctor::operator()(amrex::MultiFab& mf_dst, const int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // For staggered and nodal calculations, divE is computed on the nodes.
    // The temporary divE MultiFab is generated to comply with the location of divE.
    const BoxArray& ba = amrex::convert(warpx.boxArray(m_lev),IntVect::TheUnitVector());
    MultiFab divE(ba, warpx.DistributionMap(m_lev), 2*warpx.n_rz_azimuthal_modes-1, ng );
    warpx.ComputeDivE(divE, m_lev);

#ifdef WARPX_DIM_RZ
    if (m_convertRZmodes2cartesian) {
        // In cylindrical geometry, sum real part of all modes of divE in
        // temporary multifab mf_dst_stag, and cell-center it to mf_dst.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nComp()==1,
            "The RZ averaging over modes must write into 1 single component");
        MultiFab mf_dst_stag(divE.boxArray(), warpx.DistributionMap(m_lev), 1, divE.nGrowVect());
        // Mode 0
        MultiFab::Copy(mf_dst_stag, divE, 0, 0, 1, divE.nGrowVect());
        for (int ic=1 ; ic < divE.nComp() ; ic += 2) {
            // Real part of all modes > 0
            MultiFab::Add(mf_dst_stag, divE, ic, 0, 1, divE.nGrowVect());
        }
        CoarsenIO::Coarsen( mf_dst, mf_dst_stag, dcomp, 0, nComp(), 0,  m_crse_ratio);
    } else {
        CoarsenIO::Coarsen( mf_dst, divE, dcomp, 0, nComp(), 0, m_crse_ratio);
    }
#else
    // In cartesian geometry, coarsen and interpolate from simulation MultiFab, divE,
    // to output diagnostic MultiFab, mf_dst.
    CoarsenIO::Coarsen( mf_dst, divE, dcomp, 0, nComp(), 0, m_crse_ratio);
#endif
}
