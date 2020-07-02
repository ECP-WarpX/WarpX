#include "BackTransformFunctor.H"

BackTransformFunctor::BackTransformFunctor(amrex::MultiFab const * mf_src, int lev,
                                           const int ncomp, const amrex::IntVect crse_ratio)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev)
{
    // Get the right slice of each field in the CC MultiFab, transform it and store it in the output.
}
