#ifndef WARPX_COMPUTEDIAGFUNCTOR_H_
#define WARPX_COMPUTEDIAGFUNCTOR_H_

#include <AMReX_MultiFab.H>

/**
 * \brief Functor to compute a diagnostic and store the result in existing
 * MultiFab
 */
class
ComputeDiagFunctor
{
public:
    /** Compute a field and store the result in mf_dst
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which the result is written
     * \param[im] ncomp number of components of mf_dst to update
     */
    virtual void operator() (amrex::MultiFab& mf_dst, int dcomp, int ncomp=1) const = 0;
};

#endif // WARPX_COMPUTEDIAGFUNCTOR_H_
