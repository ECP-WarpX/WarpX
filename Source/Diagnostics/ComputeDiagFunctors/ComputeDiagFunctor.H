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
    ComputeDiagFunctor( int ncomp, amrex::IntVect diag_crse_ratio) :
                        m_ncomp(ncomp), m_crse_ratio(diag_crse_ratio) {};
    /** Compute a field and store the result in mf_dst
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which the result is written
     * \param[in] ncomp number of components of mf_dst to update
     * \param[in] m_crse_ratio Source MultiFab is coarsened by this ratio to dst
     */
    virtual void operator() (amrex::MultiFab& mf_dst, int dcomp) const = 0;
    /** Number of component from the input multifab to write to the output
     * multifab */
    int nComp () const { return m_ncomp; };
private:
    int m_ncomp;
protected:
    amrex::IntVect m_crse_ratio;
};

#endif // WARPX_COMPUTEDIAGFUNCTOR_H_
