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
     * \param[im] ncomp number of components of mf_dst to update
     */
    virtual void operator()(amrex::MultiFab& mf_dst, int dcomp) const = 0;
    /** Number of component from the input multifab to write to the output
     * multifab */
    int nComp () const { return m_ncomp; };
    /** coarsening ratio for fine-to-coarse averaging */
    amrex::IntVect d_crse_ratio () const {return m_crse_ratio; };
    /** boolean for coarsening */
    bool do_coarsen () const {
        bool coarsen_diag = false;
        for (int idim = 0; idim < AMREX_SPACEDIM; idim ++) {
            if (m_crse_ratio[idim] > 1) coarsen_diag = true;
        }
        return coarsen_diag;
    };
private:
    int m_ncomp;
    amrex::IntVect m_crse_ratio;
};

#endif // WARPX_COMPUTEDIAGFUNCTOR_H_
