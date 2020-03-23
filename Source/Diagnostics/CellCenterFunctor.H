#ifndef WARPX_CELLCENTERFUNCTOR_H_
#define WARPX_CELLCENTERFUNCTOR_H_

#include "ComputeDiagFunctor.H"

/**
 * \brief Functor to cell-center MF and store result in mf_out.
 */
class
CellCenterFunctor : public ComputeDiagFunctor
{
public:
    CellCenterFunctor(amrex::MultiFab* mf_src, int lev, bool mode_avg=true, int ncomp=1);
    /** Cell-center m_mf_src and write the result in mf_dst
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which the result is written
     * \param[im] ncomp number of components of mf_dst to update
     */
    virtual void operator()(amrex::MultiFab& mf_dst, int dcomp) const override final;
private:
    amrex::MultiFab* m_mf_src;
    int m_lev;
    bool m_mode_avg;
};

#endif // WARPX_CELLCENTERFUNCTOR_H_
