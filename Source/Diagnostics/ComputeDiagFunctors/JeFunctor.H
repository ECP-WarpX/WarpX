#ifndef WARPX_JEFUNCTOR_H
#define WARPX_JEFUNCTOR_H

#include "ComputeDiagFunctor.H"

#include <AMReX_BaseFwd.H>

/**
 * \brief Functor to cell-center MF for electron current density from Kinetic-Fluid
 * Hybrid model and store result in mf_out.
 */
class
JeFunctor final: public ComputeDiagFunctor
{
public:
    /**
     *
     * \param[in] dir director o vector field to operate on
     * \param[in] lev level of multifab. Used for averaging in rz.
     * \param[in] crse_ratio for interpoalting field values from the simulation Multifab, src_mf,
     *            to the output diagnostic MultiFab, mf_dst.
     * \param[in] ncomp Number of component of mf_src to cell-center in dst multifab.
    */
    JeFunctor (int dir, int lev,
              amrex::IntVect crse_ratio,
              bool deposit_current, int ncomp=1);
    /** \brief Cell-center m_mf_src and write the result in mf_dst.
     *
     * Currently, the Kinetic-Fluid Hybrid model only supports the m=0 mode, so only
     * that component is written.
     *
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which cell-centered
     *            data is stored
    */
    virtual void operator()(amrex::MultiFab& mf_dst, int dcomp, int /*i_buffer=0*/) const override;
private:
    /** direction of the electron current density to save */
    const int m_dir;
    /** level on which mf_src is defined */
    int m_lev;
};

#endif // WARPX_JEFUNCTOR_H_
