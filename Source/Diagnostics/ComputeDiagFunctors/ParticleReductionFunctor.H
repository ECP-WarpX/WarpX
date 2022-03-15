#ifndef WARPX_PARTICLEREDUCTIONFUNCTOR_H_
#define WARPX_PARTICLEREDUCTIONFUNCTOR_H_

#include "ComputeDiagFunctor.H"

#include <AMReX_Parser.H>

#include <AMReX_BaseFwd.H>

#include <memory>
#include <string>

/**
 * \brief Functor to calculate per-cell averages of particle properties.
 */
class
ParticleReductionFunctor final : public ComputeDiagFunctor
{
public:
    /** Constructor.
     * \param[in] mf_src source multifab. Must be nullptr as no source MF is needed
     *            to compute the number of particles per cell, banane.
     * \param[in] lev level of multifab.
     * \param[in] crse_ratio for interpolating field values from simulation MultiFabs
                  to the output diagnostic MultiFab mf_dst.
     * \param[in] fn_str parser string that describes the function to apply to particles
     *            and then average over each cell
     * \param[in] ispec index of the species over which to calculate the average
     * \param[in] ncomp Number of component of mf_src to cell-center in dst multifab.
     */
    ParticleReductionFunctor(const amrex::MultiFab * const mf_src, const int lev,
                       const amrex::IntVect crse_ratio, const std::string fn_str, const int ispec, const int ncomp=1);

    /** \brief Compute the average of the function m_map_fn over each grid cell.
     *
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which cell-centered
     *            data is stored
     */
    virtual void operator()(amrex::MultiFab& mf_dst, const int dcomp, const int /*i_buffer=0*/) const override;
private:
    int const m_lev; /**< level on which mf_src is defined */
    int const m_ispec; /**< index of species to average over */
    /** Parser function to be averaged by the functor. Arguments: x, y, z, ux, uy, uz */
    std::unique_ptr<amrex::Parser> m_map_fn_parser;
    /** Compiled #m_map_fn_parser */
    amrex::ParserExecutor<6> m_map_fn;
};

#endif // WARPX_PARTICLEREDUCTIONFUNCTOR_H_
