#ifndef WARPX_FLUSHFORMATASCENT_H_
#define WARPX_FLUSHFORMATASCENT_H_

#include "FlushFormat.H"
#include "Diagnostics/ParticleDiag/ParticleDiag.H"

#ifdef AMREX_USE_ASCENT
#   include <ascent.hpp>
#   include <AMReX_Conduit_Blueprint.H>
#endif

/**
 * \brief This class aims at dumping performing in-situ diagnostics with ASCENT.
 * In particular, function WriteToFile takes fields and particles as input arguments,
 * and calls amrex functions to do the in-situ visualization.
 */
class FlushFormatAscent : public FlushFormat
{
public:
    /** Do in-situ visualization for field and particle data */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev, const std::string prefix,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool plot_raw_rho, bool plot_raw_F,
        bool isBTD = false, int snapshotID = -1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false) const override;

    /** \brief Do in-situ visualization for particle data.
     * \param[in] particle_diags Each element of this vector handles output of 1 species.
     * Only compile if AMREX_USE_ASCENT because we need to pass a conduit class
     */
#ifdef AMREX_USE_ASCENT
    void WriteParticles(const amrex::Vector<ParticleDiag>& particle_diags, conduit::Node& a_bp_mesh) const;
#endif

    ~FlushFormatAscent() {}
};

#endif // WARPX_FLUSHFORMATASCENT_H_
