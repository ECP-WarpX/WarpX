#ifndef WARPX_FLUSHFORMATSENSEI_H_
#define WARPX_FLUSHFORMATSENSEI_H_

#include "FlushFormat.H"

#include <AMReX_AmrMesh.H>
#if defined(BL_USE_SENSEI_INSITU)
#  include <AMReX_AmrMeshInSituBridge.H>
#else
namespace amrex {
using AmrMeshInSituBridge = void;
}
#endif

/**
 * \brief This class aims at dumping performing in-situ diagnostics with
 * SENSEI.  In particular, function WriteToFile takes fields and particles as
 * input arguments, and calls amrex functions to do the in-situ visualization.
 *
 * To use the SENSEI flush format one must compile WarpX with the GNUMakefile
 * setting:  USE_SENSEI_INSITU=TRUE
 *
 * In addition to the standard controls, the following inputs file fields are
 * supported:
 *
 *      sensei_config - the path to a SENSEI XML configuration (required)
 *      sensei_pin_mesh - integer 0 or 1 forcing the moving mesh to be fixed
 *                        at 0,0,0 (optional)
 *
 * Fri 29 May 2020 11:19:38 AM PDT : Tested with SENSEI version 3.2.0
 */
class FlushFormatSensei : public FlushFormat
{
public:
    FlushFormatSensei();
    ~FlushFormatSensei();

    FlushFormatSensei(const FlushFormatSensei &) = delete;
    void operator=(const FlushFormatSensei &) = delete;

    /** \breif Construct.
     * \parm[in] amr_mesh an AmrMesh instance
     * \parm[in] diag_name ParmParse scope string.
     */
    FlushFormatSensei (amrex::AmrMesh *amr_mesh, std::string diag_name);

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
     */
    void WriteParticles(const amrex::Vector<ParticleDiag>& particle_diags) const;

private:
    std::string m_insitu_config;
    int m_insitu_pin_mesh = 0;
    amrex::AmrMeshInSituBridge * m_insitu_bridge = nullptr;
    amrex::AmrMesh * m_amr_mesh = nullptr;
};

#endif // WARPX_FLUSHFORMATSENSEI_H_
