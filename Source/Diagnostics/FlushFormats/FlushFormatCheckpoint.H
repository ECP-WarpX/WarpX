#ifndef WARPX_FLUSHFORMATCHECKPOINT_H_
#define WARPX_FLUSHFORMATCHECKPOINT_H_

#include "FlushFormatPlotfile.H"

class FlushFormatCheckpoint final : public FlushFormatPlotfile
{
    /** Flush fields and particles to plotfile */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev, const std::string prefix,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool plot_raw_rho, bool plot_raw_F) const override final;

    void CheckpointParticles(const std::string& dir,
                             const amrex::Vector<ParticleDiag>& particle_diags) const;
};

#endif // WARPX_FLUSHFORMATCHECKPOINT_H_
