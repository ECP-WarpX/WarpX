#ifndef WARPX_FLUSHFORMATCHECKPOINT_H_
#define WARPX_FLUSHFORMATCHECKPOINT_H_

#include "FlushFormatPlotfile.H"

#include "Diagnostics/ParticleDiag/ParticleDiag_fwd.H"

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <string>

class FlushFormatCheckpoint final : public FlushFormatPlotfile
{
    /** Flush fields and particles to plotfile */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
        const std::string prefix, int file_min_digits,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool isBTD = false, int snapshotID = -1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false,
        const amrex::Vector<int>& totalParticlesFlushedAlready = amrex::Vector<int>() ) const override final;

    void CheckpointParticles (const std::string& dir,
                              const amrex::Vector<ParticleDiag>& particle_diags) const;

    void WriteDMaps (const std::string& dir, int nlev) const;
};

#endif // WARPX_FLUSHFORMATCHECKPOINT_H_
