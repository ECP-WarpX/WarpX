#ifndef WARPX_FLUSHFORMAT_H_
#define WARPX_FLUSHFORMAT_H_

#include "FlushFormat_fwd.H"

#include "Diagnostics/ParticleDiag/ParticleDiag.H"
#include "Particles/MultiParticleContainer.H"
#include "WarpX.H"

class FlushFormat
{
public:
    /** Flush fields and particles to file */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
        const std::string prefix, const int file_min_digits,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        const bool use_pinned_pc = false,
        bool isBTD = false, int snapshotID = -1,
        int bufferID = 1, int numBuffers = 1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false,
        const amrex::Vector<int>& totalParticlesFlushedAlready = amrex::Vector<int>() ) const = 0;

     virtual ~FlushFormat() {}
};

#endif // WARPX_FLUSHFORMAT_H_
