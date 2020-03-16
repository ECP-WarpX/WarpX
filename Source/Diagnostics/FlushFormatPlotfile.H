#ifndef WARPX_FLUSHFORMATPLOTFILE_H_
#define WARPX_FLUSHFORMATPLOTFILE_H_

#include "FlushFormat.H"

class FlushFormatPlotfile : public FlushFormat
{
public:
    virtual void WriteToFile (const amrex::Vector<std::string> varnames,
                              const amrex::Vector<const amrex::MultiFab*> mf,
                              amrex::Vector<amrex::Geometry>& geom,
                              const amrex::Vector<int> iteration, const double time,
                              MultiParticleContainer& mpc, int nlev) const override final;

    void WriteJobInfo(const std::string& dir) const;
    
    void WriteWarpXHeader(const std::string& dir) const;

};

#endif // WARPX_FLUSHFORMATPLOTFILE_H_
