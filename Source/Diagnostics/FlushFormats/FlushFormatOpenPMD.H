#ifndef WARPX_FLUSHFORMATOPENPMD_H_
#define WARPX_FLUSHFORMATOPENPMD_H_

#include "FlushFormat.H"
#include "Diagnostics/WarpXOpenPMD.H"

class FlushFormatOpenPMD : public FlushFormat
{
public:

    FlushFormatOpenPMD (const std::string& diag_name);
    
    /** Flush fields and particles to plotfile */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<const amrex::MultiFab*> mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        MultiParticleContainer& mpc, int nlev, const std::string prefix,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool plot_rho, bool plot_F) const override final;

    ~FlushFormatOpenPMD ();

private:
    std::string openpmd_backend {"default"};
    bool openpmd_tspf = true; //!< one file per timestep (or one file for all steps)
    WarpXOpenPMDPlot* m_OpenPMDPlotWriter = nullptr;
};

#endif // WARPX_FLUSHFORMATOPENPMD_H_
