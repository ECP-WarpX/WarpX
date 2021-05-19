#ifndef WARPX_FLUSHFORMATOPENPMD_H_
#define WARPX_FLUSHFORMATOPENPMD_H_

#include "FlushFormat.H"
#include "Diagnostics/WarpXOpenPMD.H"

#include <memory>


/**
 * \brief This class aims at dumping diags data to disk using the openPMD standard.
 * In particular, function WriteToFile takes fields and particles as input arguments,
 * and writes data to file.
 */
class FlushFormatOpenPMD final : public FlushFormat
{
public:

    /** Constructor takes name of diagnostics to set the output directory */
    explicit FlushFormatOpenPMD (const std::string& diag_name);

    /** Flush fields and particles to plotfile */
    void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int /*nlev*/, const std::string prefix,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool plot_raw_rho, bool plot_raw_F,
        bool isBTD = false, int snapshotID = -1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false) const override;

    ~FlushFormatOpenPMD () override = default;

private:
    /** This is responsible for dumping to file */
    std::unique_ptr< WarpXOpenPMDPlot > m_OpenPMDPlotWriter;
};

#endif // WARPX_FLUSHFORMATOPENPMD_H_
