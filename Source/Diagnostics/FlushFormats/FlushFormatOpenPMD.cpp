#include "FlushFormatOpenPMD.H"
#include "WarpX.H"
#include "Utils/Interpolate.H"

#include <AMReX_buildInfo.H>

using namespace amrex;

FlushFormatOpenPMD::FlushFormatOpenPMD (const std::string& diag_name)
{
    ParmParse pp(diag_name);
    // Which backend to use (ADIOS, ADIOS2 or HDF5). Default depends on what is available
    std::string openpmd_backend {"default"};
    // one file per timestep (or one file for all steps)
    bool openpmd_tspf = true;
    pp.query("openpmd_backend", openpmd_backend);
    pp.query("openpmd_tspf", openpmd_tspf);
    auto & warpx = WarpX::GetInstance();
    m_OpenPMDPlotWriter = new WarpXOpenPMDPlot(
        openpmd_tspf, openpmd_backend, warpx.getPMLdirections()
        );
}

void
FlushFormatOpenPMD::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<amrex::MultiFab>& mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    const amrex::Vector<ParticleDiag>& particle_diags, int /*nlev*/,
    const std::string prefix, bool plot_raw_fields,
    bool plot_raw_fields_guards, bool plot_raw_rho, bool plot_raw_F,
    bool isBTD, int snapshotID, const amrex::Geometry& full_BTD_snapshot,
    bool isLastBTDFlush) const
{
    WARPX_PROFILE("FlushFormatOpenPMD::WriteToFile()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !plot_raw_fields && !plot_raw_fields_guards && !plot_raw_rho && !plot_raw_F,
        "Cannot plot raw data with OpenPMD output format. Use plotfile instead.");

    // we output at full steps of the coarsest level
    int output_iteration = iteration[0];
    // in backtransformed diagnostics (BTD), we dump into a series of labframe
    // snapshots
    if( isBTD )
        output_iteration = snapshotID;

    // Set step and output directory name.
    m_OpenPMDPlotWriter->SetStep(output_iteration, prefix, isBTD);

    // fields: only dumped for coarse level
    m_OpenPMDPlotWriter->WriteOpenPMDFields(
        varnames, mf[0], geom[0], output_iteration, time, isBTD, full_BTD_snapshot);

    // particles: all (reside only on locally finest level)
    m_OpenPMDPlotWriter->WriteOpenPMDParticles(particle_diags);

    // signal that no further updates will be written to this iteration
    m_OpenPMDPlotWriter->CloseStep(isBTD, isLastBTDFlush);
}

FlushFormatOpenPMD::~FlushFormatOpenPMD (){
    delete m_OpenPMDPlotWriter;
}
