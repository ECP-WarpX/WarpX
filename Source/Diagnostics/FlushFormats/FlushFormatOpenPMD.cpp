#include "FlushFormatOpenPMD.H"
#include "WarpX.H"

#include <memory>

using namespace amrex;


FlushFormatOpenPMD::FlushFormatOpenPMD (const std::string& diag_name)
{
    ParmParse pp_diag_name(diag_name);
    // Which backend to use (ADIOS, ADIOS2 or HDF5). Default depends on what is available
    std::string openpmd_backend {"default"};
    // one file per timestep (or one file for all steps)
    std::string  openpmd_encoding {"f"};
    pp_diag_name.query("openpmd_backend", openpmd_backend);
    bool encodingDefined = pp_diag_name.query("openpmd_encoding", openpmd_encoding);

    openPMD::IterationEncoding encoding = openPMD::IterationEncoding::groupBased;
    if ( 0 == openpmd_encoding.compare("v") )
#if OPENPMDAPI_VERSION_GE(0, 14, 0)
      encoding = openPMD::IterationEncoding::variableBased;
#else
      encoding = openPMD::IterationEncoding::groupBased;
#endif
    else if ( 0 == openpmd_encoding.compare("f") )
      encoding = openPMD::IterationEncoding::fileBased;

    std::string diag_type_str;
    pp_diag_name.get("diag_type", diag_type_str);
    if (diag_type_str == "BackTransformed")
    {
      if ( ( openPMD::IterationEncoding::fileBased != encoding ) &&
           ( openPMD::IterationEncoding::groupBased != encoding ) )
      {
        std::string warnMsg = diag_name+" Unable to support BTD with streaming. Using GroupBased ";
        amrex::Warning(warnMsg);
        encoding = openPMD::IterationEncoding::groupBased;
      }
    }

  //
  // if no encoding is defined, then check to see if tspf is defined.
  // (backward compatibility)
  //
  if ( !encodingDefined )
    {
      bool openpmd_tspf = false;
      bool tspfDefined = pp_diag_name.query("openpmd_tspf", openpmd_tspf);
      if ( tspfDefined && openpmd_tspf )
          encoding = openPMD::IterationEncoding::fileBased;
    }

  auto & warpx = WarpX::GetInstance();
  m_OpenPMDPlotWriter = std::make_unique<WarpXOpenPMDPlot>(
                               encoding, openpmd_backend, warpx.getPMLdirections()
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
    m_OpenPMDPlotWriter->WriteOpenPMDFieldsAll(
        varnames, mf, geom, output_iteration, time, isBTD, full_BTD_snapshot);

    // particles: all (reside only on locally finest level)
    m_OpenPMDPlotWriter->WriteOpenPMDParticles(particle_diags);

    // signal that no further updates will be written to this iteration
    m_OpenPMDPlotWriter->CloseStep(isBTD, isLastBTDFlush);
}
