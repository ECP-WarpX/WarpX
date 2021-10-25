#include "FlushFormatOpenPMD.H"

#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <map>
#include <memory>
#include <set>
#include <string>

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
      encoding = openPMD::IterationEncoding::variableBased;
    else if ( 0 == openpmd_encoding.compare("g") )
      encoding = openPMD::IterationEncoding::groupBased;
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
        WarpX::GetInstance().RecordWarning("Diagnostics", warnMsg);
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

  // ADIOS2 operator type & parameters
  std::string operator_type;
  pp_diag_name.query("adios2_operator.type", operator_type);
  std::string const prefix = diag_name + ".adios2_operator.parameters";
  ParmParse pp;
  auto entr = pp.getEntries(prefix);

  std::map< std::string, std::string > operator_parameters;
  auto const prefix_len = prefix.size() + 1;
  for (std::string k : entr) {
    std::string v;
    pp.get(k.c_str(), v);
    k.erase(0, prefix_len);
    operator_parameters.insert({k, v});
  }

  auto & warpx = WarpX::GetInstance();
  m_OpenPMDPlotWriter = std::make_unique<WarpXOpenPMDPlot>(
    encoding, openpmd_backend,
    operator_type, operator_parameters,
    warpx.getPMLdirections()
  );

  // Temporarily adding Abort for adios filetype if species is selected for BTD output
  bool species_output = true;
  int write_species = 1;
  std::vector< std::string > output_species_names;
  bool species_specified = pp_diag_name.queryarr("species", output_species_names);
  if (species_specified == true and output_species_names.size() > 0) {
      species_output = true;
  } else {
      // By default species output is computed for all diagnostics, if write_species is not set to 0
      species_output = true;
  }
  // Check user-defined option to turn off species output
  pp_diag_name.query("write_species", write_species);
  if (write_species == 0) species_output = false;
  if (diag_type_str == "BackTransformed" and species_output == true) {
      if (m_OpenPMDPlotWriter->OpenPMDFileType() == "bp") {
          amrex::Abort(" Currently BackTransformed diagnostics type does not support species output for ADIOS backend. Please select h5 as openpmd backend");
      }
  }
}

void
FlushFormatOpenPMD::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<amrex::MultiFab>& mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    const amrex::Vector<ParticleDiag>& particle_diags, int /*nlev*/,
    const std::string prefix, int file_min_digits, bool plot_raw_fields,
    bool plot_raw_fields_guards,
    bool isBTD, int snapshotID, const amrex::Geometry& full_BTD_snapshot,
    bool isLastBTDFlush, const amrex::Vector<int>& totalParticlesFlushedAlready) const
{
    WARPX_PROFILE("FlushFormatOpenPMD::WriteToFile()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !plot_raw_fields && !plot_raw_fields_guards,
        "Cannot plot raw data with OpenPMD output format. Use plotfile instead.");

    // we output at full steps of the coarsest level
    int output_iteration = iteration[0];
    // in backtransformed diagnostics (BTD), we dump into a series of labframe
    // snapshots
    if( isBTD )
        output_iteration = snapshotID;

    // Set step and output directory name.
    m_OpenPMDPlotWriter->SetStep(output_iteration, prefix, file_min_digits, isBTD);

    // fields: only dumped for coarse level
    m_OpenPMDPlotWriter->WriteOpenPMDFieldsAll(
        varnames, mf, geom, output_iteration, time, isBTD, full_BTD_snapshot);

    // particles: all (reside only on locally finest level)
    m_OpenPMDPlotWriter->WriteOpenPMDParticles(particle_diags, isBTD, totalParticlesFlushedAlready);

    // signal that no further updates will be written to this iteration
    m_OpenPMDPlotWriter->CloseStep(isBTD, isLastBTDFlush);
}
