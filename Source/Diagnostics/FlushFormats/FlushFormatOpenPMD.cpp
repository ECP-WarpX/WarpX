#include "FlushFormatOpenPMD.H"
#include "WarpX.H"
#include "Utils/Interpolate.H"

#include <AMReX_buildInfo.H>

using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

FlushFormatOpenPMD::FlushFormatOpenPMD (const std::string& diag_name)
{
    ParmParse pp(diag_name);
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
    const amrex::Vector<const amrex::MultiFab*> mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    MultiParticleContainer& mpc, int nlev,
    const std::string prefix, bool plot_raw_fields,
    bool plot_raw_fields_guards, bool plot_rho, bool plot_F) const
{
    WARPX_PROFILE("FlushFormatOpenPMD::WriteToFile()");

    //auto & warpx = WarpX::GetInstance();
    const auto step = iteration[0];

    m_OpenPMDPlotWriter->SetStep(step, prefix);

    /*
    Vector<std::string> varnames; // Name of the written fields
    Vector<MultiFab> mf_avg; // contains the averaged, cell-centered fields
    Vector<const MultiFab*> output_mf; // will point to the data to be written
    Vector<Geometry> output_geom;

    prepareFields(step, varnames, mf_avg, output_mf, output_geom);
    */
    // fields: only dumped for coarse level
    m_OpenPMDPlotWriter->WriteOpenPMDFields(
        varnames, *mf[0], geom[0], step, iteration[0]);
    // particles: all (reside only on locally finest level)
    m_OpenPMDPlotWriter->WriteOpenPMDParticles(mpc);
}

FlushFormatOpenPMD::~FlushFormatOpenPMD (){
    delete m_OpenPMDPlotWriter;
}
