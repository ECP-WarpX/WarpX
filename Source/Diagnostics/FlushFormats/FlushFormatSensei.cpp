#include "FlushFormatSensei.H"
#include "WarpX.H"

#ifdef BL_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif

using namespace amrex;

FlushFormatSensei::FlushFormatSensei (std::string diag_name)
{
#ifdef BL_USE_SENSEI_INSITU

    ParmParse pp(diag_name);

    pp.query("sensei_config", m_insitu_config);
    pp.query("sensei_pin_mesh", m_insitu_pin_mesh);

    m_insitu_bridge = new amrex::AmrMeshInSituBridge;
    m_insitu_bridge->setEnabled(true);
    m_insitu_bridge->setConfig(m_insitu_config);
    m_insitu_bridge->setPinMesh(m_insitu_pin_mesh);
    if (m_insitu_bridge->initialize())
    {
        amrex::ErrorStream()
            << "WarpX::InitData : Failed to initialize the in situ bridge."
            << std::endl;

        amrex::Abort();
    }
    m_insitu_bridge->setFrequency(1);

#endif
}

FlushFormatSensei::~FlushFormatSensei ()
{
#ifdef BL_USE_SENSEI_INSITU
    delete m_insitu_bridge;
#endif
}

void
FlushFormatSensei::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<const amrex::MultiFab*> mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
    const std::string prefix, bool plot_raw_fields,
    bool plot_raw_fields_guards, bool plot_raw_rho, bool plot_raw_F) const
{
#ifdef BL_USE_SENSEI_INSITU
    if (insitu_bridge->update(iteration, time,
        dynamic_cast<amrex::AmrMesh*>(const_cast<WarpX*>(this)),
        {&mf}, {varnames}))
    {
        amrex::ErrorStream()
            << "WarpXIO::UpdateInSitu : Failed to update the in situ bridge."
            << std::endl;

        amrex::Abort();
    }
#endif
}
