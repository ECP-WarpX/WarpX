#include "FlushFormatSensei.H"
#include "WarpX.H"

#ifdef BL_USE_SENSEI_INSITU
# include <AMReX_AmrMeshInSituBridge.H>
#endif

FlushFormatSensei::FlushFormatSensei () :
    m_insitu_config(), m_insitu_pin_mesh(0), m_insitu_bridge(nullptr),
    m_amr_mesh(nullptr)
{}

FlushFormatSensei::FlushFormatSensei (amrex::AmrMesh *amr_mesh,
    std::string diag_name) :
    m_insitu_config(), m_insitu_pin_mesh(0), m_insitu_bridge(nullptr),
    m_amr_mesh(amr_mesh)
{
#ifndef BL_USE_SENSEI_INSITU
    (void)amr_mesh;
    (void)diag_name;
#else
    amrex::ParmParse pp(diag_name);

    pp.query("sensei_config", m_insitu_config);
    pp.query("sensei_pin_mesh", m_insitu_pin_mesh);

    m_insitu_bridge = new amrex::AmrMeshInSituBridge;
    m_insitu_bridge->setEnabled(true);
    m_insitu_bridge->setConfig(m_insitu_config);
    m_insitu_bridge->setPinMesh(m_insitu_pin_mesh);
    if (!m_amr_mesh || m_insitu_bridge->initialize())
    {
        amrex::ErrorStream() << "FlushFormtSensei::FlushFormatSensei : "
            "Failed to initialize the in situ bridge." << std::endl;

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
    const amrex::Vector<amrex::MultiFab>& mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
    const std::string prefix, bool plot_raw_fields,
    bool plot_raw_fields_guards, bool plot_raw_rho, bool plot_raw_F) const
{
#ifndef BL_USE_SENSEI_INSITU
    (void)varnames;
    (void)mf;
    (void)geom;
    (void)iteration;
    (void)time;
    (void)particle_diags;
    (void)nlev;
    (void)prefix;
    (void)plot_raw_fields;
    (void)plot_raw_fields_guards;
    (void)plot_raw_rho;
    (void)plot_raw_F;
#else
    amrex::Vector<amrex::MultiFab> *mf_ptr =
        const_cast<amrex::Vector<amrex::MultiFab>*>(&mf);

    if (m_insitu_bridge->update(iteration[0], time, m_amr_mesh,
        {mf_ptr}, {varnames}))
    {
        amrex::ErrorStream() << "FlushFormatSensei::WriteToFile : "
            "Failed to update the in situ bridge." << std::endl;

        amrex::Abort();
    }
#endif
}

void
FlushFormatSensei::WriteParticles (
    const amrex::Vector<ParticleDiag>& particle_diags) const
{
#ifndef BL_USE_SENSEI_INSITU
    (void)particle_diags;
#else
    amrex::ErrorStream() << "FlushFormatSensei::WriteParticles : "
        "Not yet implemented." << std::endl;

    amrex::Abort();
#endif
}
