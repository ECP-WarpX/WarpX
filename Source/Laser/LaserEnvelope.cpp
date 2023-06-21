#include "LaserEnvelope.H"
#include <AMReX.H>

using namespace amrex;

void LaserEnvelope::AllocateLaserEnvelope (int lev, const BoxArray& ba, const DistributionMapping& dm, const IntVect& ngA, const bool aux_is_nodal)
{
    A_nodal_flag = IntVect::TheNodeVector();
    A_ncomps = ptr.2
    AllocInitMultiFab(rho_fp[lev], amrex::convert(ba, A_nodal_flag), dm, A_ncomps, ngA, tag("A_fp"), 0.0_rt);
}

public:
LaserEnvelope::LaserEnvelope(AmrCore* amr_core, int ispecies, const std::string& name)
    : LaserEnvelope(amr_core, ispecies),
      m_laser_name{name}
{
    charge = 1.0;
    mass = std::numeric_limits<Real>::max();

    const ParmParse pp_laser_name(m_laser_name);

    // Parse the type of laser profile and set the corresponding flag `profile`
    std::string laser_type_s;
    pp_laser_name.get("profile", laser_type_s);
    std::transform(laser_type_s.begin(), laser_type_s.end(), laser_type_s.begin(), ::tolower);

    utils::parser::getWithParser(pp_laser_name, "wavelength", m_wavelength);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_wavelength > 0, "The laser wavelength must be >0.");
    const bool e_max_is_specified =
        utils::parser::queryWithParser(pp_laser_name, "e_max", m_e_max);
    Real a0;
    const bool a0_is_specified =
        utils::parser::queryWithParser(pp_laser_name, "a0", a0);
    if (a0_is_specified){
        const Real omega = 2._rt*MathConst::pi*PhysConst::c/m_wavelength;
        m_e_max = PhysConst::m_e * omega * PhysConst::c * a0 / PhysConst::q_e;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        e_max_is_specified ^ a0_is_specified,
        "Exactly one of e_max or a0 must be specified for the laser.\n"
        );

    pp_laser_name.query("do_continuous_injection", do_continuous_injection);
    utils::parser::queryWithParser(pp_laser_name,
        "min_particles_per_mode", m_min_particles_per_mode);

    if (m_e_max == amrex::Real(0.)){
        ablastr::warn_manager::WMRecordWarning("Laser",
            m_laser_name + " with zero amplitude disabled.",
            ablastr::warn_manager::WarnPriority::low);
        m_enabled = false;
        return; // Disable laser if amplitude is 0
    }

    m_laser_envelope_model = std::make_unique<HybridPICModel>(nlevs_max);
    m_laser_envelope_model->AllocateLevelMFs(
            lev, ba, dm, ncomps, ngJ, ngRho, jx_nodal_flag, jy_nodal_flag,
            jz_nodal_flag, rho_nodal_flag
        );
}

