#include "LaserEnvelope.H"

#include "Evolve/WarpXDtType.H"
#include "Laser/LaserProfiles.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"

using namespace amrex;
using namespace WarpXLaserProfiles;

LaserEnvelope::LaserEnvelope (const int nlevs_max){
    ReadParameters();
    AllocateMFs(nlevs_max);
}

void LaserEnvelope::ReadParameters ()
{
    const ParmParse pp_laser_envelopes("laser_envelopes");
    std::string laser_name;
    pp_laser_envelopes.query("names", laser_name);

    const ParmParse pp_laser_name(laser_name);
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

    amrex::Print() << "The wavelength of the laser is " << m_wavelength << " nm\n";
    amrex::Print() << "The profile waist of the laser is " << m_profile_waist << " s\n";
    amrex::Print() << "The profile peak of the laser is " << m_profile_t_peak << " s\n";
    amrex::Print() << "We are considering the laser " << laser_name << " \n";
}

void LaserEnvelope::AllocateMFs (const int nlevs_max)
{
    A_laser_envelope.resize(nlevs_max);
}

void LaserEnvelope::AllocateLevelMFs (
    const int lev,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    const int ncomps,
    const amrex::IntVect& ngA,
    const amrex::IntVect& A_nodal_flag)
{
    auto & warpx = WarpX::GetInstance();

    warpx.AllocInitMultiFab(A_laser_envelope[lev], amrex::convert(ba, A_nodal_flag),
        dm, ncomps, ngA, lev, "A_laser_envelope", 0.0_rt);
}

void LaserEnvelope::ClearLevel (const int lev)
{
    A_laser_envelope[lev].reset();
}
