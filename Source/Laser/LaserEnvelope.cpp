#include "LaserEnvelope.H"

#include "Laser/LaserProfiles.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpX_Complex.H"

#include <AMReX_BLassert.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <cstdlib>
#include <numeric>
#include <vector>


using namespace amrex;

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

    // Read laser parameters:
    // * wavelength
    // * e_max or a0
    utils::parser::getWithParser(pp_laser_name, "wavelength", m_wavelength);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_wavelength > 0._rt, "The laser wavelength must be > 0.");
    const bool e_max_is_specified =
        utils::parser::queryWithParser(pp_laser_name, "e_max", m_e_max);
    amrex::Real a0;
    const bool a0_is_specified =
        utils::parser::queryWithParser(pp_laser_name, "a0", a0);
    if (a0_is_specified)
    {
        const amrex::Real omega = 2._rt*MathConst::pi*PhysConst::c/m_wavelength;
        m_e_max = PhysConst::m_e * omega * PhysConst::c * a0 / PhysConst::q_e;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        e_max_is_specified || a0_is_specified,
        "Exactly one of e_max or a0 must be specified for the laser.\n");

    if (laser_type_s == "gaussian")
    {
        utils::parser::getWithParser(pp_laser_name, "profile_waist", m_waist);
        utils::parser::getWithParser(pp_laser_name, "profile_duration", m_duration);
        utils::parser::getWithParser(pp_laser_name, "profile_z_peak", m_z_peak);
        utils::parser::getWithParser(pp_laser_name, "profile_focal_distance", m_focal_distance);
        utils::parser::queryWithParser(pp_laser_name, "zeta", m_zeta);
        utils::parser::queryWithParser(pp_laser_name, "beta", m_beta);
        utils::parser::queryWithParser(pp_laser_name, "phi2", m_phi2);
        utils::parser::queryWithParser(pp_laser_name, "phi0", m_phi0);
        utils::parser::queryArrWithParser(pp_laser_name, "direction", m_direction);
    }

    amrex::Print() << "The wavelength of the laser is " << m_wavelength << " nm\n";
    amrex::Print() << "The profile waist of the laser is " << m_waist << " s\n";
    amrex::Print() << "The profile peak of the laser is " << m_z_peak << " s\n";
    amrex::Print() << "We are considering the laser " << laser_name << " \n";
    amrex::Print() << "Tue duration " << m_profile_duration << " \n";
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

Complex LaserEnvelope::FillAmplitude (const amrex::Real x, const amrex::Real y, const amrex::Real z)
{
    Complex amplitude = Complex{0._rt, 0._rt};

    const Complex I(0,1);
    // Calculate a few factors which are independent of the macroparticle
    const Real k0 = 2._rt*MathConst::pi/m_wavelength;
    const Real inv_tau2 = 1._rt /(m_duration * m_duration);
    const Real oscillation_phase = k0 * ( - m_z_peak ) + m_phi0;
    // The coefficients below contain info about Gouy phase,
    // laser diffraction, and phase front curvature
    const Complex diffract_factor =
        1._rt + I * m_focal_distance * 2._rt/
        ( k0 * m_waist * m_waist );
    const Complex inv_complex_waist_2 =
        1._rt /(m_waist*m_waist * diffract_factor );

    // Time stretching due to STCs and phi2 complex envelope
    // (1 if zeta=0, beta=0, phi2=0)
    const Complex stretch_factor = 1._rt + 4._rt *
        (m_zeta+m_beta*m_focal_distance*inv_tau2)
        * (m_zeta+m_beta*m_focal_distance*inv_complex_waist_2)
        + 2._rt*I*(m_phi2-m_beta*m_beta*k0*m_focal_distance)*inv_tau2;

    // Amplitude and monochromatic oscillations
    const Complex t_prefactor =
        m_e_max * amrex::exp( I * oscillation_phase );

    // Because diffract_factor is a complex, the code below takes into
    // account the impact of the dimensionality on both the Gouy phase
    // and the amplitude of the laser
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    const Complex prefactor = t_prefactor / diffract_factor;
#elif defined(WARPX_DIM_XZ)
    const Complex prefactor = t_prefactor / amrex::sqrt(diffract_factor);
#else
    const Complex prefactor = t_prefactor;
#endif

#if defined(WARPX_DIM_3D)
    auto arg = m_direction[0]*m_p_X[0] +
        m_direction[1]*m_p_X[1] +
        m_direction[2]*m_p_X[2];

    if (arg < -1.0_rt || arg > 1.0_rt)
        m_theta_stc = 0._rt;
    else
        m_theta_stc = std::acos(arg);
#else
    m_theta_stc = 0.;
#endif

    // Copy member variables to tmp copies for GPU runs.
    auto const tmp_profile_z_peak = m_z_peak;
    auto const tmp_beta = m_beta;
    auto const tmp_zeta = m_zeta;
    auto const tmp_theta_stc = m_theta_stc;
    auto const tmp_profile_focal_distance = m_focal_distance;
    // Loop through the macroparticle to calculate the proper amplitude
    const Complex stc_exponent = 1._rt / stretch_factor * inv_tau2 *
        amrex::pow((- tmp_profile_z_peak / PhysConst::c -
        tmp_beta*k0*(x*std::cos(tmp_theta_stc) + y*std::sin(tmp_theta_stc)) -
        2._rt *I*(x*std::cos(tmp_theta_stc) + y*std::sin(tmp_theta_stc))
        *( tmp_zeta - tmp_beta*tmp_profile_focal_distance ) * inv_complex_waist_2),2);
    // stcfactor = everything but complex transverse envelope
    const Complex stcfactor = prefactor * amrex::exp( - stc_exponent );
    // Exp argument for transverse envelope
    const Complex exp_argument = - (x*x + y*y) * inv_complex_waist_2;
    // stcfactor + transverse envelope
    amplitude = ( stcfactor  * amrex::exp( exp_argument )).real();

    return amplitude;
}

void LaserEnvelope::InitData (const int finestLevel)
{
    // TODO
    // Put all the body function inside a loop over levels
    for (int lev = 0; lev <= finestLevel; ++lev) {
    // Loop through the grids, and over the tiles within each grid
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*A_laser_envelope[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        WarpX &warpx = WarpX::GetInstance();
        const amrex::Geometry &geom = warpx.Geom(lev);
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        // Extract field data for this grid/tile
        amrex::Array4<amrex::Real> const& A_laser_envelope_arr = A_laser_envelope[lev]->array(mfi);

        // Extract tileboxes for which to loop
        const amrex::Box& tilebox  = mfi.tilebox();

        ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            #if defined(WARPX_DIM_3D)
                amrex::Real x = problo[0] + i * dx[0];
                amrex::Real y = problo[1] + j * dx[1];
                amrex::Real z = problo[2] + k * dx[2];
            #elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                amrex::Real x = problo[0] + i * dx[0];
                amrex::Real y = 0.0_rt;
                amrex::Real z = problo[1] + j * dx[1];
            #else
                amrex::Real x = 0.0_rt;
                amrex::Real y = 0.0_rt;
                amrex::Real z = problo[0] + i * dx[0];
            #endif

            A_laser_envelope_arr(i,j,k,0) = FillAmplitude(x, y, z).real();
            A_laser_envelope_arr(i,j,k,1) = FillAmplitude(x, y, z).imag();

            //w_z =
            //Lz = std
            //a_T(i, j, k) = std::exp(-(std::pow(x, 2)+ std::pow(y, 2))/std::pow(w_z,2));
            //a_L(i, j, k) = a_0/(1+std::pow(z,2)) * std::exp(-std::pow((z-zf),2)/std::pow(Lz, 2))
            //A_arr(i, j, k) = a_T(i, j, k) * a_L(i, j, k);
        });
    }
}
}
