/* Copyright 2022 Ryan Sandberg
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
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

#include <openPMD/openPMD.hpp>

#include <cmath>
#include <cstdlib>
#include <memory>   // std::shared_ptr
#include <numeric>
#include <vector>

using namespace amrex;

void
WarpXLaserProfiles::SpectralLaserProfile::init (
    const amrex::ParmParse& ppl,
    CommonLaserParameters params)
{

    // Parse the spectral file name
    ppl.get("spectral_file_name", m_params.spectral_file_name);
    if(m_params.spectral_file_name.empty())
    {
        Abort("spectral_file_name must be provided for spectral-file laser profile!");
    }
    parse_spectral_file(m_params.spectral_file_name);

    // May want to interpolate and FT the spectral field here
    // May want to also determine waveleneth from spectral data and over-write

    //Copy common params
    m_common_params = params;

    // Parse the properties of the Gaussian profile
    utils::parser::getWithParser(ppl, "profile_waist", m_params.waist);
    // getWithParser(ppl, "profile_duration", m_params.duration);
    utils::parser::getWithParser(ppl, "profile_t_peak", m_params.t_peak);
    utils::parser::getWithParser(ppl, "profile_focal_distance", m_params.focal_distance);
    // queryWithParser(ppl, "zeta", m_params.zeta);
    // queryWithParser(ppl, "beta", m_params.beta);
    // queryWithParser(ppl, "phi2", m_params.phi2);
    // queryWithParser(ppl, "phi0", m_params.phi0);

    // m_params.stc_direction = m_common_params.p_X;
    // queryArrWithParser(ppl, "stc_direction", m_params.stc_direction);
    // auto const s = 1.0_rt / std::sqrt(
    //     m_params.stc_direction[0]*m_params.stc_direction[0] +
    //     m_params.stc_direction[1]*m_params.stc_direction[1] +
    //     m_params.stc_direction[2]*m_params.stc_direction[2]);
    // m_params.stc_direction = {
    //     m_params.stc_direction[0]*s,
    //     m_params.stc_direction[1]*s,
    //     m_params.stc_direction[2]*s };
    // auto const dp2 =
    //     std::inner_product(
    //         m_common_params.nvec.begin(),
    //         m_common_params.nvec.end(),
    //         m_params.stc_direction.begin(), 0.0);
    // WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(dp2) < 1.0e-14,
    //     "stc_direction is not perpendicular to the laser plane vector");

    // Get angle between p_X and stc_direction
    // in 2d, stcs are in the simulation plane
// #if defined(WARPX_DIM_3D )
//     auto arg = m_params.stc_direction[0]*m_common_params.p_X[0] +
//         m_params.stc_direction[1]*m_common_params.p_X[1] +
//         m_params.stc_direction[2]*m_common_params.p_X[2];

//     if (arg < -1.0_rt || arg > 1.0_rt)
//         m_params.theta_stc = 0._rt;
//     else
//         m_params.theta_stc = std::acos(arg);
// #else
    // m_params.theta_stc = 0.;
// #endif

}

void
WarpXLaserProfiles::SpectralLaserProfile::parse_spectral_file(std::string spectral_file_name)
{

    openPMD::Series series =
        openPMD::Series(spectral_file_name, openPMD::Access::READ_ONLY);
    auto i = series.iterations[1];

    // record
    auto E = i.meshes["E"];

    // record components
    auto E_x = E["x"];

    std::shared_ptr< std::complex<double> > x_data =
            E_x.loadChunk< std::complex<double> >();

    series.flush();
    auto extent = E_x.getExtent();
    auto offset = E.gridGlobalOffset();
    auto spacing = E.gridSpacing<Real>();

    m_params.nt = extent[0];
    m_params.tmin = offset[0];
    m_params.dt = spacing[0];
    m_params.tmax = m_params.tmin + m_params.dt * (m_params.nt-1);

    m_params.E_data.resize(extent[0]);

    for (size_t ii = 0; ii < extent[0]; ++ii) {
        m_params.E_data[ii] = Complex(x_data.get()[ii].real(), x_data.get()[ii].imag());
    }
}


void
WarpXLaserProfiles::SpectralLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    Complex I(0,1);
    // Calculate a few factors which are independent of the macroparticle
    const Real k0 = 2._rt*MathConst::pi/m_common_params.wavelength;
    // const Real inv_tau2 = 1._rt /(m_params.duration * m_params.duration);

    Complex longitudinal;

    const Real dt = m_params.dt;
    const int nt = m_params.nt;
    const Real t_min = m_params.tmin;
    const Real t_max = m_params.tmax;

    const auto delta_t =  m_params.t_peak - t; // this is negative of what I think it should be and I can't explain why this works
    if (delta_t < t_min || delta_t >= t_max)  {
        longitudinal = 0.0;
    } else {

        const auto idx_t_left_temp = static_cast<int>((delta_t-t_min) / dt);
        const auto t_left = t_min + dt * idx_t_left_temp;
        const int idx_t_left = max(min(idx_t_left_temp, nt-2),0);
        // interpolate between idx_t_left and idx_t_right onto t to get E there
        if (delta_t < t_left || delta_t >= t_left + dt) {
            longitudinal = 0;
            std::cout << "Out of bounds t slipped through!\n";
        } else {
            longitudinal = -1._rt *  // this minus sign seems necessary but I can't explain why
            ((delta_t-t_left)*m_params.E_data[idx_t_left]
            + (t_left+dt - delta_t)*m_params.E_data[idx_t_left+1])
            / dt;
        }
    }

    // The coefficients below contain info about Gouy phase,
    // laser diffraction, and phase front curvature
    const Complex diffract_factor =
        1._rt + I * m_params.focal_distance * 2._rt/
        ( k0 * m_params.waist * m_params.waist );
    const Complex inv_complex_waist_2 =
        1._rt /(m_params.waist*m_params.waist * diffract_factor );

    // Time stretching due to STCs and phi2 complex envelope
    // (1 if zeta=0, beta=0, phi2=0)
    // const Complex stretch_factor = 1._rt;

    // Amplitude and monochromatic oscillations
    Complex emax_and_diffract =
        m_common_params.e_max;

    // Because diffract_factor is a complex, the code below takes into
    // account the impact of the dimensionality on both the Gouy phase
    // and the amplitude of the laser
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    emax_and_diffract = emax_and_diffract / diffract_factor;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    emax_and_diffract = emax_and_diffract / amrex::sqrt(diffract_factor);
#endif

    // Copy member variables to tmp copies for GPU runs.
    // auto const tmp_profile_t_peak = m_params.t_peak;
    // auto const tmp_beta = m_params.beta;
    // auto const tmp_zeta = m_params.zeta;
    // auto const tmp_theta_stc = m_params.theta_stc;
    // auto const tmp_profile_focal_distance = m_params.focal_distance;

    auto const tmp_emax_diffract = emax_and_diffract;
    auto const tmp_long = longitudinal;

    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            // Exp argument for transverse envelope
            const Complex exp_argument = - ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_complex_waist_2;
            // stcfactor + transverse envelope
            amplitude[i] = ( tmp_emax_diffract * tmp_long * amrex::exp( exp_argument ) ).real();
        }
        );
}


/*
void
WarpXLaserProfiles::SpectralLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    Complex I(0,1);
    // Calculate a few factors which are independent of the macroparticle
    const Real k0 = 2._rt*MathConst::pi/m_common_params.wavelength;
    const Real inv_tau2 = 1._rt /(m_params.duration * m_params.duration);
    const Real oscillation_phase = k0 * PhysConst::c * ( t - m_params.t_peak ) + m_params.phi0;
    const Complex simple_oscillation = amrex::exp( I * oscillation_phase );
    const Complex stc_exponent = inv_tau2 *
        (t - m_params.t_peak)*(t - m_params.t_peak);

    bool use_gaussian = false;
    Complex longitudinal;
    if (use_gaussian) {
        // stcfactor = everything but complex transverse envelope
        longitudinal = simple_oscillation * amrex::exp( - stc_exponent );
    } else {
        bool use_external = true;
        int nt, nt_2;
        Real dt, t_min, t_max;
        std::vector<Complex> e_arr;

        if (use_external) {
            dt = m_params.dt;
            nt = m_params.nt;
            nt_2 = int(nt/2);
            t_min = m_params.tmin;
            t_max = m_params.tmax;
            for (int ii = 0; ii < nt; ii++) {
                e_arr.push_back(m_params.E_data[ii]);
            }
        } else {
            nt = 2001;
            nt_2 = 1000;
            std::vector<Real> tarr(nt);
            e_arr = std::vector<Complex>(nt);
            // std::iota(tarr.begin(),tarr.end(), -1000);
            dt = 1e-16;
            for (int ii = 0; ii < nt; ++ii) {
                tarr[ii] = dt * (ii - nt_2);
            }
            t_min = tarr[0];
            t_max = tarr[nt-1];
            // std::transform(tarr.begin(), tarr.end(), [](){})
            std::transform(tarr.begin(), tarr.end(), e_arr.begin(), [inv_tau2, simple_oscillation](Real t){ Complex env_exp = t*t*inv_tau2;
                return simple_oscillation * amrex::exp(-env_exp);
            // return amrex::exp(-env_exp);
            });
        }

        const auto delta_t =  m_params.t_peak - t; // this is negative of what I think it should be and I can't explain why this works
        if (delta_t < t_min || delta_t >= t_max)  {
            longitudinal = 0.0;
        } else {

            const auto idx_t_left_temp = static_cast<int>((delta_t-t_min) / dt);
            // std::cout << "t_left < delta_t < t_right\n";
            const auto t_left = t_min + dt * idx_t_left_temp;
            // std::cout << t_left << " <? " << delta_t << " <? " << t_left + m_params.dt;
            const int idx_t_left = max(min(idx_t_left_temp, nt-2),0);
            // interpolate between idx_t_left and idx_t_right onto t to get E there
            if (delta_t < t_left || delta_t >= t_left + dt) {
                longitudinal = 0;
                std::cout << "Out of bounds t slipped through!\n";
            } else {
                longitudinal = -1._rt *  // this minus sign seems necessary but I can't explain why
                ((delta_t-t_left)*e_arr[idx_t_left]
                + (t_left+dt - delta_t)*e_arr[idx_t_left+1])
                / dt;
            }
        }

    }
    // The coefficients below contain info about Gouy phase,
    // laser diffraction, and phase front curvature
    const Complex diffract_factor =
        1._rt + I * m_params.focal_distance * 2._rt/
        ( k0 * m_params.waist * m_params.waist );
    const Complex inv_complex_waist_2 =
        1._rt /(m_params.waist*m_params.waist * diffract_factor );

    // Time stretching due to STCs and phi2 complex envelope
    // (1 if zeta=0, beta=0, phi2=0)
    // const Complex stretch_factor = 1._rt;

    // Amplitude and monochromatic oscillations
    Complex prefactor =
        m_common_params.e_max;

    // Because diffract_factor is a complex, the code below takes into
    // account the impact of the dimensionality on both the Gouy phase
    // and the amplitude of the laser
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    prefactor = prefactor / diffract_factor;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    prefactor = prefactor / amrex::sqrt(diffract_factor);
#endif

    // Copy member variables to tmp copies for GPU runs.
    auto const tmp_profile_t_peak = m_params.t_peak;
    auto const tmp_beta = m_params.beta;
    auto const tmp_zeta = m_params.zeta;
    auto const tmp_theta_stc = m_params.theta_stc;
    auto const tmp_profile_focal_distance = m_params.focal_distance;

    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            // Exp argument for transverse envelope
            const Complex exp_argument = - ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_complex_waist_2;
            // stcfactor + transverse envelope
            amplitude[i] = ( prefactor *  longitudinal * amrex::exp( exp_argument ) ).real();
        }
        );
}
*/
