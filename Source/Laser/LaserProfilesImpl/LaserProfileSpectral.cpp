/* Copyright 2022 Ryan Sandberg
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"

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
    const amrex::ParmParse& ppc,
    CommonLaserParameters params)
{
    std::cout << "initializing spectral file\n";

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
    getWithParser(ppl, "profile_waist", m_params.waist);
    getWithParser(ppl, "profile_duration", m_params.duration);
    getWithParser(ppl, "profile_t_peak", m_params.t_peak);
    getWithParser(ppl, "profile_focal_distance", m_params.focal_distance);
    queryWithParser(ppl, "zeta", m_params.zeta);
    queryWithParser(ppl, "beta", m_params.beta);
    queryWithParser(ppl, "phi2", m_params.phi2);
    queryWithParser(ppl, "phi0", m_params.phi0);

    m_params.stc_direction = m_common_params.p_X;
    queryArrWithParser(ppl, "stc_direction", m_params.stc_direction);
    auto const s = 1.0_rt / std::sqrt(
        m_params.stc_direction[0]*m_params.stc_direction[0] +
        m_params.stc_direction[1]*m_params.stc_direction[1] +
        m_params.stc_direction[2]*m_params.stc_direction[2]);
    m_params.stc_direction = {
        m_params.stc_direction[0]*s,
        m_params.stc_direction[1]*s,
        m_params.stc_direction[2]*s };
    auto const dp2 =
        std::inner_product(
            m_common_params.nvec.begin(),
            m_common_params.nvec.end(),
            m_params.stc_direction.begin(), 0.0);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(dp2) < 1.0e-14,
        "stc_direction is not perpendicular to the laser plane vector");

    // Get angle between p_X and stc_direction
    // in 2d, stcs are in the simulation plane
#if defined(WARPX_DIM_3D )
    auto arg = m_params.stc_direction[0]*m_common_params.p_X[0] +
        m_params.stc_direction[1]*m_common_params.p_X[1] +
        m_params.stc_direction[2]*m_common_params.p_X[2];

    if (arg < -1.0_rt || arg > 1.0_rt)
        m_params.theta_stc = 0._rt;
    else
        m_params.theta_stc = std::acos(arg);
#else
    m_params.theta_stc = 0.;
#endif

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

        // int nt = 2001;
        // int nt_2 = 1001;
        // std::vector<Real> tarr(nt);
        // std::vector<Complex> env_arr(nt);
        // // std::iota(tarr.begin(),tarr.end(), -1000);
        // Real dt = 1e-16;
        // for (size_t ii = 0; ii < tarr.size(); ++ii) {
        //     tarr[ii] = dt * (ii - nt_2);
        // }
        // // std::transform(tarr.begin(), tarr.end(), [](){})
        // std::transform(tarr.begin(), tarr.end(), env_arr.begin(), [inv_tau2, simple_oscillation](Real t){ Complex env_exp = t*t*inv_tau2; return simple_oscillation * amrex::exp(-env_exp); });

        // t_min = tarr[0];
        // t_max = tarr[nt-1];


        // std::cout << "array:\n nt = " << nt << ", dt = " << dt << ", t_min = " << t_min << ", t_max = " << t_max << "\n";
        // std::cout << "loaded:\n nt = " << m_params.nt << ", dt = " << m_params.dt << ", t_min = " << m_params.tmin << ", t_max = " << m_params.tmax << "\n";

    for (size_t ii = 0; ii < extent[0]; ++ii) {
        m_params.E_data[ii] = Complex(x_data.get()[ii].real(), x_data.get()[ii].imag());
    }



        // Real error = 0;
        // for (int ii = 0; ii < nt; ++ii) {
        //     Real difference = abs(env_arr[ii] - m_params.E_data[ii]);
        //     if (difference > 0.1) {
        //         std::cout << "discrepancy at index " << ii << "\nE_arr = " << env_arr[ii] << ", Eload = " << m_params.E_data[ii] << "\n"; 
        //     }
        //     error += difference * difference;
        // }
        // std::cout << "E error = " << error << "\n";
}

// std::pair<int,int>
// WarpXLaserProfiles::FromTXYEFileLaserProfile::find_left_right_time_indices(amrex::Real t) const
// {
//     int idx_t_right;
//     // if(m_params.is_grid_uniform){
//         const auto t_min = m_params.t_coords.front();
//         const auto t_max = m_params.t_coords.back();
//         const auto temp_idx_t_right = static_cast<int>(
//             std::ceil( (m_params.nt-1)*(t-m_params.t_peak-t_min)/(t_max-t_min)));
//         idx_t_right = max(min(temp_idx_t_right, m_params.nt-1),1);
//     // }
//     // else{
//     //     idx_t_right = std::distance(m_params.t_coords.begin(),
//     //     std::upper_bound(m_params.t_coords.begin(),
//     //         m_params.t_coords.end(), t));
//     // }
//     return std::make_pair(idx_t_right-1, idx_t_right);
// }


// create fine spectral grid
// interpolate to fine spectral grid
// set fine temporal grid
// inverse FT to fine temporal grid

/* \brief compute field amplitude for a laser with user-provided spectrum, at particles' position
 *
 * Both Xp and Yp are given in laser plane coordinate.
 * For each particle with position Xp and Yp, this routine computes the
 * amplitude of the laser electric field, stored in array amplitude.
 *
 * \param np: number of laser particles
 * \param Xp: pointer to first component of positions of laser particles
 * \param Yp: pointer to second component of positions of laser particles
 * \param t: Current physical time
 * \param amplitude: pointer to array of field amplitude.
 */

// void
// WarpXLaserProfiles::SpectralLaserProfile::fill_amplitude (
//     const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
//     Real t, Real * AMREX_RESTRICT const amplitude) const
// {
//     // std::cout << "Made it into spectral laser fill_amplitude\n";
//     Complex I(0,1);
//     // Calculate a few factors which are independent of the macroparticle
//     const Real k0 = 2._rt*MathConst::pi/m_common_params.wavelength;
//     const Real inv_tau2 = 1._rt /(m_params.duration * m_params.duration);
//     const Real oscillation_phase = k0 * PhysConst::c * ( t - m_params.t_peak ) + m_params.phi0;
//     // The coefficients below contain info about Gouy phase,
//     // laser diffraction, and phase front curvature
//     const Complex diffract_factor =
//         1._rt + I * m_params.focal_distance * 2._rt /
//         ( k0 * m_params.waist * m_params.waist );
//     const Complex inv_complex_waist_2 =
//         1._rt /(m_params.waist*m_params.waist * diffract_factor );
//     // std::cout << "parsed some parameters\n";

// // interpolate and apply


// // gaussian transverse profile * longitudinal profile provided by user
// // making the Paraxial approximation
//     const auto t_min = m_params.tmin;
//     const auto t_max = m_params.tmax;
//     const auto delta_t = t - m_params.t_peak;
//     if (delta_t < t_min || delta_t >= t_max) {
//         amrex::ParallelFor(
//             np,
//             [=] AMREX_GPU_DEVICE (int i) {
//                 amplitude[i] = 0.0_rt;
//             }
//         );
//     } else {
//         // const auto temp_idx_t_right = static_cast<int>(
//         //     std::ceil( (m_params.nt-1)*(delta_t-t_min)/(t_max-t_min)));
//         const auto idx_t_left_temp = static_cast<int>((delta_t-t_min) / m_params.dt);
//         // std::cout << "t_left < delta_t < t_right\n";
//         const auto t_left = t_min + m_params.dt * idx_t_left_temp;
//         // std::cout << t_left << " <? " << delta_t << " <? " << t_left + m_params.dt;
//         const int idx_t_left = max(min(idx_t_left_temp, m_params.nt-2),0);
//         // interpolate between idx_t_left and idx_t_right onto t to get E there
//         if (delta_t < t_left || delta_t >= t_left + m_params.dt) {
//             std::cout << "Out of bounds t slipped through!\n";
//         } 

//         // std::cout << "diffract factor: " << diffract_factor << "\n";
        
//         amrex::ParallelFor(
//             np,
//             [=] AMREX_GPU_DEVICE (int i) {
//                 if (idx_t_left < 0 || idx_t_left > m_params.nt -1) {
//                     std::cout << "out of range index slipped through!\n";
//                 } 
//                 const Complex longitudinal_profile = 
//                     ((delta_t-t_left)*m_params.E_data[idx_t_left]
//                      + (t_left+m_params.dt - delta_t)*m_params.E_data[idx_t_left+1])
//                      / m_params.dt;

                
//                 // const Complex stc_exponent = 1._rt / stretch_factor * inv_tau2 *
//                 //     amrex::pow((t - tmp_profile_t_peak -
//                 //         tmp_beta*k0*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc)) -
//                 //         2._rt *I*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc))
//                 //         *( tmp_zeta - tmp_beta*tmp_profile_focal_distance ) * inv_complex_waist_2),2);
                        
//                 // // stcfactor = everything but complex transverse envelope
//                 // const Complex stcfactor = prefactor * amrex::exp( - stc_exponent );
                
//                 // Exp argument for transverse envelope
//                 const Complex exp_argument = - ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_complex_waist_2;
//                 // stcfactor + transverse envelope
//                 // amplitude[i] = ( stcfactor * amrex::exp( exp_argument ) ).real();


//     #if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
//                 const Complex transverse_profile = amrex::exp( exp_argument ) / diffract_factor;
//     #elif defined(WARPX_DIM_XZ)
//                 const Complex transverse_profile = amrex::exp( exp_argument ) / amrex::sqrt(diffract_factor);
//     #endif
//                 amplitude[i] = (m_common_params.e_max * longitudinal_profile * transverse_profile).real();
//             }
//         );
//     }

// }

// void
// WarpXLaserProfiles::SpectralLaserProfile::fill_amplitude (
//     const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
//     Real t, Real * AMREX_RESTRICT const amplitude) const
// {
//     Complex I(0,1);
//     // Calculate a few factors which are independent of the macroparticle
//     const Real k0 = 2._rt*MathConst::pi/m_common_params.wavelength;
//     const Real inv_tau2 = 1._rt /(m_params.duration * m_params.duration);
//     const Real oscillation_phase = k0 * PhysConst::c * ( t - m_params.t_peak ) + m_params.phi0;
//     // The coefficients below contain info about Gouy phase,
//     // laser diffraction, and phase front curvature
//     const Complex diffract_factor =
//         1._rt + I * m_params.focal_distance * 2._rt/
//         ( k0 * m_params.waist * m_params.waist );
//     const Complex inv_complex_waist_2 =
//         1._rt /(m_params.waist*m_params.waist * diffract_factor );

//     // Time stretching due to STCs and phi2 complex envelope
//     // (1 if zeta=0, beta=0, phi2=0)
//     const Complex stretch_factor = 1._rt + 4._rt *
//         (m_params.zeta+m_params.beta*m_params.focal_distance*inv_tau2)
//         * (m_params.zeta+m_params.beta*m_params.focal_distance*inv_complex_waist_2)
//         + 2._rt*I*(m_params.phi2-m_params.beta*m_params.beta*k0*m_params.focal_distance)*inv_tau2;

//     // Amplitude and monochromatic oscillations
//     Complex prefactor =
//         // m_common_params.e_max * amrex::exp( I * oscillation_phase );
//         m_common_params.e_max;

//     // Because diffract_factor is a complex, the code below takes into
//     // account the impact of the dimensionality on both the Gouy phase
//     // and the amplitude of the laser
// #if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
//     prefactor = prefactor / diffract_factor;
// #elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
//     prefactor = prefactor / amrex::sqrt(diffract_factor);
// #endif

//     // Copy member variables to tmp copies for GPU runs.
//     auto const tmp_profile_t_peak = m_params.t_peak;
//     auto const tmp_beta = m_params.beta;
//     auto const tmp_zeta = m_params.zeta;
//     auto const tmp_theta_stc = m_params.theta_stc;
//     auto const tmp_profile_focal_distance = m_params.focal_distance;
//     const auto t_min = m_params.tmin;
//     const auto t_max = m_params.tmax;
//     const auto delta_t = t - m_params.t_peak;
//     if (delta_t < t_min || delta_t >= t_max) {
//         amrex::ParallelFor(
//             np,
//             [=] AMREX_GPU_DEVICE (int i) {
//                 amplitude[i] = 0.0_rt;
//             }
//         );
//     } else {

//         // const auto temp_idx_t_right = static_cast<int>(
//         //     std::ceil( (m_params.nt-1)*(delta_t-t_min)/(t_max-t_min)));
//         const auto idx_t_left_temp = static_cast<int>((delta_t-t_min) / m_params.dt);
//         // std::cout << "t_left < delta_t < t_right\n";
//         const auto t_left = t_min + m_params.dt * idx_t_left_temp;
//         // std::cout << t_left << " <? " << delta_t << " <? " << t_left + m_params.dt;
//         const int idx_t_left = max(min(idx_t_left_temp, m_params.nt-2),0);
//         // interpolate between idx_t_left and idx_t_right onto t to get E there
//         if (delta_t < t_left || delta_t >= t_left + m_params.dt) {
//             std::cout << "Out of bounds t slipped through!\n";
//         } 
//     // Loop through the macroparticle to calculate the proper amplitude
//         amrex::ParallelFor(
//             np,
//             [=] AMREX_GPU_DEVICE (int i) {
                
//                 const Complex stc_exponent = 1._rt / stretch_factor * inv_tau2 *
//                     amrex::pow((t - tmp_profile_t_peak -
//                         tmp_beta*k0*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc)) -
//                         2._rt *I*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc))
//                         *( tmp_zeta - tmp_beta*tmp_profile_focal_distance ) * inv_complex_waist_2),2);
//                 // stcfactor = everything but complex transverse envelope
//                 const Complex stcfactor = prefactor * amrex::exp( - stc_exponent );
//                 // const Complex stcfactor = 
//                 //     ((delta_t-t_left)*m_params.E_data[idx_t_left]
//                 //      + (t_left+m_params.dt - delta_t)*m_params.E_data[idx_t_left+1])
//                 //      / m_params.dt;
//                 // Exp argument for transverse envelope
//                 const Complex exp_argument = - ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_complex_waist_2;
//                 // stcfactor + transverse envelope
//                 amplitude[i] = ( stcfactor * amrex::exp( exp_argument ) ).real();
//             }
//         );
//     }
// }

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
    std::cout << "phi0 = " << m_params.phi0 << "\n";
    std::cout << "wavelen = " << m_common_params.wavelength << "\n";
    // The coefficients below contain info about Gouy phase,
    // laser diffraction, and phase front curvature
    const Complex diffract_factor =
        1._rt + I * m_params.focal_distance * 2._rt/
        ( k0 * m_params.waist * m_params.waist );
    const Complex inv_complex_waist_2 =
        1._rt /(m_params.waist*m_params.waist * diffract_factor );

    // Time stretching due to STCs and phi2 complex envelope
    // (1 if zeta=0, beta=0, phi2=0)
    const Complex stretch_factor = 1._rt + 4._rt *
        (m_params.zeta+m_params.beta*m_params.focal_distance*inv_tau2)
        * (m_params.zeta+m_params.beta*m_params.focal_distance*inv_complex_waist_2)
        + 2._rt*I*(m_params.phi2-m_params.beta*m_params.beta*k0*m_params.focal_distance)*inv_tau2;

    // Amplitude and monochromatic oscillations
    Complex simple_oscillation = amrex::exp( I * oscillation_phase );
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
    bool use_external = true;
    int nt, nt_2;
    std::vector<Complex> env_arr;
    Real dt, t_min, t_max;
    if (use_external) {
        dt = m_params.dt;
        nt = m_params.nt;
        nt_2 = int(nt/2);
        t_min = m_params.tmin;
        t_max = m_params.tmax;
        for (int ii = 0; ii < nt; ii++) {
            env_arr.push_back(m_params.E_data[ii]);
        }

    } else {
        nt = 2001;
        nt_2 = 1000;
        std::vector<Real> tarr(nt);
        env_arr = std::vector<Complex>(nt);
        // std::iota(tarr.begin(),tarr.end(), -1000);
        dt = 1e-16;
        for (int ii = 0; ii < nt; ++ii) {
            tarr[ii] = dt * (ii - nt_2);
        }
        // std::transform(tarr.begin(), tarr.end(), [](){})
        std::transform(tarr.begin(), tarr.end(), env_arr.begin(), [inv_tau2, simple_oscillation](Real t){ Complex env_exp = t*t*inv_tau2; 
        return simple_oscillation * amrex::exp(-env_exp); 
        // return amrex::exp(-env_exp); 
        });

        t_min = tarr[0];
        t_max = tarr[nt-1];


        // std::cout << "array:\n nt = " << nt << ", dt = " << dt << ", t_min = " << t_min << ", t_max = " << t_max << "\n";
        // std::cout << "loaded:\n nt = " << m_params.nt << ", dt = " << m_params.dt << ", t_min = " << m_params.tmin << ", t_max = " << m_params.tmax << "\n";

        // Real error = 0;
        // for (int ii = 0; ii < nt; ++ii) {
        //     Real difference = abs(env_arr[ii] - m_params.E_data[ii]);
        //     if (difference > 0.1) {
        //         std::cout << "discrepancy at index " << ii << "\nE_arr = " << env_arr[ii] << ", Eload = " << m_params.E_data[ii] << "\n"; 
        //     }
        //     error += difference * difference;
        // }
        // std::cout << "E error = " << error << "\n";
    }

    const auto delta_t =  m_params.t_peak - t; // this is negative of what I think it should be and I can't explain why this works
    if (delta_t < t_min || delta_t >= t_max) {
        amrex::ParallelFor(
            np,
            [=] AMREX_GPU_DEVICE (int i) {
                amplitude[i] = 0.0_rt;
            }
        );
    } else {
        // const auto temp_idx_t_right = static_cast<int>(
        //     std::ceil( (m_params.nt-1)*(delta_t-t_min)/(t_max-t_min)));
        const auto idx_t_left_temp = static_cast<int>((delta_t-t_min) / dt);
        // std::cout << "t_left < delta_t < t_right\n";
        const auto t_left = t_min + dt * idx_t_left_temp;
        // std::cout << t_left << " <? " << delta_t << " <? " << t_left + m_params.dt;
        const int idx_t_left = max(min(idx_t_left_temp, nt-2),0);
        // interpolate between idx_t_left and idx_t_right onto t to get E there
        if (delta_t < t_left || delta_t >= t_left + dt) {
            std::cout << "Out of bounds t slipped through!\n";
        } 

        // std::cout << "diffract factor: " << diffract_factor << "\n";
        

    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            // const Complex stc_exponent = 1._rt / stretch_factor * inv_tau2 *
                // amrex::pow((t - tmp_profile_t_peak -
                //     tmp_beta*k0*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc)) -
                //     2._rt *I*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc))
                //     *( tmp_zeta - tmp_beta*tmp_profile_focal_distance ) * inv_complex_waist_2),2);
            const Complex stc_exponent = inv_tau2 * (t - tmp_profile_t_peak) * (t - tmp_profile_t_peak);

            // const Complex t_envelope = amrex::exp( - stc_exponent );
            const Complex t_envelope = -1._rt *  // this minus sign seems necessary but I can't explain why
                        ((delta_t-t_left)*env_arr[idx_t_left]
                     + (t_left+dt - delta_t)*env_arr[idx_t_left+1])
                     / dt;

            // stcfactor = everything but complex transverse envelope
            // const Complex stcfactor = simple_oscillation * amrex::exp( - stc_exponent );
            // const Complex longitudinal = simple_oscillation * t_envelope;
            const Complex longitudinal = t_envelope;
            // Exp argument for transverse envelope
            const Complex exp_argument = - ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_complex_waist_2;
            // stcfactor + transverse envelope
            amplitude[i] = ( prefactor * longitudinal * amrex::exp( exp_argument ) ).real();
        }
        );
    }
}
