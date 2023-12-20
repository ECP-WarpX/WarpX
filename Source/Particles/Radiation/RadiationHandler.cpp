/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "RadiationHandler.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/TextMsg.H"

#include <ablastr/constant.H>
#include <ablastr/math/Utils.H>

#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#include <tuple>
#include <vector>

using namespace amrex;
using namespace ablastr::math;

namespace
{
    auto compute_detector_positions(
        const amrex::Array<amrex::Real,3>& center,
        const amrex::Array<amrex::Real,3>& direction,
        const amrex::Real distance,
        const amrex::Array<amrex::Real,3>& orientation,
        const amrex::Array<int,2>& det_points,
        const amrex::Array<amrex::Real,2>& theta_range)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            direction[0]*orientation[0] +
            direction[1]*orientation[1] +
            direction[2]*orientation[2] != 0,
            "Radiation detector orientation cannot be aligned with detector direction");

        auto u = amrex::Array<amrex::Real,3>{
            orientation[0] - direction[0]*orientation[0],
            orientation[1] - direction[1]*orientation[1],
            orientation[2] - direction[2]*orientation[2]};
        const auto one_over_u = 1.0_rt/std::sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        u[0] *= one_over_u;
        u[1] *= one_over_u;
        u[2] *= one_over_u;

        auto v = amrex::Array<amrex::Real,3>{
            direction[1]*u[2]-direction[2]*u[1],
            direction[2]*u[0]-direction[0]*u[2],
            direction[0]*u[1]-direction[1]*u[0]};
        const auto one_over_v = 1.0_rt/std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        v[0] *= one_over_v;
        v[1] *= one_over_v;
        v[2] *= one_over_v;

        auto us = amrex::Vector<amrex::Real>(det_points[0]);
        const auto ulim = distance*std::tan(theta_range[0]*0.5_rt);
        ablastr::math::LinSpaceFill(-ulim,ulim, det_points[0], us.begin());

        auto vs = amrex::Vector<amrex::Real>(det_points[1]);
        const auto vlim = distance*std::tan(theta_range[1]*0.5_rt);
        ablastr::math::LinSpaceFill(-vlim, vlim, det_points[1], vs.begin());

        const auto how_many = det_points[0]*det_points[1];

        auto det_x = amrex::Gpu::DeviceVector<amrex::Real>(how_many);
        auto det_y = amrex::Gpu::DeviceVector<amrex::Real>(how_many);
        auto det_z = amrex::Gpu::DeviceVector<amrex::Real>(how_many);

        const auto one_over_direction = 1.0_rt/std::sqrt(
            direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
        const auto norm_direction = amrex::Array<amrex::Real,3>{
            direction[0]*one_over_direction,
            direction[1]*one_over_direction,
            direction[2]*one_over_direction};

        for (int i = 0; i < det_points[0]; ++i)
        {
            for (int j = 0; j < det_points[1]; ++j)
            {
                det_x[i*det_points[1] + j] = center[0] + distance * norm_direction[0] + us[i]*u[0] + vs[j]*v[0];
                det_x[i*det_points[1] + j] = center[1] + distance * norm_direction[1] + us[i]*u[1] + vs[j]*v[1];
                det_x[i*det_points[1] + j] = center[2] + distance * norm_direction[2] + us[i]*u[2] + vs[j]*v[2];
            }
        }

        amrex::Gpu::synchronize();

        return std::make_tuple(det_x, det_y, det_z);

    }
}


RadiationHandler::RadiationHandler(const amrex::Array<amrex::Real,3>& center)
{
#if defined  WARPX_DIM_RZ
    WARPX_ABORT_WITH_MESSAGE("Radiation is not supported yet with RZ.");
#endif

    // Read in radiation input
    const amrex::ParmParse pp_radiation("radiation");

    //Resolution in frequency of the detector
    auto omega_range = std::vector<amrex::Real>(2);
    pp_radiation.getarr("omega_range", omega_range);
    std::copy(omega_range.begin(), omega_range.end(), m_omega_range.begin());
    pp_radiation.get("omega_points", m_omega_points);

    //Angle theta AND phi
    auto theta_range = std::vector<amrex::Real>(2);
    pp_radiation.getarr("angle_aperture", theta_range);
    std::copy(theta_range.begin(), theta_range.end(), m_theta_range.begin());

    //Detector parameters
    auto det_pts = std::vector<int>(2);
    pp_radiation.getarr("detector_number_points", det_pts);
    std::copy(det_pts.begin(), det_pts.end(), m_det_pts.begin());

    auto det_direction = std::vector<amrex::Real>(3);
    pp_radiation.getarr("detector_direction", det_direction);
    std::copy(det_direction.begin(), det_direction.end(), m_det_direction.begin());

    auto det_orientation = std::vector<amrex::Real>(3);
    pp_radiation.getarr("detector_orientation", det_orientation);
    std::copy(det_orientation.begin(), det_orientation.end(), m_det_orientation.begin());

    pp_radiation.get("detector_distance", m_det_distance);

    std::tie(det_pos_x, det_pos_y, det_pos_z) = compute_detector_positions(
        center, m_det_direction, m_det_distance,
        m_det_orientation, m_det_pts, m_theta_range);

    constexpr auto ncomp = 3;
    m_radiation_data = amrex::Gpu::DeviceVector<ablastr::math::Complex>(m_det_pts[0]*m_det_pts[1]*m_omega_points*ncomp);

    m_omegas = amrex::Gpu::DeviceVector<amrex::Real>(m_omega_points);
    ablastr::math::LinSpaceFill(m_omega_range[0], m_omega_range[1], m_omega_points, m_omegas.begin());
}


void RadiationHandler::add_radiation_contribution
    (const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc, amrex::Real current_time)
    {

        for (int lev = 0; lev <= pc->finestLevel(); ++lev)
        {
#ifdef AMREX_USE_OMP
            #pragma omp parallel
#endif
            {
                for (WarpXParIter pti(*pc, lev); pti.isValid(); ++pti)
                {
                    long const np = pti.numParticles();
                    const auto& attribs = pti.GetAttribs();
                    const auto& p_w = attribs[PIdx::w];
                    const auto& p_ux = attribs[PIdx::ux];
                    const auto& p_uy = attribs[PIdx::uy];
                    const auto& p_uz = attribs[PIdx::uz];

                    const auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                    auto& part=pc->GetParticles(lev)[index];
                    auto& soa = part.GetStructOfArrays();

                    const auto* p_ux_old = soa.GetRealData(pc->GetRealCompIndex("old_u_x")).data();
                    const auto* p_uy_old = soa.GetRealData(pc->GetRealCompIndex("old_u_y")).data();
                    const auto* p_uz_old = soa.GetRealData(pc->GetRealCompIndex("old_u_z")).data();

                    auto GetPosition = GetParticlePosition<PIdx>(pti);
                    auto const q = pc->getCharge();

                    const auto how_many_det_pos = static_cast<int>(det_pos_x.size());

                    const auto p_omegas = m_omegas.dataPtr();

                    const auto p_det_pos_x = det_pos_x.dataPtr();
                    const auto p_det_pos_y = det_pos_y.dataPtr();
                    const auto p_det_pos_z = det_pos_z.dataPtr();

                    auto p_radiation_data = m_radiation_data.dataPtr();

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip){
                        amrex::ParticleReal xp, yp, zp;
                        GetPosition.AsStored(ip,xp, yp, zp);

                        const auto ux = 0.5_prt*(p_ux[ip] + p_ux_old[ip]);
                        const auto uy = 0.5_prt*(p_uy[ip] + p_uy_old[ip]);
                        const auto uz = 0.5_prt*(p_uz[ip] + p_uz_old[ip]);
                        auto const u2 = ux*ux + uy*uy + uz*uz;

                        constexpr auto inv_c2 = 1._prt/(PhysConst::c*PhysConst::c);
                        auto const one_over_gamma = 1._prt/std::sqrt(1.0_rt + u2*inv_c2);
                        auto const one_over_gamma_c = one_over_gamma/PhysConst::c;
                        const auto bx = ux*one_over_gamma_c;
                        const auto by = uy*one_over_gamma_c;
                        const auto bz = uz*one_over_gamma_c;

                        const auto one_over_dt_gamma_c = one_over_gamma_c/dt;
                        const auto bpx = (p_ux[ip] - p_ux_old[ip])*one_over_dt_gamma_c;
                        const auto bpy = (p_uy[ip] - p_uy_old[ip])*one_over_dt_gamma_c;
                        const auto bpz = (p_uz[ip] - p_uz_old[ip])*one_over_dt_gamma_c;

                        const auto tot_q = q*p_w[ip];

                        for(int i_om=0; i_om < m_omega_points; ++i_om){

                            const auto i_omega_over_c = Complex{0.0_prt, 1.0_prt}*p_omegas[i_om]/PhysConst::c;

                            for (int i_det = 0; i_det < how_many_det_pos; ++i_det){

                                const auto part_det_x = xp - p_det_pos_x[i_det];
                                const auto part_det_y = yp - p_det_pos_y[i_det];
                                const auto part_det_z = zp - p_det_pos_z[i_det];
                                const auto d_part_det = std::sqrt(
                                    part_det_x*part_det_x + part_det_y*part_det_y + part_det_z*part_det_z);
                                const auto one_over_d_part_det = 1.0_prt/d_part_det;

                                const auto nx = part_det_x*one_over_d_part_det;
                                const auto ny = part_det_y*one_over_d_part_det;
                                const auto nz = part_det_z*one_over_d_part_det;

                                //Calculation of 1_beta.n, n corresponds to m_det_direction, the direction of the normal
                                const auto one_minus_b_dot_n = 1.0_prt - (bx*nx + by*ny + bz*nz);

                                const auto n_minus_beta_x = nx - bx;
                                const auto n_minus_beta_y = ny - by;
                                const auto n_minus_beta_z = nz - bz;

                                //Calculation of nxbeta
                                const auto n_minus_beta_cross_bp_x = n_minus_beta_y*bpz - n_minus_beta_z*bpy;
                                const auto n_minus_beta_cross_bp_y = n_minus_beta_z*bpx - n_minus_beta_x*bpz;
                                const auto n_minus_beta_cross_bp_z = n_minus_beta_x*bpy - n_minus_beta_y*bpx;

                                //Calculation of nxnxbeta
                                const auto n_cross_n_minus_beta_cross_bp_x =ny*n_minus_beta_cross_bp_z-nz*n_minus_beta_cross_bp_y;
                                const auto n_cross_n_minus_beta_cross_bp_y =nz*n_minus_beta_cross_bp_x-nx*n_minus_beta_cross_bp_z;
                                const auto n_cross_n_minus_beta_cross_bp_z =nx*n_minus_beta_cross_bp_y-ny*n_minus_beta_cross_bp_x;

                                const auto phase_term = amrex::exp(i_omega_over_c*(PhysConst::c*current_time - d_part_det));

                                const auto coeff = tot_q*phase_term/(one_minus_b_dot_n*one_minus_b_dot_n);

                                const auto cx = coeff*n_cross_n_minus_beta_cross_bp_x;
                                const auto cy = coeff*n_cross_n_minus_beta_cross_bp_y;
                                const auto cz = coeff*n_cross_n_minus_beta_cross_bp_z;

                                const int ncomp = 3;
                                const int idx0 = (i_om*how_many_det_pos + i_det)*ncomp;
                                const int idx1 = idx0 + 1;
                                const int idx2 = idx0 + 2;

#ifdef AMREX_USE_OMP
                                #pragma omp atomic
                                p_radiation_data[idx0].m_real += cx.m_real;
                                #pragma omp atomic
                                p_radiation_data[idx0].m_imag += cx.m_imag;
                                #pragma omp atomic
                                p_radiation_data[idx1].m_real += cy.m_real;
                                #pragma omp atomic
                                p_radiation_data[idx1].m_imag += cy.m_imag;
                                #pragma omp atomic
                                p_radiation_data[idx2].m_real += cz.m_real;
                                #pragma omp atomic
                                p_radiation_data[idx2].m_imag += cz.m_imag;
#else
                                p_radiation_data[idx0] += cx;
                                p_radiation_data[idx1] += cy;
                                p_radiation_data[idx2] += cz;

#endif
                            }
                        }
                    });

                }
            }
        }
    }


void RadiationHandler::gather_and_write_radiation(const std::string& filename)
{
    auto radiation_data_cpu = amrex::Vector<amrex::Real>(m_det_pts[0]*m_det_pts[1]*m_omega_points);
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost,
        m_radiation_data.begin(), m_radiation_data.end(), m_radiation_data.begin());
    amrex::Gpu::streamSynchronize();

    amrex::ParallelDescriptor::ReduceRealSum(radiation_data_cpu.data(), radiation_data_cpu.size());

    if (amrex::ParallelDescriptor::IOProcessor() ){
        auto of = std::ofstream(filename, std::ios::binary);

        for (const auto& d : radiation_data_cpu){
            of << static_cast<double>(d);
        }

        of.close();
    }
}

/*void RadiationHandler::Integral_overtime()
{
    for(int i_om=0; i_om<m_omega_points;i_om++){
                        for(int i_x=0; i_x<m_det_pts[0]; i_x++){
                            for(int i_y=0; i_y<m_det_pts[1]; i_y++){
                                 m_radiation_calculation[i_x,i_y,i_om]=std::pow(m_radiation_data[i_x,i_y,i_om,0],2)+std::pow(m_radiation_data[i_x,i_y,i_om,1],2)+std::pow(m_radiation_data[i_x,i_y,i_om,2],2)+std::pow(m_radiation_data[i_x,i_y,i_om,3],2)+std::pow(m_radiation_data[i_x,i_y,i_om,4],2)+std::pow(m_radiation_data[i_x,i_y,i_om,5],2);
                            }
                        }
    }
}*/
