/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "RadiationHandler.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Utils/TextMsg.H"

#include <ablastr/constant.H>
#include <ablastr/math/Utils.H>

#include <AMReX_ParmParse.H>

#include <vector>

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
        ablastr::math::LinSpaceFill(-ul,ul, det_points[0], us.being());

        auto vs = amrex::Vector<amrex::Real>(det_points[1]);
        const auto vlim = distance*std::tan(theta_range[1]*0.5_rt);
        ablastr::math::LinSpaceFill(-vlim, vlim, det_points[1], vs.being());

        const auto how_many = det_points[0]*det_points[1];

        auto det_x = amrex::GPU::DeviceVector<amrex::Real>(how_many);
        auto det_y = amrex::GPU::DeviceVector<amrex::Real>(how_many);
        auto det_z = amrex::GPU::DeviceVector<amrex::Real>(how_many);

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

        amrex::GPU::synchronize();

        return {det_x, det_y, det_z}

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
    auto omega_range std::vector<amrex::Real>(2);
    pp_radiation.getarr("omega_range", omega_range);
    std::copy(omega_range.begin(), omega_range.end(), m_omega_range.begin());
    pp_radiation.get("omega_points", m_omega_points);

    //Angle theta AND phi
    auto theta_range std::vector<amrex::Real>(2);
    pp_radiation.get("angle_aperture", theta_range);
    std::copy(theta_range.begin(), theta_range.end(), m_theta_range.begin());

    //Detector parameters
    auto det_pts std::vector<int>(2);
    pp_radiation.getarr("detector_number_points", det_pts);
    std::copy(det_pts.begin(), det_pts.end(), m_det_pts.begin());

    auto det_direction std::vector<amrex::Real>(3);
    pp_radiation.getarr("detector_direction", det_direction);
    std::copy(det_direction.begin(), det_direction.end(), m_det_direction.begin());

    auto det_orientation std::vector<amrex::Real>(3);
    pp_radiation.getarr("detector_orientation", det_orientation);
    std::copy(det_orientation.begin(), det_orientation.end(), m_det_orientation.begin());

    pp_radiation.get("detector_distance", m_det_distance);

    [det_pos_x, det_pos_y, det_pos_z] = compute_detector_positions(
        center, det_direction, m_det_distance,
        m_det_orientation, m_det_pts, m_theta_range);

    constexpr auto ncomp = 3;
    m_radiation_data = amrex::Gpu::DeviceVector<ablastr::math::Complex>(m_det_pts[0]*m_det_pts[1]*m_omega_points*ncomp);

    m_omegas = amrex::Gpu::DeviceVector<amrex::Real>(m_omega_points);
    ablastr::math::LinSpaceFill(m_omega_range[0], m_omega_range[1], m_omega_points, omega_calc.being());
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
                    auto& attribs = pti.GetAttribs();
                    auto& p_w = attribs[PIdx::w];
                    auto& p_ux = attribs[PIdx::ux];
                    auto& p_uy = attribs[PIdx::uy];
                    auto& p_uz = attribs[PIdx::uz];

                    auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                    auto& part=pc->GetParticles(lev)[index];
                    auto& soa = part.GetStructOfArrays();

                    auto* p_ux_old = soa.GetRealData(pc->GetRealCompIndex("old_u_x")).data();
                    auto* p_uy_old = soa.GetRealData(pc->GetRealCompIndex("old_u_y")).data();
                    auto* p_uz_old = soa.GetRealData(pc->GetRealCompIndex("old_u_z")).data();

                    auto GetPosition = GetParticlePosition<PIdx>(pti);
                    amrex::ParticleReal const q = pc->getCharge();

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip){
                        amrex::ParticleReal xp, yp, zp;
                        GetPosition.AsStored(ip,xp, yp, zp);


                    amrex::ParticleReal u2 = p_ux[ip]*p_ux[ip]+)p_uy[ip]*p_uy[ip]+p_uz[ip]*p_uz[ip];
                    amrex::ParticleReal gamma = sdt::sqrt(1.0_rt + u2);

                    //Calculation of ei(omegat-n.r)
                    n.resize(3);
                    //omega effective calculé
                    for(int i_om=0; i_om<m_omega_points;i_om++){
                        for(int i_x=0; i_x<m_det_pts[0]; i_x++){
                            for(int i_y=0; i_y<m_det_pts[1]; i_y++){
                                for(int idimo; idimo<3; idimo++){
                                    if(m_det_direction[idimo]==1){
                                        n[idimo]=(Part_pos[idimo]-m_det_distance)/m_det_distance;

                                    }
                                    else{
                                        n[idimo]=(Part_pos[idimo]-det_pos[(idimo+1)%2][i_y])/m_det_distance;
                                        n[idimo]=(Part_pos[idimo]-det_pos[(idimo+2)%2][i_x])/m_det_distance;
                                        }
                                }
                                    //Calculation of 1_beta.n, n corresponds to m_det_direction, the direction of the normal
                                    amrex::ParticleReal un_betan = 1-(p_ux_unit*n[0]+p_uy_unit*n[1]+p_uz_unit*n[2]);

                                    //Calculation of betapoint
                                    amrex::ParticleReal betapointx=(p_ux_unit-p_ux_old_unit)/dt/ablastr::constant::SI::c;
                                    amrex::ParticleReal betapointy=(p_uy_unit-p_uy_old_unit)/dt/ablastr::constant::SI::c;
                                    amrex::ParticleReal betapointz=(p_uz_unit-p_uz_old_unit)/dt/ablastr::constant::SI::c;

                                    //Calculation of nxbeta
                                    amrex::ParticleReal ncrossBetapointx=(n[1]-p_uy_unit/gamma)*betapointz - (n[2]-p_uz_unit/gamma)*betapointy;
                                    amrex::ParticleReal ncrossBetapointy=(n[2]-p_uz_unit/gamma)*betapointx - (n[0]-p_ux_unit/gamma)*betapointz;
                                    amrex::ParticleReal ncrossBetapointz=(n[0]-p_ux_unit/gamma)*betapointy - (n[1]-p_uy_unit/gamma)*betapointx;

                                    //Calculation of nxnxbeta
                                    amrex::ParticleReal ncrossncrossBetapointx=n[2]*ncrossBetapointz-n[2]*ncrossBetapointy;
                                    amrex::ParticleReal ncrossncrossBetapointy=n[0]*ncrossBetapointx-n[0]*ncrossBetapointz;
                                    amrex::ParticleReal ncrossncrossBetapointz=n[1]*ncrossBetapointy-n[1]*ncrossBetapointx;

                                     //function cospi and not cos
                                    Complex eiomega=(amrex::Math::cospi(dephas/ablastr::constant::math::pi), amrex::Math::cospi((ablastr::constant::math::pi/2-dephas)/ablastr::constant::math::pi));
                                    Complex Term_x= q*dt/(16*pow(ablastr::constant::math::pi,3)*ablastr::constant::SI::ep0*ablastr::constant::SI::c)*ncrossncrossBetapointx/std::pow(un_betan,2)*eiomega;
                                    Complex Term_y= q*dt/(16*pow(ablastr::constant::math::pi,3)*ablastr::constant::SI::ep0*ablastr::constant::SI::c)*ncrossncrossBetapointy/std::pow(un_betan,2)*eiomega;
                                    Complex Term_z= q*dt/(16*pow(ablastr::constant::math::pi,3)*ablastr::constant::SI::ep0*ablastr::constant::SI::c)*ncrossncrossBetapointz/std::pow(un_betan,2)*eiomega;

                                    //Add the contributions
                                    m_radiation_data[i_x,i_y,i_om,0]+=Term_x.m_real;
                                    m_radiation_data[i_x,i_y,i_om,1]+=Term_x.m_imag;
                                    m_radiation_data[i_x,i_y,i_om,2]+=Term_y.m_real;
                                    m_radiation_data[i_x,i_y,i_om,3]+=Term_y.m_imag;
                                    m_radiation_data[i_x,i_y,i_om,4]+=Term_z.m_real;
                                    m_radiation_data[i_x,i_y,i_om,5]+=Term_z.m_imag;

                         }
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
        m_radiation_calculation.begin(), m_radiation_calculation.end(), radiation_data_cpu.begin());
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

void RadiationHandler::Integral_overtime()
{
    for(int i_om=0; i_om<m_omega_points;i_om++){
                        for(int i_x=0; i_x<m_det_pts[0]; i_x++){
                            for(int i_y=0; i_y<m_det_pts[1]; i_y++){
                                 m_radiation_calculation[i_x,i_y,i_om]=std::pow(m_radiation_data[i_x,i_y,i_om,0],2)+std::pow(m_radiation_data[i_x,i_y,i_om,1],2)+std::pow(m_radiation_data[i_x,i_y,i_om,2],2)+std::pow(m_radiation_data[i_x,i_y,i_om,3],2)+std::pow(m_radiation_data[i_x,i_y,i_om,4],2)+std::pow(m_radiation_data[i_x,i_y,i_om,5],2);
                            }
                        }
    }
}
