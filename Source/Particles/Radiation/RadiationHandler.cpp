/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "RadiationHandler.H"

#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>

#include <vector>
using namespace warpx::utils;

RadiationHandler::RadiationHandler()
{   
    // Read in radiation input      
    const amrex::ParmParse pp_radiations("radiations");

    // Verify if there is a detector 
    pp_radiations.query("put_a_detector", m_get_a_detector);


#if defined  WARPX_DIM_RZ
    WARPX_ABORT_WITH_MESSAGE("Radiation is not supported yet with RZ.");
#endif    

    //Compute the radiation
    if(m_get_a_detector){
        //Resolution in frequency of the detector
        pp_radiations.queryarr("omega_range", m_omega_range);
        pp_radiations.query("omega_points", m_omega_points);
        //Angle theta AND phi
        pp_radiations.queryarr("angle_aperture", m_theta_range);

        //Type of detector
        pp_radiations.getarr("detector_number_points", m_det_pts,0,2);
        pp_radiations.query("detector_direction", m_det_direction);
        pp_radiations.query("detector_distance", m_det_distance);
        add_detector();
    }

    //


    /* Initialize the Fab with the field in function of angle and frequency */
    //int numcomps = 4;
    //BaseFab<Complex> fab(bx,numcomps);
}

void RadiationHandler::add_radiation_contribution
    (const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc){
        const auto level0=0;
        const auto part=pc->GetParticles(level0);
        
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
{
            for (WarpXParIter pti(*pc, level0); pti.isValid(); ++pti) {
                    long const np = pti.numParticles();
                    auto& attribs = pti.GetAttribs();
                    auto& p_w = attribs[PIdx::w];
                    auto& p_ux = attribs[PIdx::ux];
                    auto& p_uy = attribs[PIdx::uy];
                    auto& p_uz = attribs[PIdx::uz];
                    auto& attribut = attribs[PIdx::nattribs+1];
                    amrex::Print() << attribut[-1] << "\n";
                    //auto& ux = src.m_rdata[PIdx::ux][i_src];
                    //auto& uy = src.m_rdata[PIdx::uy][i_src];
                    //auto& uz = src.m_rdata[PIdx::uz][i_src];

                    //auto GetPosition = GetParticlePosition<PIdx>(pti);
                amrex::ParallelFor(np,
                 [=] AMREX_GPU_DEVICE(int ip)
                 {   //amrex::ParticleReal xp, yp, zp;
                     //GetPosition.AsStored(ip, xp, yp, zp);
                    const amrex::ParticleReal p_ux_unit=p_ux[ip];
                    const amrex::ParticleReal p_uy_unit=p_uy[ip];
                    const amrex::ParticleReal p_uz_unit=p_uz[ip];

                    amrex::ParticleReal p_u = std::sqrt(std::pow(p_ux_unit,2)+std::pow(p_uy_unit,2)+std::pow(p_uz_unit,2));

                        //const amrex::ParticleReal* p_pos0 = src_data.pos(0);

                 });

//                 const int total_partdiag_size = amrex::Scan::ExclusiveSum(np,Flag,IndexLocation);
//                 auto& ptile_dst = pc_dst.DefineAndReturnParticleTile(lev, pti.index(), pti.LocalTileIndex() );
//                 auto old_size = ptile_dst.numParticles();
//                 ptile_dst.resize(old_size + total_partdiag_size);
//                 amrex::filterParticles(ptile_dst, ptile_src, GetParticleFilter, 0, old_size, np);
//                 auto dst_data = ptile_dst.getParticleTileData();
//                 amrex::ParallelFor(np,
//                 [=] AMREX_GPU_DEVICE(int i)
//                 {
//                    if (Flag[i] == 1) GetParticleLorentzTransform(dst_data, src_data, i,
//                                                                  old_size + IndexLocation[i]);
//                 });
//                 amrex::Gpu::synchronize();
        }
    }
}

void RadiationHandler::add_detector
    (){

    // Box of angle and frequency
    const amrex::Box detect_box({0,0,0}, {m_det_pts[0], m_det_pts[1], m_omega_points});
    ncomp = 2; //Real and complex part
    amrex::FArrayBox fab_detect(detect_box, ncomp);

    //Calculation of angle resolution 
     m_d_theta.resize(2);
    for(int i=0; i<2; i++){
        m_d_theta[i] = 2*m_theta_range[i]/static_cast<double>(m_det_pts[i]);
    }
    //Set the resolution of the detector
    m_d_d.resize(2);
    for(int idim = 0; idim<3; ++idim){
        if(m_det_direction[idim]==1){
            m_d_d[(idim+1)%3]=2*m_det_distance*tan(m_d_theta[0]/2);
            m_d_d[(idim+2)%3]=2*m_det_distance*tan(m_d_theta[1]/2);
        }

    m_d_omega=m_omega_range[1]-m_omega_range[0]/static_cast<double>(m_omega_points);
    }
}