/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "RadiationHandler.H"

#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>

#include <ablastr/constant.H>


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
    (const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc, amrex::Real current_time){
        const auto level0=0;        
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

                    auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                    auto& part=pc->GetParticles(level0)[index];
                    auto& soa = part.GetStructOfArrays();

                    const auto prev_u_x_idx =pc->GetRealCompIndex("prev_u_x");
                    const auto prev_u_y_idx =pc->GetRealCompIndex("prev_u_y");
                    const auto prev_u_z_idx =pc->GetRealCompIndex("prev_u_z");

                    auto* p_ux_old = soa.GetRealData(prev_u_x_idx).data();
                    auto* p_uy_old = soa.GetRealData(prev_u_y_idx).data();
                    auto* p_uz_old = soa.GetRealData(prev_u_z_idx).data();

                    auto GetPosition = GetParticlePosition<PIdx>(pti);
                amrex::ParallelFor(np,
                 [=] AMREX_GPU_DEVICE(int ip)
                 {   amrex::ParticleReal xp, yp, zp;
                     GetPosition.AsStored(ip, xp, yp, zp);

                                     std::cout << "######### " <<  ip << " "
                            << xp << " " << yp << " " << zp << " "
                            << p_ux[ip] << " " << p_uy[ip] << " " << p_uz[ip] << " "
                            << p_ux_old[ip] << " " << p_uy_old[ip] << " " << p_uz_old[ip] << std::endl;

                    amrex::ParticleReal p_ux_old_unit=p_ux_old[ip];
                    amrex::ParticleReal p_uy_old_unit=p_uy_old[ip];
                    amrex::ParticleReal p_uz_old_unit=p_uz_old[ip];

                    amrex::ParticleReal p_ux_unit=p_ux[ip];
                    amrex::ParticleReal p_uy_unit=p_uy[ip];
                    amrex::ParticleReal p_uz_unit=p_uz[ip];
                            
                    amrex::ParticleReal const q = pc.getCharge();

                    amrex::ParticleReal p_u = std::sqrt(std::pow(p_ux_unit,2)+std::pow(p_uy_unit,2)+std::pow(p_uz_unit,2));
                    //Calculation of 1_beta.n, n corresponds to m_det_direction, the direction of the normal
                    amrex::ParticleReal un_betan = 1-(p_ux_unit*static_cast<double>(m_det_direction[0])+p_uy_unit*static_cast<double>(m_det_direction[1])+p_uz_unit*static_cast<double>(m_det_direction[2]));

                    //Calculation of gamma 
                    amrex::ParticleReal gamma = 1/(1-std::pow(p_u,2)/std::pow(ablastr::constant::SI::c,2));

                    //Calculation of betapoint
                    amrex::ParticleReal betapointx=(p_ux_unit-p_ux_old_unit)/dt/ablastr::constant::SI::c;
                    amrex::ParticleReal betapointy=(p_uy_unit-p_uy_old_unit)/dt/ablastr::constant::SI::c;
                    amrex::ParticleReal betapointz=(p_uz_unit-p_uz_old_unit)/dt/ablastr::constant::SI::c;

                    //Calculation of nxbeta
                    amrex::ParticleReal ncrossBetapointx=(static_cast<double>(m_det_direction[1])-p_uy_unit/gamma)*betapointz - (static_cast<double>(m_det_direction[2])-p_uz_unit/gamma)*betapointy;
                    amrex::ParticleReal ncrossBetapointy=(static_cast<double>(m_det_direction[2])-p_uz_unit/gamma)*betapointx - (static_cast<double>(m_det_direction[0])-p_ux_unit/gamma)*betapointz;
                    amrex::ParticleReal ncrossBetapointz=(static_cast<double>(m_det_direction[0])-p_ux_unit/gamma)*betapointy - (static_cast<double>(m_det_direction[1])-p_uy_unit/gamma)*betapointx;

                    //Calculation of nxnxbeta
                    amrex::ParticleReal ncrossncrossBetapointx=static_cast<double>(m_det_direction[2])*ncrossBetapointz-static_cast<double>(m_det_direction[2])*ncrossBetapointy;
                    amrex::ParticleReal ncrossncrossBetapointy=static_cast<double>(m_det_direction[0])*ncrossBetapointx-static_cast<double>(m_det_direction[0])*ncrossBetapointz;
                    amrex::ParticleReal ncrossncrossBetapointz=static_cast<double>(m_det_direction[1])*ncrossBetapointy-static_cast<double>(m_det_direction[1])*ncrossBetapointx;

                    //Calculation of ei(omegat-n.r)
                    
                    //omega effective calculé
                    amrex::Real omega_calc = 0
                    for(int i_om=0, i_om<m_omega_points,i_om++){
                        for(int i_x=0, i_x<m_det_pts[0], i_x++){
                            for(int i_y=0, i_y<m_det_pts[1], i_y++){
                        omega_calc = omega_calc + m_d_omega

                        for(int idim = 0; idim<3; ++idim){
                            if(m_det_direction[idim]==1){
                        amrex::Real dephas= omega_calc*(current_time-(static_cast<double>(m_det_direction[(idim+1)%3])*(xp-m_det_distance*static_cast<double>(m_det_direction[0])-pos_det_x[0])-static_cast<double>(m_det_direction[1])*(yp-m_det_distance*static_cast<double>(m_det_direction[1]))-static_cast<double>(m_det_direction[2])*(zp-m_det_distance*static_cast<double>(m_det_direction[2])))/ablastr::constant::SI::c);
                            m_d_d[(idim+1)%3]=2*m_det_distance*tan(m_d_theta[0]/2);
                            m_d_d[(idim+2)%3]=2*m_det_distance*tan(m_d_theta[1]/2);
                        }
                        amrex::Complex eiomega=(amrex::cos(dephas), amrex::sin(dephas));
                        amrex::Complex Term_x=ncrossncrossBetapointx/un_betan**2*eiomega;
                        amrex::Complex Term_y=ncrossncrossBetapointy/un_betan**2*eiomega;
                        amrex::Complex Term_z=ncrossncrossBetapointz/un_betan**2*eiomega;

                        //Calcul du terme direct de contribution
                        amrex::Real Term_int = q*dt/(16*ablastr::constant::SI::pi*ablastr::constant::SI::mu0*ablastr::constant::SI::cx)*std(std::pow(qpTerm_x,2)+std::pow(Term_y,2)+std::pow(Term_z,2));
                        fab_detect.saxpy(Real 1, fab_detect);
                        gather_and_write_radiation(const std::string& filename)    
                            }
                        }
                    }
                 }
            }
        }       
    
    });

void RadiationHandler::gather_and_write_radiation(const std::string& filename)
{
    auto radiation_data_cpu = amrex::Vector<amrex::Real>(m_I*m_J*m_W*2);
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost,
        m_radiation_data.begin(), m_radiation_data.end(), radiation_data_cpu.begin());
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
void RadiationHandler::add_detector
    (){

    // Box of angle and frequency
    const amrex::Box detect_box({0,0,0}, {m_det_pts[0], m_det_pts[1], m_omega_points});
    ncomp = 2; //Real and complex part
    amrex::FArrayBox fab_detect(detect_box, ncomp);
    fab_detect.setval(0,0)
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