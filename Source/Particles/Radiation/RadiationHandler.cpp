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
        pp_radiations.getarr("detector_direction", m_det_direction,0,3);
        pp_radiations.query("detector_distance", m_det_distance);
        add_detector();
         /* Initialize the Fab with the field in function of angle and frequency */
         //int numcomps = 4;
        //BaseFab<Complex> fab(bx,numcomps);
        // Box of angle and frequency
        //const amrex::Box detect_box({0,0,0}, {m_det_pts[0], m_det_pts[1], m_omega_points});
        ncomp = 6; //Real and complex part
        //amrex::FArrayBox fab_detect(detect_box, ncomp);
        //fab_detect.setval(0,0)
        m_radiation_data = amrex::Gpu::DeviceVector<amrex::Real>(m_det_pts[0]*m_det_pts[1]*m_omega_points*ncomp);
        m_radiation_calculation = amrex::Gpu::DeviceVector<amrex::Real>(m_det_pts[0]*m_det_pts[1]*m_omega_points);

        det_pos.resize(2);
        det_pos[0].resize(m_det_pts[0]);
        det_pos[1].resize(m_det_pts[1]);
        omega_calc.resize(m_omega_points);

        //Initialization :
        for(int i=0; i<m_det_pts[0];i++){
            det_pos[0][i]=det_bornes[0][0]+i*m_d_d[0];
                    }
        for(int i=0; i<m_det_pts[1];i++){
            det_pos[1][i]=det_bornes[1][0]+i*m_d_d[1];
                    }
        for(int i=0; i<m_omega_points;i++){
            omega_calc[i]=m_omega_range[0]+i*m_d_omega;
                    }

    }
}

void RadiationHandler::add_radiation_contribution
    (const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc, amrex::Real current_time){
        const auto level0=0;        
#ifdef AMREX_USE_OMP
#pragma omp parallel 
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
                    amrex::ParticleReal const q = pc->getCharge();

                amrex::ParallelFor(np,
                 [=] AMREX_GPU_DEVICE(int ip)
                 {   amrex::ParticleReal xp, yp, zp;

                    GetPosition.AsStored(ip,xp, yp, zp);
                    amrex::GpuArray<amrex::Real, 3> Part_pos{xp,yp,zp};

                    amrex::ParticleReal p_ux_old_unit=p_ux_old[ip];
                    amrex::ParticleReal p_uy_old_unit=p_uy_old[ip];
                    amrex::ParticleReal p_uz_old_unit=p_uz_old[ip];

                    amrex::ParticleReal p_ux_unit=p_ux[ip];
                    amrex::ParticleReal p_uy_unit=p_uy[ip];
                    amrex::ParticleReal p_uz_unit=p_uz[ip];
                    
                    amrex::ParticleReal p_u = std::sqrt(std::pow(p_ux_unit,2)+std::pow(p_uy_unit,2)+std::pow(p_uz_unit,2));

                    //Calculation of gamma 
                    amrex::ParticleReal gamma = 1/(1-std::pow(p_u,2)/std::pow(ablastr::constant::SI::c,2));

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
                                        n[idimo]=(Part_pos[idimo]-det_pos[(idimo+1)%3][i_y])/m_det_distance;
                                        n[idimo]=(Part_pos[idimo]-det_pos[(idimo+2)%3][i_x])/m_det_distance;                                                   
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
                                    m_radiation_data[i_x,i_y,i_om,0]+=std::real(Term_x);
                                    m_radiation_data[i_x,i_y,i_om,1]+=std::imag(Term_x);
                                    m_radiation_data[i_x,i_y,i_om,2]+=std::real(Term_y);
                                    m_radiation_data[i_x,i_y,i_om,3]+=std::imag(Term_y);
                                    m_radiation_data[i_x,i_y,i_om,4]+=std::real(Term_z);
                                    m_radiation_data[i_x,i_y,i_om,5]+=std::imag(Term_z);

                         }
                        }
                    }
                
                });
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
void RadiationHandler::add_detector
    (){
    const auto level0=0;
    amrex::Geometry const& geom = WarpX::GetInstance().Geom(level0);
    //Calculation of angle resolution 
    m_d_theta.resize(2);
    for(int i=0; i<2; i++){
        m_d_theta[i] = 2*m_theta_range[i]/static_cast<double>(m_det_pts[i]);
    }
    //Set the resolution of the detector
    m_d_d.resize(3);
    amrex::Vector<amrex::Real> center;
    center.resize(3);
    for(int idim=0; idim<AMREX_SPACEDIM; idim++){
        center[idim] = (geom.ProbHi()[idim]+geom.ProbLo()[idim])/2;
    }
    for(int idi = 0; idi<3; ++idi){
        if(m_det_direction[idi]==1){
            m_d_d[(idi+1)%3]=2*m_det_distance*tan(m_d_theta[0]/2);
            m_d_d[(idi+2)%3]=2*m_det_distance*tan(m_d_theta[1]/2);
        }

    m_d_omega=m_omega_range[1]-m_omega_range[0]/static_cast<double>(m_omega_points);

    //Calculate the sides of the detector 
    det_bornes.resize(2);
    det_bornes[0].resize(2);
    det_bornes[1].resize(2);

    for(int idim = 0; idim<2; ++idim){
        det_bornes[idim][0]=center[idim]-std::tan(m_theta_range[idim]/2)*m_det_distance;
        det_bornes[idim][1]=center[idim]+std::tan(m_theta_range[idim]/2)*m_det_distance;
    }
    //fillWithConsecutiveReal(pos_det_x)
    //pos_det_y
    //omega_calc
    
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