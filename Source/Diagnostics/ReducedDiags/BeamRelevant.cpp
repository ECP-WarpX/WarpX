/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BeamRelevant.H"

#include "Diagnostics/ReducedDiags/ReducedDiags.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/IntervalsParser.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#include <AMReX_REAL.H>
#include <AMReX_Tuple.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <vector>

using namespace amrex;

// constructor
BeamRelevant::BeamRelevant (std::string rd_name)
: ReducedDiags{rd_name}
{
    // read beam name
    ParmParse pp_rd_name(rd_name);
    pp_rd_name.get("species",m_beam_name);

    // resize data array
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
    //  0, 1, 2: mean x,y,z
    //  3, 4, 5: mean px,py,pz
    //        6: gamma
    //  7, 8, 9: rms x,y,z
    // 10,11,12: rms px,py,pz
    //       13: rms gamma
    // 14,15,16: emittance x,y,z
    //       17: charge
    m_data.resize(18, 0.0_rt);
#elif (defined WARPX_DIM_XZ)
    //     0, 1: mean x,z
    //  2, 3, 4: mean px,py,pz
    //        5: gamma
    //     6, 7: rms x,z
    //  8, 9,10: rms px,py,pz
    //       11: rms gamma
    //    12,13: emittance x,z
    //       14: charge
    m_data.resize(15, 0.0_rt);
#endif

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension, std::ofstream::out};
            // write header row
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
            int c = 0;
            ofs << "#";
            ofs << "[" << c++ << "]step()";           ofs << m_sep;
            ofs << "[" << c++ << "]time(s)";          ofs << m_sep;
            ofs << "[" << c++ << "]x_mean(m)";        ofs << m_sep;
            ofs << "[" << c++ << "]y_mean(m)";        ofs << m_sep;
            ofs << "[" << c++ << "]z_mean(m)";        ofs << m_sep;
            ofs << "[" << c++ << "]px_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[" << c++ << "]py_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[" << c++ << "]pz_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[" << c++ << "]gamma_mean()";     ofs << m_sep;
            ofs << "[" << c++ << "]x_rms(m)";         ofs << m_sep;
            ofs << "[" << c++ << "]y_rms(m)";         ofs << m_sep;
            ofs << "[" << c++ << "]z_rms(m)";         ofs << m_sep;
            ofs << "[" << c++ << "]px_rms(kg*m/s)";   ofs << m_sep;
            ofs << "[" << c++ << "]py_rms(kg*m/s)";   ofs << m_sep;
            ofs << "[" << c++ << "]pz_rms(kg*m/s)";   ofs << m_sep;
            ofs << "[" << c++ << "]gamma_rms()";      ofs << m_sep;
            ofs << "[" << c++ << "]emittance_x(m)";   ofs << m_sep;
            ofs << "[" << c++ << "]emittance_y(m)";   ofs << m_sep;
            ofs << "[" << c++ << "]emittance_z(m)";   ofs << m_sep;
            ofs << "[" << c++ << "]charge(C)";        ofs << std::endl;
#elif (defined WARPX_DIM_XZ)
            int c = 0;
            ofs << "#";
            ofs << "[" << c++ << "]step()";           ofs << m_sep;
            ofs << "[" << c++ << "]time(s)";          ofs << m_sep;
            ofs << "[" << c++ << "]x_mean(m)";        ofs << m_sep;
            ofs << "[" << c++ << "]z_mean(m)";        ofs << m_sep;
            ofs << "[" << c++ << "]px_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[" << c++ << "]py_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[" << c++ << "]pz_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[" << c++ << "]gamma_mean()";     ofs << m_sep;
            ofs << "[" << c++ << "]x_rms(m)";         ofs << m_sep;
            ofs << "[" << c++ << "]z_rms(m)";         ofs << m_sep;
            ofs << "[" << c++ << "]px_rms(kg*m/s)";   ofs << m_sep;
            ofs << "[" << c++ << "]py_rms(kg*m/s)";   ofs << m_sep;
            ofs << "[" << c++ << "]pz_rms(kg*m/s)";   ofs << m_sep;
            ofs << "[" << c++ << "]gamma_rms()";      ofs << m_sep;
            ofs << "[" << c++ << "]emittance_x(m)";   ofs << m_sep;
            ofs << "[" << c++ << "]emittance_z(m)";   ofs << m_sep;
            ofs << "[" << c++ << "]charge(C)";        ofs << std::endl;
#endif
            // close file
            ofs.close();
        }
    }
}
// end constructor

// function that compute beam relevant quantities
void BeamRelevant::ComputeDiags (int step)
{
    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get MultiParticleContainer class object
    const auto & mypc = WarpX::GetInstance().GetPartContainer();

    // get number of species
    int const nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    auto const species_names = mypc.GetSpeciesNames();

    // inverse of speed of light squared
    Real constexpr inv_c2 = 1.0_rt / (PhysConst::c * PhysConst::c);

    // If 2D-XZ, p.pos(1) is z, rather than p.pos(2).
#if (defined WARPX_DIM_3D)
    int const index_z = 2;
#elif (defined WARPX_DIM_XZ || defined WARPX_DIM_RZ)
    int const index_z = 1;
#endif

    // loop over species
    for (int i_s = 0; i_s < nSpecies; ++i_s)
    {
        // only beam species does
        if (species_names[i_s] != m_beam_name) { continue; }

        // get WarpXParticleContainer class object
        auto const & myspc = mypc.GetParticleContainer(i_s);

        // get mass and charge
        ParticleReal const m = myspc.getMass();
        ParticleReal const q = myspc.getCharge();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        amrex::ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,
        ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_ops;
        auto r = amrex::ParticleReduce<amrex::ReduceData<ParticleReal,ParticleReal,
        ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple
            <ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,
            ParticleReal,ParticleReal,ParticleReal>
            {
                const ParticleReal p_ux = p.rdata(PIdx::ux);
                const ParticleReal p_uy = p.rdata(PIdx::uy);
                const ParticleReal p_uz = p.rdata(PIdx::uz);
                const ParticleReal p_us = p_ux*p_ux + p_uy*p_uy + p_uz*p_uz;
                const ParticleReal p_pos0 = p.pos(0);
                const ParticleReal p_w = p.rdata(PIdx::w);

#if (defined WARPX_DIM_RZ)
                const ParticleReal p_theta = p.rdata(PIdx::theta);
                const ParticleReal p_x_mean = p_pos0*std::cos(p_theta)*p_w;
                const ParticleReal p_y_mean = p_pos0*std::sin(p_theta)*p_w;
#else
                const ParticleReal p_pos1 = p.pos(1);
                const ParticleReal p_x_mean = p_pos0*p_w;
                const ParticleReal p_y_mean = p_pos1*p_w;
#endif
                const ParticleReal p_z_mean = p.pos(index_z)*p_w;

                const ParticleReal p_ux_mean = p_ux*p_w;
                const ParticleReal p_uy_mean = p_uy*p_w;
                const ParticleReal p_uz_mean = p_uz*p_w;
                const ParticleReal p_gm_mean = std::sqrt(1.0_rt+p_us*inv_c2)*p_w;

                return {p_w,
                        p_x_mean, p_y_mean, p_z_mean,
                        p_ux_mean, p_uy_mean, p_uz_mean,
                        p_gm_mean};
            },
            reduce_ops);

        std::vector<ParticleReal> values_per_rank_1st = {
            amrex::get<0>(r), // w
            amrex::get<1>(r), // x_mean
            amrex::get<2>(r), // y_mean
            amrex::get<3>(r), // z_mean
            amrex::get<4>(r), // ux_mean
            amrex::get<5>(r), // uy_mean
            amrex::get<6>(r), // uz_mean
            amrex::get<7>(r), // gm_mean
        };

        // reduced sum over mpi ranks (allreduce)
        amrex::ParallelAllReduce::Sum
        ( values_per_rank_1st.data(), values_per_rank_1st.size(), ParallelDescriptor::Communicator());

        ParticleReal w_sum   = values_per_rank_1st.at(0);
        ParticleReal x_mean  = values_per_rank_1st.at(1) /= w_sum;
        ParticleReal y_mean  = values_per_rank_1st.at(2) /= w_sum;
        ParticleReal z_mean  = values_per_rank_1st.at(3) /= w_sum;
        ParticleReal ux_mean = values_per_rank_1st.at(4) /= w_sum;
        ParticleReal uy_mean = values_per_rank_1st.at(5) /= w_sum;
        ParticleReal uz_mean = values_per_rank_1st.at(6) /= w_sum;
        ParticleReal gm_mean = values_per_rank_1st.at(7) /= w_sum;

        if (w_sum < std::numeric_limits<Real>::min() )
        {
            for (int i = 0; i < static_cast<int>(m_data.size()); ++i){
                m_data[i] = 0.0_rt;
            }
            return;
        }

        amrex::ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,
        ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_ops2;

        auto r2 = amrex::ParticleReduce<amrex::ReduceData<ParticleReal,ParticleReal,
        ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,
        ParticleReal,ParticleReal,ParticleReal>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple
            <ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal,
            ParticleReal,ParticleReal,ParticleReal,ParticleReal,ParticleReal>
            {
                const ParticleReal p_ux = p.rdata(PIdx::ux);
                const ParticleReal p_uy = p.rdata(PIdx::uy);
                const ParticleReal p_uz = p.rdata(PIdx::uz);
                const ParticleReal p_us = p_ux*p_ux + p_uy*p_uy + p_uz*p_uz;
                const ParticleReal p_gm = std::sqrt(1.0_rt+p_us*inv_c2);
                const ParticleReal p_pos0 = p.pos(0);
                const ParticleReal p_w = p.rdata(PIdx::w);

#if (defined WARPX_DIM_RZ)
                const ParticleReal p_theta = p.rdata(PIdx::theta);
                const ParticleReal p_x = p_pos0*std::cos(p_theta);
                const ParticleReal p_y = p_pos0*std::sin(p_theta);
#else
                const ParticleReal p_pos1 = p.pos(1);
                const ParticleReal p_x = p_pos0;
                const ParticleReal p_y = p_pos1;
#endif
                const ParticleReal p_z = p.pos(index_z);

                const ParticleReal p_x_ms = (p_x-x_mean)*(p_x-x_mean)*p_w;
                const ParticleReal p_y_ms = (p_y-y_mean)*(p_y-y_mean)*p_w;
                const ParticleReal p_z_ms = (p_z-z_mean)*(p_z-z_mean)*p_w;

                const ParticleReal p_ux_ms = (p_ux-ux_mean)*(p_ux-ux_mean)*p_w;
                const ParticleReal p_uy_ms = (p_uy-uy_mean)*(p_uy-uy_mean)*p_w;
                const ParticleReal p_uz_ms = (p_uz-uz_mean)*(p_uz-uz_mean)*p_w;
                const ParticleReal p_gm_ms = (p_gm-gm_mean)*(p_gm-gm_mean)*p_w;

                const ParticleReal p_xux = (p_x-x_mean)*(p_ux-ux_mean)*p_w;
                const ParticleReal p_yuy = (p_y-y_mean)*(p_uy-uy_mean)*p_w;
                const ParticleReal p_zuz = (p_z-z_mean)*(p_uz-uz_mean)*p_w;

                const ParticleReal p_charge = q*p_w;

                return {p_x_ms, p_y_ms, p_z_ms,
                        p_ux_ms, p_uy_ms, p_uz_ms,
                        p_gm_ms,
                        p_xux, p_yuy, p_zuz,
                        p_charge};
            },
            reduce_ops2);

        std::vector<ParticleReal> values_per_rank_2nd = {
            amrex::get<0>(r2), // x_ms
            amrex::get<1>(r2), // y_ms
            amrex::get<2>(r2), // z_ms
            amrex::get<3>(r2), // ux_ms
            amrex::get<4>(r2), // uy_ms
            amrex::get<5>(r2), // uz_ms
            amrex::get<6>(r2), // gm_ms
            amrex::get<7>(r2), // xux
            amrex::get<8>(r2), // yuy
            amrex::get<9>(r2), // zuz
            amrex::get<10>(r2) // charge
        };

        // reduced sum over mpi ranks (reduce to IO rank)
        ParallelDescriptor::ReduceRealSum
        ( values_per_rank_2nd.data(), values_per_rank_2nd.size(), ParallelDescriptor::IOProcessorNumber());

        ParticleReal x_ms   = values_per_rank_2nd.at(0) /= w_sum;
        ParticleReal y_ms   = values_per_rank_2nd.at(1) /= w_sum;
        ParticleReal z_ms   = values_per_rank_2nd.at(2) /= w_sum;
        ParticleReal ux_ms  = values_per_rank_2nd.at(3) /= w_sum;
        ParticleReal uy_ms  = values_per_rank_2nd.at(4) /= w_sum;
        ParticleReal uz_ms  = values_per_rank_2nd.at(5) /= w_sum;
        ParticleReal gm_ms  = values_per_rank_2nd.at(6) /= w_sum;
        ParticleReal xux    = values_per_rank_2nd.at(7) /= w_sum;
        ParticleReal yuy    = values_per_rank_2nd.at(8) /= w_sum;
        ParticleReal zuz    = values_per_rank_2nd.at(9) /= w_sum;
        ParticleReal charge = values_per_rank_2nd.at(10);

        // save data
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
        m_data[0]  = x_mean;
        m_data[1]  = y_mean;
        m_data[2]  = z_mean;
        m_data[3]  = ux_mean * m;
        m_data[4]  = uy_mean * m;
        m_data[5]  = uz_mean * m;
        m_data[6]  = gm_mean;
        m_data[7]  = std::sqrt(x_ms);
        m_data[8]  = std::sqrt(y_ms);
        m_data[9]  = std::sqrt(z_ms);
        m_data[10] = std::sqrt(ux_ms) * m;
        m_data[11] = std::sqrt(uy_ms) * m;
        m_data[12] = std::sqrt(uz_ms) * m;
        m_data[13] = std::sqrt(gm_ms);
        m_data[14] = std::sqrt(x_ms*ux_ms-xux*xux) / PhysConst::c;
        m_data[15] = std::sqrt(y_ms*uy_ms-yuy*yuy) / PhysConst::c;
        m_data[16] = std::sqrt(z_ms*uz_ms-zuz*zuz) / PhysConst::c;
        m_data[17] = charge;
#elif (defined WARPX_DIM_XZ)
        m_data[0]  = x_mean;
        m_data[1]  = z_mean;
        m_data[2]  = ux_mean * m;
        m_data[3]  = uy_mean * m;
        m_data[4]  = uz_mean * m;
        m_data[5]  = gm_mean;
        m_data[6]  = std::sqrt(x_ms);
        m_data[7]  = std::sqrt(z_ms);
        m_data[8]  = std::sqrt(ux_ms) * m;
        m_data[9]  = std::sqrt(uy_ms) * m;
        m_data[10] = std::sqrt(uz_ms) * m;
        m_data[11] = std::sqrt(gm_ms);
        m_data[12] = std::sqrt(x_ms*ux_ms-xux*xux) / PhysConst::c;
        m_data[13] = std::sqrt(z_ms*uz_ms-zuz*zuz) / PhysConst::c;
        m_data[14] = charge;
        amrex::ignore_unused(y_mean, y_ms, yuy);
#endif
    }
    // end loop over species
}
// end void BeamRelevant::ComputeDiags
