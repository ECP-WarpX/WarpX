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

        amrex::ReduceOps<ReduceOpSum, ReduceOpSum> reduce_ops;
        auto r = amrex::ParticleReduce<amrex::ReduceData<Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real>
            {
                const ParticleReal ux = p.rdata(PIdx::ux);
                const ParticleReal uy = p.rdata(PIdx::uy);
                const ParticleReal uz = p.rdata(PIdx::uz);
                const ParticleReal us = ux*ux + uy*uy + uz*uz;
                const ParticleReal pos0 = p.pos(0);
                const ParticleReal pos1 = p.pos(1);
                const ParticleReal w = p.rdata(PIdx::w);

#if (defined WARPX_DIM_RZ)
                const ParticleReal theta = p.rdata(PIdx::theta);
                const ParticleReal x_mean = pos0*std::cos(theta)*w;
                const ParticleReal y_mean = pos0*std::sin(theta)*w;
#else
                const ParticleReal x_mean = pos0*w;
                const ParticleReal y_mean = pos1*w;
#endif
                const ParticleReal z_mean = p.pos(index_z)*w;

                const ParticleReal ux_mean = ux*w;
                const ParticleReal uy_mean = uy*w;
                const ParticleReal uz_mean = uz*w;
                const ParticleReal gm_mean = std::sqrt(1.0_rt+us*inv_c2)*w;

                return {w,x_mean,y_mean,z_mean,ux_mean,uy_mean,uz_mean,gm_mean};
            },
            reduce_ops);
            
        ParticleReal w       = amrex::get<0>(r);
        ParticleReal x_mean  = amrex::get<1>(r);
        ParticleReal y_mean  = amrex::get<2>(r);
        ParticleReal z_mean  = amrex::get<3>(r);
        ParticleReal ux_mean = amrex::get<4>(r);
        ParticleReal uy_mean = amrex::get<5>(r);
        ParticleReal uz_mean = amrex::get<6>(r);
        ParticleReal gm_mean = amrex::get<7>(r);

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum(w);
        ParallelDescriptor::ReduceRealSum(x_mean);  x_mean  /= w;
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
        ParallelDescriptor::ReduceRealSum(y_mean);  y_mean  /= w;
#endif
        ParallelDescriptor::ReduceRealSum(z_mean);  z_mean  /= w;
        ParallelDescriptor::ReduceRealSum(ux_mean); ux_mean /= w;
        ParallelDescriptor::ReduceRealSum(uy_mean); uy_mean /= w;
        ParallelDescriptor::ReduceRealSum(uz_mean); uz_mean /= w;
        ParallelDescriptor::ReduceRealSum(gm_mean); gm_mean /= w;

        if (w < std::numeric_limits<Real>::min() )
        {
            for (int i = 0; i < static_cast<int>(m_data.size()); ++i){
                m_data[i] = 0.0_rt;
            }
            return;
        }

        auto r2 = amrex::ParticleReduce<amrex::ReduceData<Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real>
            {
                const ParticleReal ux = p.rdata(PIdx::ux);
                const ParticleReal uy = p.rdata(PIdx::uy);
                const ParticleReal uz = p.rdata(PIdx::uz);
                const ParticleReal us = ux*ux + uy*uy + uz*uz;
                const ParticleReal gm = std::sqrt(1.0_rt+us*inv_c2);
                const ParticleReal pos0 = p.pos(0);
                const ParticleReal pos1 = p.pos(1);
                const ParticleReal w = p.rdata(PIdx::w);

#if (defined WARPX_DIM_RZ)
                const ParticleReal theta = p.rdata(PIdx::theta);
                const ParticleReal x = pos0*std::cos(theta);
                const ParticleReal y = pos0*std::sin(theta);
#else
                const ParticleReal x = pos0;
                const ParticleReal y = pos1;
#endif
                const ParticleReal z = p.pos(index_z);

                const ParticleReal x_ms = (x-x_mean)*(x-x_mean)*w;
                const ParticleReal y_ms = (y-y_mean)*(y-y_mean)*w;
                const ParticleReal z_ms = (z-z_mean)*(z-z_mean)*w;

                const ParticleReal ux_ms = (ux-ux_mean)*(ux-ux_mean)*w;
                const ParticleReal uy_ms = (uy-uy_mean)*(uy-uy_mean)*w;
                const ParticleReal uz_ms = (uz-uz_mean)*(uz-uz_mean)*w;
                const ParticleReal gm_ms = (gm-gm_mean)*(gm-gm_mean)*w;

                const ParticleReal xux = (x-x_mean)*(ux-ux_mean)*w;
                const ParticleReal yuy = (y-y_mean)*(uy-uy_mean)*w;
                const ParticleReal zuz = (z-z_mean)*(uz-uz_mean)*w;

                const ParticleReal charge = q*w;

                return {x_ms,y_ms,z_ms,ux_ms,uy_ms,uz_ms,gm_ms,xux,yuy,zuz,charge};
            },
            reduce_ops);

        ParticleReal x_ms   = amrex::get<0>(r2);
        ParticleReal y_ms   = amrex::get<1>(r2);
        ParticleReal z_ms   = amrex::get<2>(r2);
        ParticleReal ux_ms  = amrex::get<3>(r2);
        ParticleReal uy_ms  = amrex::get<4>(r2);
        ParticleReal uz_ms  = amrex::get<5>(r2);
        ParticleReal gm_ms  = amrex::get<6>(r2);
        ParticleReal xux    = amrex::get<7>(r2);
        ParticleReal yuy    = amrex::get<8>(r2);
        ParticleReal zuz    = amrex::get<9>(r2);
        ParticleReal charge = amrex::get<10>(r2);

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum
        ( x_ms, ParallelDescriptor::IOProcessorNumber()); x_ms /= w;
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
        ParallelDescriptor::ReduceRealSum
        ( y_ms, ParallelDescriptor::IOProcessorNumber()); y_ms /= w;
#endif
        ParallelDescriptor::ReduceRealSum
        ( z_ms, ParallelDescriptor::IOProcessorNumber()); z_ms /= w;
        ParallelDescriptor::ReduceRealSum
        (ux_ms, ParallelDescriptor::IOProcessorNumber()); ux_ms /= w;
        ParallelDescriptor::ReduceRealSum
        (uy_ms, ParallelDescriptor::IOProcessorNumber()); uy_ms /= w;
        ParallelDescriptor::ReduceRealSum
        (uz_ms, ParallelDescriptor::IOProcessorNumber()); uz_ms /= w;
        ParallelDescriptor::ReduceRealSum
        (gm_ms, ParallelDescriptor::IOProcessorNumber()); gm_ms /= w;
        ParallelDescriptor::ReduceRealSum
        (   xux, ParallelDescriptor::IOProcessorNumber()); xux /= w;
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
        ParallelDescriptor::ReduceRealSum
        (   yuy, ParallelDescriptor::IOProcessorNumber()); yuy /= w;
#endif
        ParallelDescriptor::ReduceRealSum
        (   zuz, ParallelDescriptor::IOProcessorNumber()); zuz /= w;
        ParallelDescriptor::ReduceRealSum
        ( charge, ParallelDescriptor::IOProcessorNumber());

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
#endif
    }
    // end loop over species
}
// end void BeamRelevant::ComputeDiags
