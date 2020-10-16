/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BeamRelevant.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>

#include <iostream>
#include <cmath>
#include <limits>


using namespace amrex;

// constructor
BeamRelevant::BeamRelevant (std::string rd_name)
: ReducedDiags{rd_name}
{

    // read beam name
    ParmParse pp(rd_name);
    pp.get("species",m_beam_name);

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
    m_data.resize(18,0.0);
#elif (defined WARPX_DIM_XZ)
    //     0, 1: mean x,z
    //  2, 3, 4: mean px,py,pz
    //        5: gamma
    //     6, 7: rms x,z
    //  8, 9,10: rms px,py,pz
    //       11: rms gamma
    //    12,13: emittance x,z
    //       14: charge
    m_data.resize(15,0.0);
#endif

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app};
            // write header row
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
            ofs << "#";
            ofs << "[1]step()";           ofs << m_sep;
            ofs << "[2]time(s)";          ofs << m_sep;
            ofs << "[3]x_mean(m)";        ofs << m_sep;
            ofs << "[4]y_mean(m)";        ofs << m_sep;
            ofs << "[5]z_mean(m)";        ofs << m_sep;
            ofs << "[6]px_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[7]py_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[8]pz_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[9]gamma_mean()";     ofs << m_sep;
            ofs << "[10]x_rms(m)";        ofs << m_sep;
            ofs << "[11]y_rms(m)";        ofs << m_sep;
            ofs << "[12]z_rms(m)";        ofs << m_sep;
            ofs << "[13]px_rms(kg*m/s)";  ofs << m_sep;
            ofs << "[14]py_rms(kg*m/s)";  ofs << m_sep;
            ofs << "[15]pz_rms(kg*m/s)";  ofs << m_sep;
            ofs << "[16]gamma_rms()";     ofs << m_sep;
            ofs << "[17]emittance_x(m)";  ofs << m_sep;
            ofs << "[18]emittance_y(m)";  ofs << m_sep;
            ofs << "[19]emittance_z(m)";  ofs << m_sep;
            ofs << "[20]charge(C)";       ofs << std::endl;
#elif (defined WARPX_DIM_XZ)
            ofs << "#";
            ofs << "[1]step()";           ofs << m_sep;
            ofs << "[2]time(s)";          ofs << m_sep;
            ofs << "[3]x_mean(m)";        ofs << m_sep;
            ofs << "[4]z_mean(m)";        ofs << m_sep;
            ofs << "[5]px_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[6]py_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[7]pz_mean(kg*m/s)";  ofs << m_sep;
            ofs << "[8]gamma_mean()";     ofs << m_sep;
            ofs << "[9]x_rms(m)";         ofs << m_sep;
            ofs << "[10]z_rms(m)";        ofs << m_sep;
            ofs << "[11]px_rms(kg*m/s)";  ofs << m_sep;
            ofs << "[12]py_rms(kg*m/s)";  ofs << m_sep;
            ofs << "[13]pz_rms(kg*m/s)";  ofs << m_sep;
            ofs << "[14]gamma_rms()";     ofs << m_sep;
            ofs << "[15]emittance_x(m)";  ofs << m_sep;
            ofs << "[16]emittance_z(m)";  ofs << m_sep;
            ofs << "[17]charge(C)";       ofs << std::endl;
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
    Real constexpr inv_c2 = 1.0 / (PhysConst::c * PhysConst::c);

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

        // get mass and charge (Real), FIXME actually all here are ParticleReal
        Real const m = myspc.getMass();
        Real const q = myspc.getCharge();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        // weight sum
        Real w_sum = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::w); });

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum(w_sum);

        if (w_sum < std::numeric_limits<Real>::min() )
        {
            for (int i = 0; i < static_cast<int>(m_data.size()); ++i){
                m_data[i] = 0.0;
            }
            return;
        }

#if (defined WARPX_DIM_RZ)
        // x mean
        Real x_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0)*std::cos(p.rdata(PIdx::theta)) * p.rdata(PIdx::w); });
#else
        // x mean
        Real x_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0) * p.rdata(PIdx::w); });
#endif

#if (defined WARPX_DIM_3D)
        // y mean
        Real y_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(1) * p.rdata(PIdx::w); });
#elif (defined WARPX_DIM_RZ)
        // y mean
        Real y_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0)*std::sin(p.rdata(PIdx::theta)) * p.rdata(PIdx::w); });
#endif

        // z mean
        Real z_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(index_z) * p.rdata(PIdx::w); });

        // ux mean
        Real ux_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::ux) * p.rdata(PIdx::w); });

        // uy mean
        Real uy_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uy) * p.rdata(PIdx::w); });

        // uz mean
        Real uz_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uz) * p.rdata(PIdx::w); });

        // gamma mean
        Real gm_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real ux = p.rdata(PIdx::ux);
            Real uy = p.rdata(PIdx::uy);
            Real uz = p.rdata(PIdx::uz);
            Real us = ux*ux + uy*uy + uz*uz;
            return std::sqrt(1.0 + us*inv_c2) * p.rdata(PIdx::w);
        });

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum(x_mean);  x_mean  /= w_sum;
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
        ParallelDescriptor::ReduceRealSum(y_mean);  y_mean  /= w_sum;
#endif
        ParallelDescriptor::ReduceRealSum(z_mean);  z_mean  /= w_sum;
        ParallelDescriptor::ReduceRealSum(ux_mean); ux_mean /= w_sum;
        ParallelDescriptor::ReduceRealSum(uy_mean); uy_mean /= w_sum;
        ParallelDescriptor::ReduceRealSum(uz_mean); uz_mean /= w_sum;
        ParallelDescriptor::ReduceRealSum(gm_mean); gm_mean /= w_sum;

#if (defined WARPX_DIM_RZ)
        // x mean square
        Real x_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const x = p.pos(0)*std::cos(p.rdata(PIdx::theta));
            Real const a = (x-x_mean) * (x-x_mean);
            return a * p.rdata(PIdx::w);
        });
#else
        // x mean square
        Real x_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.pos(0)-x_mean) * (p.pos(0)-x_mean);
            return a * p.rdata(PIdx::w);
        });
#endif

#if (defined WARPX_DIM_3D)
        // y mean square
        Real y_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.pos(1)-y_mean) * (p.pos(1)-y_mean);
            return a * p.rdata(PIdx::w);
        });
#elif (defined WARPX_DIM_RZ)
        // y mean square
        Real y_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const y = p.pos(0)*std::sin(p.rdata(PIdx::theta));
            Real const a = (y-y_mean) * (y-y_mean);
            return a * p.rdata(PIdx::w);
        });
#endif

        // z mean square
        Real z_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.pos(index_z)-z_mean) * (p.pos(index_z)-z_mean);
            return a * p.rdata(PIdx::w);
        });

        // ux mean square
        Real ux_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.rdata(PIdx::ux)-ux_mean) *
                           (p.rdata(PIdx::ux)-ux_mean);
            return a * p.rdata(PIdx::w);
        });

        // uy mean square
        Real uy_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.rdata(PIdx::uy)-uy_mean) *
                           (p.rdata(PIdx::uy)-uy_mean);
            return a * p.rdata(PIdx::w);
        });

        // uz mean square
        Real uz_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.rdata(PIdx::uz)-uz_mean) *
                           (p.rdata(PIdx::uz)-uz_mean);
            return a * p.rdata(PIdx::w);
        });

        // gamma mean square
        Real gm_ms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const ux = p.rdata(PIdx::ux);
            Real const uy = p.rdata(PIdx::uy);
            Real const uz = p.rdata(PIdx::uz);
            Real const us = ux*ux + uy*uy + uz*uz;
            Real const gm = std::sqrt(1.0 + us*inv_c2);
            Real const a  = (gm - gm_mean) * (gm - gm_mean);
            return a * p.rdata(PIdx::w);
        });

#if (defined WARPX_DIM_RZ)
        // x times ux
        Real xux = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const x = p.pos(0)*std::cos(p.rdata(PIdx::theta));
            Real const a = (x-x_mean) * (p.rdata(PIdx::ux)-ux_mean);
            return a * p.rdata(PIdx::w);
        });
#else
        // x times ux
        Real xux = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.pos(0)-x_mean) * (p.rdata(PIdx::ux)-ux_mean);
            return a * p.rdata(PIdx::w);
        });
#endif

#if (defined WARPX_DIM_3D)
        // y times uy
        Real yuy = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.pos(1)-y_mean) * (p.rdata(PIdx::uy)-uy_mean);
            return a * p.rdata(PIdx::w);
        });
#elif (defined WARPX_DIM_RZ)
        // y times uy
        Real yuy = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const y = p.pos(0)*std::sin(p.rdata(PIdx::theta));
            Real const a = (y-y_mean) * (p.rdata(PIdx::uy)-uy_mean);
            return a * p.rdata(PIdx::w);
        });
#endif

        // z times uz
        Real zuz = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const a = (p.pos(index_z)-z_mean) * (p.rdata(PIdx::uz)-uz_mean);
            return a * p.rdata(PIdx::w);
        });

        // charge
        Real charge = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real const w = p.rdata(PIdx::w);
            return q * w;
        });

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum
            ( x_ms, ParallelDescriptor::IOProcessorNumber());
        x_ms /= w_sum;
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
        ParallelDescriptor::ReduceRealSum
            ( y_ms, ParallelDescriptor::IOProcessorNumber());
        y_ms /= w_sum;
#endif
        ParallelDescriptor::ReduceRealSum
            ( z_ms, ParallelDescriptor::IOProcessorNumber());
        z_ms /= w_sum;
        ParallelDescriptor::ReduceRealSum
            (ux_ms, ParallelDescriptor::IOProcessorNumber());
        ux_ms /= w_sum;
        ParallelDescriptor::ReduceRealSum
            (uy_ms, ParallelDescriptor::IOProcessorNumber());
        uy_ms /= w_sum;
        ParallelDescriptor::ReduceRealSum
            (uz_ms, ParallelDescriptor::IOProcessorNumber());
        uz_ms /= w_sum;
        ParallelDescriptor::ReduceRealSum
            (gm_ms, ParallelDescriptor::IOProcessorNumber());
        gm_ms /= w_sum;
        ParallelDescriptor::ReduceRealSum
            (   xux, ParallelDescriptor::IOProcessorNumber());
        xux /= w_sum;
#if (defined WARPX_DIM_3D || defined WARPX_DIM_RZ)
        ParallelDescriptor::ReduceRealSum
            (   yuy, ParallelDescriptor::IOProcessorNumber());
        yuy /= w_sum;
#endif
        ParallelDescriptor::ReduceRealSum
            (   zuz, ParallelDescriptor::IOProcessorNumber());
        zuz /= w_sum;
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
