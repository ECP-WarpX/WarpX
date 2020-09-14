/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleExtrema.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>

#include <iostream>
#include <cmath>
#include <limits>


using namespace amrex;

// constructor
ParticleExtrema::ParticleExtrema (std::string rd_name)
: ReducedDiags{rd_name}
{
    // read species name
    ParmParse pp(rd_name);
    pp.get("species",m_species_name);

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    // get number of species (int)
    auto nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

    // loop over species
    for (int i_s = 0; i_s < nSpecies; ++i_s)
    {
        // only chosen species does
        if (species_names[i_s] != m_species_name) { continue; }

        // get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        if (myspc.has_breit_wheeler() || myspc.has_quantum_sync())
        {
            // resize data array for QED species
            m_data.resize(18,0.0);
        } else
        {
            // resize data array for regular species
            m_data.resize(16,0.0);
        }

        if (ParallelDescriptor::IOProcessor())
        {
            if ( m_IsNotRestart )
            {
                // open file
                std::ofstream ofs;
                ofs.open(m_path + m_rd_name + "." + m_extension,
                    std::ofstream::out | std::ofstream::app);
                // write header row
                ofs << "#";
                ofs << "[1]step()";
                ofs << "[2]time(s)";
                ofs << "[3]xmin()";
                ofs << m_sep;
                ofs << "[4]xmax()";
                ofs << m_sep;
                ofs << "[5]ymin()";
                ofs << m_sep;
                ofs << "[6]ymax()";
                ofs << m_sep;
                ofs << "[7]zmin()";
                ofs << m_sep;
                ofs << "[8]zmax()";
                ofs << m_sep;
                ofs << "[9]pxmin()";
                ofs << m_sep;
                ofs << "[10]pxmax()";
                ofs << m_sep;
                ofs << "[11]pymin()";
                ofs << m_sep;
                ofs << "[12]pymax()";
                ofs << m_sep;
                ofs << "[13]pzmin()";
                ofs << m_sep;
                ofs << "[14]pzmax()";
                ofs << m_sep;
                ofs << "[15]gmin()";
                ofs << m_sep;
                ofs << "[16]gmax()";
                ofs << m_sep;
                ofs << "[17]wmin()";
                ofs << m_sep;
                ofs << "[18]wmax()";
                if (myspc.has_breit_wheeler() || myspc.has_quantum_sync())
                {
                    ofs << m_sep;
                    ofs << "[19]chimin()";
                    ofs << m_sep;
                    ofs << "[20]chimax()";
                }
                ofs << std::endl;
                // close file
                ofs.close();
            }
        }
    }
}
// end constructor

// function that computes extrema
void ParticleExtrema::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    // get MultiParticleContainer class object
    auto & mypc = WarpX::GetInstance().GetPartContainer();

    // get number of species (int)
    auto nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

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

        // only chosen species does
        if (species_names[i_s] != m_species_name) { continue; }

        // get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        // get mass (Real)
        auto m = myspc.getMass();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        // xmin
#if (defined WARPX_DIM_RZ)
        Real xmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0)*std::cos(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMin(xmin);
#else
        Real xmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0); });
        ParallelDescriptor::ReduceRealMin(xmin);
#endif

        // xmax
#if (defined WARPX_DIM_RZ)
        Real xmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0)*std::cos(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMax(xmax);
#else
        Real xmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0); });
        ParallelDescriptor::ReduceRealMax(xmax);
#endif

        // ymin
#if (defined WARPX_DIM_RZ)
        Real ymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0)*std::sin(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMin(ymin);
#elif (defined WARPX_DIM_XZ)
        Real ymin = 0.0_rt;
#else
        Real ymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(1); });
        ParallelDescriptor::ReduceRealMin(ymin);
#endif

        // ymax
#if (defined WARPX_DIM_RZ)
        Real ymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0)*std::sin(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMax(ymax);
#elif (defined WARPX_DIM_XZ)
        Real ymax = 0.0_rt;
#else
        Real ymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(1); });
        ParallelDescriptor::ReduceRealMax(ymax);
#endif

        // zmin
        Real zmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(index_z); });
        ParallelDescriptor::ReduceRealMin(zmin);

        // zmax
        Real zmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(index_z); });
        ParallelDescriptor::ReduceRealMax(zmax);

        // uxmin
        Real uxmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::ux); });
        ParallelDescriptor::ReduceRealMin(uxmin);

        // uxmax
        Real uxmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::ux); });
        ParallelDescriptor::ReduceRealMax(uxmax);

        // uymin
        Real uymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uy); });
        ParallelDescriptor::ReduceRealMin(uymin);

        // uymax
        Real uymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uy); });
        ParallelDescriptor::ReduceRealMax(uymax);

        // uzmin
        Real uzmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uz); });
        ParallelDescriptor::ReduceRealMin(uzmin);

        // uzmax
        Real uzmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uz); });
        ParallelDescriptor::ReduceRealMax(uzmax);

        // gmin
        Real gmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real ux = p.rdata(PIdx::ux);
            Real uy = p.rdata(PIdx::uy);
            Real uz = p.rdata(PIdx::uz);
            Real us = ux*ux + uy*uy + uz*uz;
            return std::sqrt(1.0 + us*inv_c2);
        });
        ParallelDescriptor::ReduceRealMin(gmin);

        // gmax
        Real gmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real ux = p.rdata(PIdx::ux);
            Real uy = p.rdata(PIdx::uy);
            Real uz = p.rdata(PIdx::uz);
            Real us = ux*ux + uy*uy + uz*uz;
            return std::sqrt(1.0 + us*inv_c2);
        });
        ParallelDescriptor::ReduceRealMax(gmax);

        // wmin
        Real wmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::w); });
        ParallelDescriptor::ReduceRealMin(wmin);

        // wmax
        Real wmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::w); });
        ParallelDescriptor::ReduceRealMax(wmax);

        // compute chimin and chimax
        Real chimin = 0.0_rt;
        Real chimax = 0.0_rt;
        if (myspc.has_breit_wheeler() || myspc.has_quantum_sync())
        {
        }

        m_data[0]  = xmin;
        m_data[1]  = xmax;
        m_data[2]  = ymin;
        m_data[3]  = ymax;
        m_data[4]  = zmin;
        m_data[5]  = zmax;
        m_data[6]  = uxmin*m;
        m_data[7]  = uxmax*m;
        m_data[8]  = uymin*m;
        m_data[9]  = uymax*m;
        m_data[10] = uzmin*m;
        m_data[11] = uzmax*m;
        m_data[12] = gmin;
        m_data[13] = gmax;
        m_data[14] = wmin;
        m_data[15] = wmax;
        if (myspc.has_breit_wheeler() || myspc.has_quantum_sync())
        {
            m_data[16] = chimin;
            m_data[17] = chimax;
        }

    }
    // end loop over species

}
// end void ParticleEnergy::ComputeDiags
