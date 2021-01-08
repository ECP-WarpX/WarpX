/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleExtrema.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#if (defined WARPX_QED)
#include "Particles/ElementaryProcess/QEDInternals/QedChiFunctions.H"
#endif

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
    const auto nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    const auto species_names = mypc.GetSpeciesNames();

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
            const int num_quantities = 18;
            m_data.resize(num_quantities,0.0);
        } else
        {
            // resize data array for regular species
            const int num_quantities = 16;
            m_data.resize(num_quantities,0.0);
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
                ofs << m_sep;
                ofs << "[2]time(s)";
                ofs << m_sep;
                ofs << "[3]xmin(m)";
                ofs << m_sep;
                ofs << "[4]xmax(m)";
                ofs << m_sep;
                ofs << "[5]ymin(m)";
                ofs << m_sep;
                ofs << "[6]ymax(m)";
                ofs << m_sep;
                ofs << "[7]zmin(m)";
                ofs << m_sep;
                ofs << "[8]zmax(m)";
                ofs << m_sep;
                ofs << "[9]pxmin(kg*m/s)";
                ofs << m_sep;
                ofs << "[10]pxmax(kg*m/s)";
                ofs << m_sep;
                ofs << "[11]pymin(kg*m/s)";
                ofs << m_sep;
                ofs << "[12]pymax(kg*m/s)";
                ofs << m_sep;
                ofs << "[13]pzmin(kg*m/s)";
                ofs << m_sep;
                ofs << "[14]pzmax(kg*m/s)";
                ofs << m_sep;
                ofs << "[15]gmin()";
                ofs << m_sep;
                ofs << "[16]gmax()";
                ofs << m_sep;
#if (defined WARPX_DIM_3D)
                ofs << "[17]wmin()";
                ofs << m_sep;
                ofs << "[18]wmax()";
#else
                ofs << "[17]wmin(1/m)";
                ofs << m_sep;
                ofs << "[18]wmax(1/m)";
#endif
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
    if (!m_intervals.contains(step+1)) { return; }

    // get MultiParticleContainer class object
    auto & mypc = WarpX::GetInstance().GetPartContainer();

    // get number of species (int)
    const auto nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    const auto species_names = mypc.GetSpeciesNames();

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

        // only chosen species does
        if (species_names[i_s] != m_species_name) { continue; }

        // get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        // get mass (Real)
        auto m = myspc.getMass();
        auto is_photon = myspc.AmIA<PhysicalSpecies::photon>();
        if ( is_photon ) {
            m = PhysConst::m_e;
        }

        using PType = typename WarpXParticleContainer::SuperParticleType;

        // xmin
#if (defined WARPX_DIM_RZ)
        Real xmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0)*std::cos(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMin(xmin);
#else
        Real xmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0); });
        ParallelDescriptor::ReduceRealMin(xmin);
#endif

        // xmax
#if (defined WARPX_DIM_RZ)
        Real xmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0)*std::cos(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMax(xmax);
#else
        Real xmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0); });
        ParallelDescriptor::ReduceRealMax(xmax);
#endif

        // ymin
#if (defined WARPX_DIM_RZ)
        Real ymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0)*std::sin(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMin(ymin);
#elif (defined WARPX_DIM_XZ)
        Real ymin = 0.0_rt;
#else
        Real ymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(1); });
        ParallelDescriptor::ReduceRealMin(ymin);
#endif

        // ymax
#if (defined WARPX_DIM_RZ)
        Real ymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0)*std::sin(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMax(ymax);
#elif (defined WARPX_DIM_XZ)
        Real ymax = 0.0_rt;
#else
        Real ymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(1); });
        ParallelDescriptor::ReduceRealMax(ymax);
#endif

        // zmin
        Real zmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(index_z); });
        ParallelDescriptor::ReduceRealMin(zmin);

        // zmax
        Real zmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(index_z); });
        ParallelDescriptor::ReduceRealMax(zmax);

        // uxmin
        Real uxmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::ux); });
        ParallelDescriptor::ReduceRealMin(uxmin);

        // uxmax
        Real uxmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::ux); });
        ParallelDescriptor::ReduceRealMax(uxmax);

        // uymin
        Real uymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::uy); });
        ParallelDescriptor::ReduceRealMin(uymin);

        // uymax
        Real uymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::uy); });
        ParallelDescriptor::ReduceRealMax(uymax);

        // uzmin
        Real uzmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::uz); });
        ParallelDescriptor::ReduceRealMin(uzmin);

        // uzmax
        Real uzmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::uz); });
        ParallelDescriptor::ReduceRealMax(uzmax);

        // gmin
        Real gmin = 0.0_rt;
        if ( is_photon ) {
            gmin = ReduceMin( myspc,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p)
            {
                Real ux = p.rdata(PIdx::ux);
                Real uy = p.rdata(PIdx::uy);
                Real uz = p.rdata(PIdx::uz);
                Real us = ux*ux + uy*uy + uz*uz;
                return std::sqrt(us*inv_c2);
            });
        } else {
            gmin = ReduceMin( myspc,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p)
            {
                Real ux = p.rdata(PIdx::ux);
                Real uy = p.rdata(PIdx::uy);
                Real uz = p.rdata(PIdx::uz);
                Real us = ux*ux + uy*uy + uz*uz;
                return std::sqrt(1.0_rt + us*inv_c2);
            });
        }
        ParallelDescriptor::ReduceRealMin(gmin);

        // gmax
        Real gmax = 0.0_rt;
        if ( is_photon ) {
            gmax = ReduceMax( myspc,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p)
            {
                Real ux = p.rdata(PIdx::ux);
                Real uy = p.rdata(PIdx::uy);
                Real uz = p.rdata(PIdx::uz);
                Real us = ux*ux + uy*uy + uz*uz;
                return std::sqrt(us*inv_c2);
            });
        } else {
            gmax = ReduceMax( myspc,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p)
            {
                Real ux = p.rdata(PIdx::ux);
                Real uy = p.rdata(PIdx::uy);
                Real uz = p.rdata(PIdx::uz);
                Real us = ux*ux + uy*uy + uz*uz;
                return std::sqrt(1.0_rt + us*inv_c2);
            });
        }
        ParallelDescriptor::ReduceRealMax(gmax);

        // wmin
        Real wmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::w); });
        ParallelDescriptor::ReduceRealMin(wmin);

        // wmax
        Real wmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.rdata(PIdx::w); });
        ParallelDescriptor::ReduceRealMax(wmax);

#if (defined WARPX_QED)
        // get number of level (int)
        const auto level_number = WarpX::GetInstance().finestLevel();

        // compute chimin and chimax
        Real chimin_f = 0.0_rt;
        Real chimax_f = 0.0_rt;
        GetExternalEField get_externalE;
        GetExternalBField get_externalB;

        if (myspc.has_breit_wheeler() || myspc.has_quantum_sync())
        {

            // declare chi arrays
            std::vector<Real> chimin, chimax;
            chimin.resize(level_number+1,0.0_rt);
            chimax.resize(level_number+1,0.0_rt);

            // define variables in preparation for field gathering
            auto & warpx = WarpX::GetInstance();
            const int n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;
            const int nox = WarpX::nox;
            const bool galerkin_interpolation = WarpX::galerkin_interpolation;
            const amrex::IntVect ngE = warpx.getngE();
            const amrex::Array<amrex::Real,3> v_galilean = myspc.get_v_galilean();
            const auto& time_of_last_gal_shift = warpx.time_of_last_gal_shift;

            // loop over refinement levels
            for (int lev = 0; lev <= level_number; ++lev)
            {
                // define variables in preparation for field gathering
                const amrex::Real cur_time = WarpX::GetInstance().gett_new(lev);
                const amrex::Real time_shift = (cur_time - time_of_last_gal_shift);
                const amrex::Array<amrex::Real,3> galilean_shift = { v_galilean[0]*time_shift, v_galilean[1]*time_shift, v_galilean[2]*time_shift };
                const std::array<amrex::Real,3>& dx = WarpX::CellSize(std::max(lev, 0));
                const GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
                const MultiFab & Ex = warpx.getEfield(lev,0);
                const MultiFab & Ey = warpx.getEfield(lev,1);
                const MultiFab & Ez = warpx.getEfield(lev,2);
                const MultiFab & Bx = warpx.getBfield(lev,0);
                const MultiFab & By = warpx.getBfield(lev,1);
                const MultiFab & Bz = warpx.getBfield(lev,2);

                // Loop over boxes
                for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
                {
                    const auto GetPosition = GetParticlePosition(pti);
                    // get particle arrays
                    amrex::ParticleReal* const AMREX_RESTRICT ux = pti.GetAttribs()[PIdx::ux].dataPtr();
                    amrex::ParticleReal* const AMREX_RESTRICT uy = pti.GetAttribs()[PIdx::uy].dataPtr();
                    amrex::ParticleReal* const AMREX_RESTRICT uz = pti.GetAttribs()[PIdx::uz].dataPtr();
                    // declare external fields
                    const int offset = 0;
                    const auto getExternalE = GetExternalEField(pti, offset);
                    const auto getExternalB = GetExternalBField(pti, offset);

                    // define variables in preparation for field gathering
                    amrex::Box box = pti.tilebox();
                    box.grow(ngE);
                    const Dim3 lo = amrex::lbound(box);
                    const std::array<amrex::Real, 3>& xyzmin = WarpX::LowerCorner(box, galilean_shift, lev);
                    const GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};
                    const auto& ex_arr = Ex[pti].array();
                    const auto& ey_arr = Ey[pti].array();
                    const auto& ez_arr = Ez[pti].array();
                    const auto& bx_arr = Bx[pti].array();
                    const auto& by_arr = By[pti].array();
                    const auto& bz_arr = Bz[pti].array();
                    const IndexType ex_type = Ex[pti].box().ixType();
                    const IndexType ey_type = Ey[pti].box().ixType();
                    const IndexType ez_type = Ez[pti].box().ixType();
                    const IndexType bx_type = Bx[pti].box().ixType();
                    const IndexType by_type = By[pti].box().ixType();
                    const IndexType bz_type = Bz[pti].box().ixType();

                    // declare reduce_op
                    ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op;
                    ReduceData<Real, Real> reduce_data(reduce_op);
                    using ReduceTuple = typename decltype(reduce_data)::Type;
                    reduce_op.eval(pti.numParticles(), reduce_data,
                    [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                    {
                        // get external fields
                        ParticleReal xp, yp, zp;
                        GetPosition(i, xp, yp, zp);
                        ParticleReal ex = 0._rt, ey = 0._rt, ez = 0._rt;
                        ParticleReal bx = 0._rt, by = 0._rt, bz = 0._rt;
                        getExternalE(i, ex, ey, ez);
                        getExternalB(i, bx, by, bz);

                        // gather E and B
                        doGatherShapeN(xp, yp, zp,
                            ex, ey, ez, bx, by, bz,
                            ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                            ex_type, ey_type, ez_type,
                            bx_type, by_type, bz_type,
                            dx_arr, xyzmin_arr, lo,
                            n_rz_azimuthal_modes, nox, galerkin_interpolation);
                        // compute chi
                        Real chi = 0.0_rt;
                        if ( is_photon ) {
                            chi = QedUtils::chi_photon(ux[i]*m, uy[i]*m, uz[i]*m,
                                             ex, ey, ez, bx, by, bz);
                        } else {
                            chi = QedUtils::chi_ele_pos(ux[i]*m, uy[i]*m, uz[i]*m,
                                             ex, ey, ez, bx, by, bz);
                        }
                        return {chi,chi};
                    });
                    chimin[lev] = get<0>(reduce_data.value());
                    chimax[lev] = get<1>(reduce_data.value());
                }
                chimin_f = *std::min_element(chimin.begin(), chimin.end());
                chimax_f = *std::max_element(chimax.begin(), chimax.end());
            }
            ParallelDescriptor::ReduceRealMin(chimin_f);
            ParallelDescriptor::ReduceRealMax(chimax_f);
        }
#endif

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
#if (defined WARPX_QED)
        if (myspc.has_breit_wheeler() || myspc.has_quantum_sync())
        {
            m_data[16] = chimin_f;
            m_data[17] = chimax_f;
        }
#endif

    }
    // end loop over species

}
// end void ParticleEnergy::ComputeDiags
