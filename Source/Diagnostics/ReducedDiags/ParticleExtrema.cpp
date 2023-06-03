/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleExtrema.H"

#include "Diagnostics/ReducedDiags/ReducedDiags.H"
#if (defined WARPX_QED)
#   include "Particles/ElementaryProcess/QEDInternals/QedChiFunctions.H"
#endif
#include "Particles/Gather/FieldGather.H"
#include "Particles/Gather/GetExternalFields.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/SpeciesPhysicalProperties.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Dim3.H>
#include <AMReX_Extension.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <vector>

using namespace amrex;

// constructor
ParticleExtrema::ParticleExtrema (std::string rd_name)
: ReducedDiags{rd_name}
{
    // read species name
    const ParmParse pp_rd_name(rd_name);
    pp_rd_name.get("species",m_species_name);

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

        auto all_diag_names = std::vector<std::string> {};
        auto add_diag = [&,c=0] (
            const std::string& name, const std::string& header) mutable {
            m_headers_indices[name] = aux_header_index{header, c++};
            all_diag_names.push_back(name);
        };

        add_diag("xmin", "xmin(m)");
        add_diag("xmax", "xmax(m)");
        add_diag("ymin", "ymin(m)");
        add_diag("ymax", "ymax(m)");
        add_diag("zmin", "zmin(m)");
        add_diag("zmax", "zmax(m)");
        add_diag("pxmin", "pxmin(kg*m/s)");
        add_diag("pxmax", "pxmax(kg*m/s)");
        add_diag("pymin", "pymin(kg*m/s)");
        add_diag("pymax", "pymax(kg*m/s)");
        add_diag("pzmin", "pzmin(kg*m/s)");
        add_diag("pzmax", "pzmax(kg*m/s)");
        add_diag("gmin", "gmin()");
        add_diag("gmax", "gmax()");

#if (defined WARPX_DIM_3D)
        add_diag("wmin", "wmin()");
        add_diag("wmax", "wmax()");
#elif (defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ))
        add_diag("wmin", "wmin(1/m)");
        add_diag("wmax", "wmax(1/m)");
#else
        add_diag("wmin", "wmin(1/m^2)");
        add_diag("wmax", "wmax(1/m^2)");
#endif
        if (myspc.DoQED()){
            add_diag("chimin", "chimin()");
            add_diag("chimax", "chimax()");
        }

        m_data.resize(all_diag_names.size());

        if (ParallelDescriptor::IOProcessor())
        {
            if ( m_IsNotRestart )
            {
                // open file
                std::ofstream ofs;
                ofs.open(m_path + m_rd_name + "." + m_extension,
                    std::ofstream::out | std::ofstream::app);
                // write header row
                int off = 0;
                ofs << "#";
                ofs << "[" << off++ << "]step()";
                ofs << m_sep;
                ofs << "[" << off++ << "]time(s)";
                for (const auto& name : all_diag_names){
                    const auto& el = m_headers_indices[name];
                    ofs << m_sep << "[" << el.idx + off << "]" << el.header;
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
#elif (defined WARPX_DIM_1D_Z)
    int const index_z = 0;
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
#elif (defined WARPX_DIM_1D_Z)
        const Real xmin = 0.0_rt;
#else
        Real xmin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0); });
        ParallelDescriptor::ReduceRealMin(xmin);
#endif

        // xmax
#if (defined WARPX_DIM_RZ)
        Real t_xmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0)*std::cos(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMax(t_xmax);
        const Real xmax = t_xmax;
#elif (defined WARPX_DIM_1D_Z)
        const Real xmax = 0.0_rt;
#else
        Real t_xmax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0); });
        ParallelDescriptor::ReduceRealMax(t_xmax);
        const Real xmax = t_xmax;
#endif

        // ymin
#if (defined WARPX_DIM_RZ)
        Real t_ymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0)*std::sin(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMin(t_ymin);
        const Real ymin = t_ymin;
#elif (defined WARPX_DIM_XZ || WARPX_DIM_1D_Z)
        const Real ymin = 0.0_rt;
#else
        Real t_ymin = ReduceMin( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(1); });
        ParallelDescriptor::ReduceRealMin(t_ymin);
        const Real ymin = t_ymin;
#endif

        // ymax
#if (defined WARPX_DIM_RZ)
        Real t_ymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(0)*std::sin(p.rdata(PIdx::theta)); });
        ParallelDescriptor::ReduceRealMax(t_ymax);
        const Real ymax = t_ymax;
#elif (defined WARPX_DIM_XZ || WARPX_DIM_1D_Z)
        const Real ymax = 0.0_rt;
#else
        Real t_ymax = ReduceMax( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        { return p.pos(1); });
        ParallelDescriptor::ReduceRealMax(t_ymax);
        const Real ymax = t_ymax;
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
                const Real ux = p.rdata(PIdx::ux);
                const Real uy = p.rdata(PIdx::uy);
                const Real uz = p.rdata(PIdx::uz);
                const Real us = ux*ux + uy*uy + uz*uz;
                return std::sqrt(us*inv_c2);
            });
        } else {
            gmin = ReduceMin( myspc,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p)
            {
                const Real ux = p.rdata(PIdx::ux);
                const Real uy = p.rdata(PIdx::uy);
                const Real uz = p.rdata(PIdx::uz);
                const Real us = ux*ux + uy*uy + uz*uz;
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
                const Real ux = p.rdata(PIdx::ux);
                const Real uy = p.rdata(PIdx::uy);
                const Real uz = p.rdata(PIdx::uz);
                const Real us = ux*ux + uy*uy + uz*uz;
                return std::sqrt(us*inv_c2);
            });
        } else {
            gmax = ReduceMax( myspc,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p)
            {
                const Real ux = p.rdata(PIdx::ux);
                const Real uy = p.rdata(PIdx::uy);
                const Real uz = p.rdata(PIdx::uz);
                const Real us = ux*ux + uy*uy + uz*uz;
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

        if (myspc.DoQED())
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
            const amrex::IntVect ngEB = warpx.getngEB();

            // loop over refinement levels
            for (int lev = 0; lev <= level_number; ++lev)
            {
                // define variables in preparation for field gathering
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
                    const auto getExternalEB = GetExternalEBField(pti, offset);
                    // define variables in preparation for field gathering
                    amrex::Box box = pti.tilebox();
                    box.grow(ngEB);
                    const Dim3 lo = amrex::lbound(box);
                    const std::array<amrex::Real, 3>& xyzmin = WarpX::LowerCorner(box, lev, 0._rt);
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
                        getExternalEB(i, ex, ey, ez, bx, by, bz);

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

        const auto get_idx = [&](const std::string& name){
            return m_headers_indices.at(name).idx;
        };

        m_data[get_idx("xmin")]  = xmin;
        m_data[get_idx("xmax")]  = xmax;
        m_data[get_idx("ymin")]  = ymin;
        m_data[get_idx("ymax")]  = ymax;
        m_data[get_idx("zmin")]  = zmin;
        m_data[get_idx("zmax")]  = zmax;
        m_data[get_idx("pxmin")]  = uxmin*m;
        m_data[get_idx("pxmax")]  = uxmax*m;
        m_data[get_idx("pymin")]  = uymin*m;
        m_data[get_idx("pymax")]  = uymax*m;
        m_data[get_idx("pzmin")] = uzmin*m;
        m_data[get_idx("pzmax")] = uzmax*m;
        m_data[get_idx("gmin")] = gmin;
        m_data[get_idx("gmax")] = gmax;
        m_data[get_idx("wmin")] = wmin;
        m_data[get_idx("wmax")] = wmax;
#if (defined WARPX_QED)
        if (myspc.DoQED())
        {
            m_data[get_idx("chimin")] = chimin_f;
            m_data[get_idx("chimax")] = chimax_f;
        }
#endif
    }
    // end loop over species
}
// end void ParticleEnergy::ComputeDiags
