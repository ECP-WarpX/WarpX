/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ColliderRelevant.H"

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
#include "Utils/TextMsg.H"
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
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>
#include <AMReX_Vector.H>

#include <ablastr/coarsen/sample.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <vector>

using namespace amrex;

// constructor
ColliderRelevant::ColliderRelevant (std::string rd_name)
: ReducedDiags{rd_name}
{
    // read colliding species names - must be 2
    ParmParse pp_rd_name(rd_name);
    pp_rd_name.getarr("species", m_beam_name);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_beam_name.size() == 2u,
        "Collider-relevant diagnostic must involve exactly two species"
    );

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    const auto & mypc =  warpx.GetPartContainer();

    // get species names (std::vector<std::string>)
    const auto species_names = mypc.GetSpeciesNames();

    // loop over species
    for (int i_s = 0; i_s < 2; ++i_s)
    {
        // get WarpXParticleContainer class object
        auto const &myspc = mypc.GetParticleContainerFromName(m_beam_name[i_s]);

        auto is_photon = myspc.AmIA<PhysicalSpecies::photon>();

        // photon number density is not available yet
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !is_photon,
            "Collider-relevant diagnostic does not work for colliding photons yet"
        );
    }

    auto all_diag_names = std::vector<std::string>{};
    auto add_diag = [&,c=0](
        const std::string& name, const std::string& header) mutable {
        m_headers_indices[name] = aux_header_index{header, c++};
        all_diag_names.push_back(name);
    };

#if (defined WARPX_DIM_3D)
    add_diag("lumi", "lumi(m^-2*s^-1)");
#elif (defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ))
    add_diag("lumi", "lumi(m^-1*s^-1)");
#else
    add_diag("lumi", "lumi(s^-1)");
#endif

    // loop over species
    for (int i_s = 0; i_s < 2; ++i_s)
    {
        // get WarpXParticleContainer class object
        auto const &myspc = mypc.GetParticleContainerFromName(m_beam_name[i_s]);

        if (myspc.DoQED()){
            add_diag("chimin_"+species_names[i_s], "chimin_"+species_names[i_s]+"()");
            add_diag("chiave_"+species_names[i_s], "chiave_"+species_names[i_s]+"()");
            add_diag("chimax_"+species_names[i_s], "chimax_"+species_names[i_s]+"()");
        }
#if (defined WARPX_DIM_3D)
        add_diag("x_std_"+species_names[i_s], "x_std_"+species_names[i_s]+"(m)");
        add_diag("y_std_"+species_names[i_s], "y_std_"+species_names[i_s]+"(m)");
#elif (defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ))
        add_diag("x_std_"+species_names[i_s], "x_std_"+species_names[i_s]+"(m)");
#endif

        m_data.resize(all_diag_names.size());
    }
    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_write_header )
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
            ofs << m_sep;
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
// end constructor

// function that compute beam relevant quantities
void ColliderRelevant::ComputeDiags (int step)
{
    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get MultiParticleContainer class object
    const auto & mypc = WarpX::GetInstance().GetPartContainer();

    // get species names (std::vector<std::string>)
    auto const species_names = mypc.GetSpeciesNames();

    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // get cell size
    amrex::Geometry const & geom = warpx.Geom(0);
#if defined(WARPX_DIM_1D_Z)
        auto dV = geom.CellSize(0);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined(WARPX_DIM_3D)
        auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

    const auto get_idx = [&](const std::string& name){
        return m_headers_indices.at(name).idx;
    };

    std::array<std::unique_ptr<amrex::MultiFab>, 2> num_dens;

    // loop over species
    for (int i_s = 0; i_s < 2; ++i_s)
    {
        // get WarpXParticleContainer class object
        auto &myspc = mypc.GetParticleContainerFromName(m_beam_name[i_s]);
        // get charge
        ParticleReal const q = myspc.getCharge();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        num_dens[i_s]=myspc.GetChargeDensity(0);
        num_dens[i_s]->mult(1./q);

        // wtot
        Real wtot = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        {
            return p.rdata(PIdx::w); });
        ParallelDescriptor::ReduceRealSum(wtot, ParallelDescriptor::IOProcessorNumber());

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        // x_ave
        Real x_ave = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        {
            const amrex::Real w  = p.rdata(PIdx::w);
            const amrex::Real x = p.pos(0);
            return w*x; });
        ParallelDescriptor::ReduceRealSum(x_ave, ParallelDescriptor::IOProcessorNumber());
        x_ave = x_ave / wtot;

        // x_std
        Real x_std = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        {
            const amrex::Real w  = p.rdata(PIdx::w);
            const amrex::Real x = p.pos(0);
            const amrex::Real tmp = (x - x_ave)*(x - x_ave)*w;
            return tmp; });
        ParallelDescriptor::ReduceRealSum(x_std, ParallelDescriptor::IOProcessorNumber());
        x_std = std::sqrt(x_std / wtot);
#elif defined(WARPX_DIM_3D)
        // x_ave
        Real x_ave = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        {
            const amrex::Real w  = p.rdata(PIdx::w);
            const amrex::Real x = p.pos(0);
            return w*x; });
        ParallelDescriptor::ReduceRealSum(x_ave, ParallelDescriptor::IOProcessorNumber());
        x_ave = x_ave / wtot;

        // x_std
        Real x_std = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        {
            const amrex::Real w  = p.rdata(PIdx::w);
            const amrex::Real x = p.pos(0);
            const amrex::Real tmp = (x - x_ave)*(x - x_ave)*w;
            return tmp; });
        ParallelDescriptor::ReduceRealSum(x_std, ParallelDescriptor::IOProcessorNumber());
        x_std = std::sqrt(x_std / wtot);

        // y_ave
        Real y_ave = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        {
            const amrex::Real w  = p.rdata(PIdx::w);
            const amrex::Real y = p.pos(1);
            return w*y; });
        ParallelDescriptor::ReduceRealSum(y_ave, ParallelDescriptor::IOProcessorNumber());
        y_ave = y_ave / wtot;

        // y_std
        Real y_std = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p)
        {
            const amrex::Real w  = p.rdata(PIdx::w);
            const amrex::Real y = p.pos(1);
            const amrex::Real tmp = (y - y_ave)*(y - y_ave)*w;
            return tmp; });
        ParallelDescriptor::ReduceRealSum(y_std, ParallelDescriptor::IOProcessorNumber());
        y_std = std::sqrt(y_std / wtot);

        //m_data[get_idx("x_ave_"+species_names[i_s])] = x_ave;
        m_data[get_idx("x_std_"+species_names[i_s])] = x_std;
        //m_data[get_idx("y_ave_"+species_names[i_s])] = y_ave;
        m_data[get_idx("y_std_"+species_names[i_s])] = y_std;
#endif

#if (defined WARPX_QED)
        // get number of level (int)
        const auto level_number = WarpX::GetInstance().finestLevel();

        // get mass
        amrex::ParticleReal m = myspc.getMass();
        auto is_photon = myspc.AmIA<PhysicalSpecies::photon>();
        if (is_photon) {
            m = PhysConst::m_e;
        }

        // compute chimin, chiave and chimax
        Real chimin_f = 0.0_rt;
        Real chimax_f = 0.0_rt;
        Real chiave_f = 0.0_rt;

        if (myspc.DoQED())
        {
            // declare chi arrays
            std::vector<Real> chimin, chiave, chimax;
            chimin.resize(level_number+1,0.0_rt);
            chimax.resize(level_number+1,0.0_rt);
            chiave.resize(level_number+1,0.0_rt);

            // define variables in preparation for field gatheeduce_data.value()ring
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
                    amrex::ParticleReal* const AMREX_RESTRICT w = pti.GetAttribs()[PIdx::w].dataPtr();
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
                    ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpSum> reduce_op;
                    ReduceData<Real, Real, Real> reduce_data(reduce_op);
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
                        //amrex::AllPrint() << "CHI DOT W " << chi << " " << w[i] << " " << chi*w[i] <<  "   \n";
                        return {chi, chi, chi*w[i]};
                    });

                    auto hv = reduce_data.value();

                    chimin[lev] = get<0>(hv);
                    chimax[lev] = get<1>(hv);
                    chiave[lev] = get<2>(hv);
                }
            }
            chimin_f = *std::min_element(chimin.begin(), chimin.end());
            chimax_f = *std::max_element(chimax.begin(), chimax.end());
            chiave_f = chiave[0]; // FIXME mesh refinement

            ParallelDescriptor::ReduceRealMin(chimin_f, ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::ReduceRealMax(chimax_f, ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::ReduceRealSum(chiave_f, ParallelDescriptor::IOProcessorNumber());

            m_data[get_idx("chimin_"+species_names[i_s])] = chimin_f;
            m_data[get_idx("chiave_"+species_names[i_s])] = chiave_f/wtot;
            m_data[get_idx("chimax_"+species_names[i_s])] = chimax_f;
        }
#endif
    } // end loop over species

    // make density MultiFabs from nodal to cell centered
    amrex::BoxArray ba = warpx.boxArray(0);
    amrex::DistributionMapping dmap = warpx.DistributionMap(0);
    constexpr int ncomp = 1;
    constexpr int ngrow = 0;
    amrex::MultiFab mf_dst1(ba.convert(amrex::IntVect::TheCellVector()), dmap, ncomp, ngrow);
    amrex::MultiFab mf_dst2(ba.convert(amrex::IntVect::TheCellVector()), dmap, ncomp, ngrow);
    ablastr::coarsen::sample::Coarsen(mf_dst1, *num_dens[0], 0, 0, ncomp, ngrow);
    ablastr::coarsen::sample::Coarsen(mf_dst2, *num_dens[1], 0, 0, ncomp, ngrow);

    // compute luminosity
    auto const n1_dot_n2 = amrex::MultiFab::Dot(mf_dst1, 0, mf_dst2, 0, 1, 0);
    auto const lumi = 2. * PhysConst::c * n1_dot_n2 * dV;

    m_data[get_idx("lumi")] = lumi;
}
// end void ColliderRelevant::ComputeDiags
