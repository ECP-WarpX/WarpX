/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Yinjian Zhao, Arianna Formenti
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
#include <ablastr/warn_manager/WarnManager.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <vector>

using namespace amrex;

ColliderRelevant::ColliderRelevant (std::string rd_name)
: ReducedDiags{rd_name}
{
    // read colliding species names - must be 2
    amrex::ParmParse pp_rd_name(rd_name);
    pp_rd_name.getarr("species", m_beam_name);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_beam_name.size() == 2u,
        "Collider-relevant diagnostics must involve exactly two species");

    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "Collider-relevant diagnostics do not work in RZ geometry.");
#endif

    ablastr::warn_manager::WMRecordWarning(
        "DIAGNOSTICS",
        "The collider-relevant reduced diagnostic is meant for \
        colliding species propagating along the z direction.",
        ablastr::warn_manager::WarnPriority::low);

    ablastr::warn_manager::WMRecordWarning(
        "DIAGNOSTICS",
        "The collider-relevant reduced diagnostic only considers the \
        coarsest level of refinement for the calculations involving chi.",
        ablastr::warn_manager::WarnPriority::low);

    // get WarpX class object
    auto& warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    const MultiParticleContainer& mypc =  warpx.GetPartContainer();

    // get species names (std::vector<std::string>)
    const std::vector<std::string> species_names = mypc.GetSpeciesNames();

    // loop over species
    for (int i_s = 0; i_s < 2; ++i_s)
    {
        // get WarpXParticleContainer class object
        const WarpXParticleContainer& myspc = mypc.GetParticleContainerFromName(m_beam_name[i_s]);

        const bool is_photon = myspc.AmIA<PhysicalSpecies::photon>();

        // photon number density is not available yet
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !is_photon,
            "Collider-relevant diagnostic does not work for colliding photons yet");
    }

    auto all_diag_names = std::vector<std::string>{};
    auto add_diag = [&,c=0](
        const std::string& name, const std::string& header) mutable {
        m_headers_indices[name] = aux_header_index{header, c++};
        all_diag_names.push_back(name);
    };

#if (defined WARPX_DIM_3D)
    add_diag("dL_dt", "dL_dt(m^-2*s^-1)");
#elif (defined WARPX_DIM_XZ)
    add_diag("dL_dt", "dL_dt(m^-1*s^-1)");
#else
    add_diag("dL_dt", "dL_dt(s^-1)");
#endif

    // loop over species
    for (int i_s = 0; i_s < 2; ++i_s)
    {
        // get WarpXParticleContainer class object
        const WarpXParticleContainer& myspc = mypc.GetParticleContainerFromName(m_beam_name[i_s]);

        if (myspc.DoQED()){
            add_diag("chimin_"+species_names[i_s], "chi_min_"+species_names[i_s]+"()");
            add_diag("chiave_"+species_names[i_s], "chi_ave_"+species_names[i_s]+"()");
            add_diag("chimax_"+species_names[i_s], "chi_max_"+species_names[i_s]+"()");
        }
#if (defined WARPX_DIM_3D)
        add_diag("x_ave_"+species_names[i_s], "x_ave_"+species_names[i_s]+"(m)");
        add_diag("x_std_"+species_names[i_s], "x_std_"+species_names[i_s]+"(m)");
        add_diag("y_ave_"+species_names[i_s], "y_ave_"+species_names[i_s]+"(m)");
        add_diag("y_std_"+species_names[i_s], "y_std_"+species_names[i_s]+"(m)");
        add_diag("thetax_min_"+species_names[i_s], "theta_x_min_"+species_names[i_s]+"(rad)");
        add_diag("thetax_ave_"+species_names[i_s], "theta_x_ave_"+species_names[i_s]+"(rad)");
        add_diag("thetax_max_"+species_names[i_s], "theta_x_max_"+species_names[i_s]+"(rad)");
        add_diag("thetax_std_"+species_names[i_s], "theta_x_std_"+species_names[i_s]+"(rad)");
        add_diag("thetay_min_"+species_names[i_s], "theta_y_min_"+species_names[i_s]+"(rad)");
        add_diag("thetay_ave_"+species_names[i_s], "theta_y_ave_"+species_names[i_s]+"(rad)");
        add_diag("thetay_max_"+species_names[i_s], "theta_y_max_"+species_names[i_s]+"(rad)");
        add_diag("thetay_std_"+species_names[i_s], "theta_y_std_"+species_names[i_s]+"(rad)");

#elif (defined WARPX_DIM_XZ)
        add_diag("x_ave_"+species_names[i_s], "x_ave_"+species_names[i_s]+"(m)");
        add_diag("x_std_"+species_names[i_s], "x_std_"+species_names[i_s]+"(m)");
        add_diag("thetax_min_"+species_names[i_s], "theta_x_min_"+species_names[i_s]+"(rad)");
        add_diag("thetax_ave_"+species_names[i_s], "theta_x_ave_"+species_names[i_s]+"(rad)");
        add_diag("thetax_max_"+species_names[i_s], "theta_x_max_"+species_names[i_s]+"(rad)");
        add_diag("thetax_std_"+species_names[i_s], "theta_x_std_"+species_names[i_s]+"(rad)");
#endif
        m_data.resize(all_diag_names.size());
    }

    if (amrex::ParallelDescriptor::IOProcessor())
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

void ColliderRelevant::ComputeDiags (int step)
{
#if defined(WARPX_DIM_RZ)
    amrex::ignore_unused(step);
#else
    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get MultiParticleContainer class object
    const MultiParticleContainer& mypc = WarpX::GetInstance().GetPartContainer();

    // get species names (std::vector<std::string>)
    const std::vector<std::string> species_names = mypc.GetSpeciesNames();

    // get a reference to WarpX instance
    auto& warpx = WarpX::GetInstance();

    // get cell size
    amrex::Geometry const & geom = warpx.Geom(0);
#if defined(WARPX_DIM_1D_Z)
        amrex::Real dV = geom.CellSize(0);
#elif defined(WARPX_DIM_XZ)
        amrex::Real dV = geom.CellSize(0) * geom.CellSize(1);
#elif defined(WARPX_DIM_3D)
        amrex::Real dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

    const auto get_idx = [&](const std::string& name){
        return m_headers_indices.at(name).idx;
    };

    std::array<std::unique_ptr<amrex::MultiFab>, 2> num_dens;

    // loop over species
    for (int i_s = 0; i_s < 2; ++i_s)
    {
        // get WarpXParticleContainer class object
        WarpXParticleContainer& myspc = mypc.GetParticleContainerFromName(m_beam_name[i_s]);
        // get charge
        amrex::ParticleReal const q = myspc.getCharge();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        num_dens[i_s] = myspc.GetChargeDensity(0);
        num_dens[i_s]->mult(1./q);

        // wtot
        amrex::Real wtot = ReduceSum( myspc,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p)
            {
                return p.rdata(PIdx::w);
            });
        amrex::ParallelDescriptor::ReduceRealSum(wtot, amrex::ParallelDescriptor::IOProcessorNumber());

#if defined(WARPX_DIM_XZ)
        // thetax_min, thetax_max
        amrex::Real thetax_min = 0.0_rt;
        amrex::Real thetax_max = 0.0_rt;

        amrex::ReduceOps<ReduceOpMin, ReduceOpMax> reduce_ops_minmax;
        auto r_minmax = amrex::ParticleReduce<amrex::ReduceData<Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real>
            {
                const amrex::Real ux = p.rdata(PIdx::ux);
                const amrex::Real uz = p.rdata(PIdx::uz);
                const amrex::Real thetax = std::atan2(ux, uz);
                return {thetax, thetax};
            },
            reduce_ops_minmax);
        thetax_min = amrex::get<0>(r_minmax);
        thetax_max = amrex::get<1>(r_minmax);
        amrex::ParallelDescriptor::ReduceRealMin(thetax_min, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::ReduceRealMax(thetax_max, amrex::ParallelDescriptor::IOProcessorNumber());

        // x_ave, thetax_ave
        amrex::Real x_ave = 0.0_rt;
        amrex::Real thetax_ave = 0.0_rt;

        amrex::ReduceOps<ReduceOpSum, ReduceOpSum> reduce_ops_ave;
        auto r_ave = amrex::ParticleReduce<amrex::ReduceData<Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real>
            {
                const amrex::Real w  = p.rdata(PIdx::w);
                const amrex::Real x = p.pos(0);
                const amrex::Real ux = p.rdata(PIdx::ux);
                const amrex::Real uz = p.rdata(PIdx::uz);
                const amrex::Real thetax = std::atan2(ux, uz);
                return {w*x, w*thetax};
            },
            reduce_ops_ave);
        x_ave = amrex::get<0>(r_ave);
        thetax_ave = amrex::get<1>(r_ave);
        amrex::ParallelDescriptor::ReduceRealSum({x_ave, thetax_ave}, amrex::ParallelDescriptor::IOProcessorNumber());
        x_ave = x_ave / wtot;
        thetax_ave = thetax_ave / wtot;

        // x_std, y_std, thetax_std, thetay_std
        amrex::Real x_std = 0.0_rt;
        amrex::Real thetax_std = 0.0_rt;

        amrex::ReduceOps<ReduceOpSum, ReduceOpSum> reduce_ops_std;
        auto r_std = amrex::ParticleReduce<amrex::ReduceData<Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real>
            {
                const amrex::Real w  = p.rdata(PIdx::w);
                const amrex::Real x = p.pos(0);
                const amrex::Real ux = p.rdata(PIdx::ux);
                const amrex::Real uz = p.rdata(PIdx::uz);
                const amrex::Real thetax = std::atan2(ux, uz);
                const amrex::Real tmp1 = (x - x_ave)*(x - x_ave)*w;
                const amrex::Real tmp2 = (thetax - thetax_ave)*(thetax - thetax_ave)*w;
                return {tmp1, tmp2};
            },
            reduce_ops_std);
        x_std = amrex::get<0>(r_std);
        thetax_std = amrex::get<1>(r_std);
        amrex::ParallelDescriptor::ReduceRealSum({x_std, thetax_std}, amrex::ParallelDescriptor::IOProcessorNumber());
        x_std = std::sqrt(x_std / wtot);
        thetax_std = std::sqrt(thetax_std / wtot);

        m_data[get_idx("x_ave_"+species_names[i_s])] = x_ave;
        m_data[get_idx("x_std_"+species_names[i_s])] = x_std;
        m_data[get_idx("thetax_min_"+species_names[i_s])] = thetax_min;
        m_data[get_idx("thetax_ave_"+species_names[i_s])] = thetax_ave;
        m_data[get_idx("thetax_max_"+species_names[i_s])] = thetax_max;
        m_data[get_idx("thetax_std_"+species_names[i_s])] = thetax_std;

#elif defined(WARPX_DIM_3D)
        // thetax_min, thetax_max, thetay_min, thetay_max
        amrex::Real thetax_min = 0.0_rt;
        amrex::Real thetax_max = 0.0_rt;
        amrex::Real thetay_min = 0.0_rt;
        amrex::Real thetay_max = 0.0_rt;

        amrex::ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpMin, ReduceOpMax> reduce_ops_minmax;
        auto r_minmax = amrex::ParticleReduce<amrex::ReduceData<Real, Real, Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real, Real, Real>
            {
                const amrex::Real ux = p.rdata(PIdx::ux);
                const amrex::Real uy = p.rdata(PIdx::uy);
                const amrex::Real uz = p.rdata(PIdx::uz);
                const amrex::Real thetax = std::atan2(ux, uz);
                const amrex::Real thetay = std::atan2(uy, uz);
                return {thetax, thetax, thetay, thetay};
            },
            reduce_ops_minmax);
        thetax_min = amrex::get<0>(r_minmax);
        thetax_max = amrex::get<1>(r_minmax);
        thetay_min = amrex::get<2>(r_minmax);
        thetay_max = amrex::get<3>(r_minmax);
        amrex::ParallelDescriptor::ReduceRealMin({thetax_min, thetay_min}, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::ReduceRealMax({thetax_max, thetay_max}, amrex::ParallelDescriptor::IOProcessorNumber());

        // x_ave, y_ave, thetax_ave, thetay_ave
        amrex::Real x_ave = 0.0_rt;
        amrex::Real y_ave = 0.0_rt;
        amrex::Real thetax_ave = 0.0_rt;
        amrex::Real thetay_ave = 0.0_rt;

        amrex::ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_ops_ave;
        auto r_ave = amrex::ParticleReduce<amrex::ReduceData<Real, Real, Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real, Real, Real>
            {
                const amrex::Real w  = p.rdata(PIdx::w);
                const amrex::Real x = p.pos(0);
                const amrex::Real ux = p.rdata(PIdx::ux);
                const amrex::Real y = p.pos(1);
                const amrex::Real uy = p.rdata(PIdx::uy);
                const amrex::Real uz = p.rdata(PIdx::uz);
                const amrex::Real thetax = std::atan2(ux, uz);
                const amrex::Real thetay = std::atan2(uy, uz);
                return {w*x, w*y, w*thetax, w*thetay};
            },
            reduce_ops_ave);
        x_ave = amrex::get<0>(r_ave);
        y_ave = amrex::get<1>(r_ave);
        thetax_ave = amrex::get<2>(r_ave);
        thetay_ave = amrex::get<3>(r_ave);
        amrex::ParallelDescriptor::ReduceRealSum({x_ave, y_ave, thetax_ave, thetay_ave}, amrex::ParallelDescriptor::IOProcessorNumber());
        x_ave = x_ave / wtot;
        y_ave = y_ave / wtot;
        thetax_ave = thetax_ave / wtot;
        thetay_ave = thetay_ave / wtot;

        // x_std, y_std, thetax_std, thetay_std
        amrex::Real x_std = 0.0_rt;
        amrex::Real y_std = 0.0_rt;
        amrex::Real thetax_std = 0.0_rt;
        amrex::Real thetay_std = 0.0_rt;

        amrex::ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_ops_std;
        auto r_std = amrex::ParticleReduce<amrex::ReduceData<Real, Real, Real, Real>>(
            myspc,
            [=] AMREX_GPU_DEVICE(const PType& p) noexcept -> amrex::GpuTuple<Real, Real, Real, Real>
            {
                const amrex::Real w  = p.rdata(PIdx::w);
                const amrex::Real x = p.pos(0);
                const amrex::Real ux = p.rdata(PIdx::ux);
                const amrex::Real y = p.pos(1);
                const amrex::Real uy = p.rdata(PIdx::uy);
                const amrex::Real uz = p.rdata(PIdx::uz);
                const amrex::Real thetax = std::atan2(ux, uz);
                const amrex::Real thetay = std::atan2(uy, uz);
                const amrex::Real tmp1 = (x - x_ave)*(x - x_ave)*w;
                const amrex::Real tmp2 = (y - y_ave)*(y - y_ave)*w;
                const amrex::Real tmp3 = (thetax - thetax_ave)*(thetax - thetax_ave)*w;
                const amrex::Real tmp4 = (thetay - thetay_ave)*(thetay - thetay_ave)*w;
                return {tmp1, tmp2, tmp3, tmp4};
            },
            reduce_ops_std);
        x_std = amrex::get<0>(r_std);
        y_std = amrex::get<1>(r_std);
        thetax_std = amrex::get<2>(r_std);
        thetay_std = amrex::get<3>(r_std);
        amrex::ParallelDescriptor::ReduceRealSum({x_std, y_std, thetax_std, thetay_std}, amrex::ParallelDescriptor::IOProcessorNumber());
        x_std = std::sqrt(x_std / wtot);
        y_std = std::sqrt(y_std / wtot);
        thetax_std = std::sqrt(thetax_std / wtot);
        thetay_std = std::sqrt(thetay_std / wtot);

        m_data[get_idx("x_ave_"+species_names[i_s])] = x_ave;
        m_data[get_idx("x_std_"+species_names[i_s])] = x_std;
        m_data[get_idx("y_ave_"+species_names[i_s])] = y_ave;
        m_data[get_idx("y_std_"+species_names[i_s])] = y_std;
        m_data[get_idx("thetax_min_"+species_names[i_s])] = thetax_min;
        m_data[get_idx("thetax_ave_"+species_names[i_s])] = thetax_ave;
        m_data[get_idx("thetax_max_"+species_names[i_s])] = thetax_max;
        m_data[get_idx("thetax_std_"+species_names[i_s])] = thetax_std;
        m_data[get_idx("thetay_min_"+species_names[i_s])] = thetay_min;
        m_data[get_idx("thetay_ave_"+species_names[i_s])] = thetay_ave;
        m_data[get_idx("thetay_max_"+species_names[i_s])] = thetay_max;
        m_data[get_idx("thetay_std_"+species_names[i_s])] = thetay_std;
#endif

#if (defined WARPX_QED)
        // get mass
        amrex::ParticleReal m = myspc.getMass();
        const bool is_photon = myspc.AmIA<PhysicalSpecies::photon>();
        if (is_photon) {
            m = PhysConst::m_e;
        }

        // compute chimin, chiave and chimax
        amrex::Real chimin_f = 0.0_rt;
        amrex::Real chimax_f = 0.0_rt;
        amrex::Real chiave_f = 0.0_rt;

        if (myspc.DoQED())
        {
            // define variables in preparation for field gatheeduce_data.value()ring
            const int n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;
            const int nox = WarpX::nox;
            const bool galerkin_interpolation = WarpX::galerkin_interpolation;
            const amrex::IntVect ngEB = warpx.getngEB();

            // TO DO: loop over refinement levels: for (int lev = 0; lev <= level_number; ++lev)
            int lev = 0;

            // define variables in preparation for field gathering
            const std::array<amrex::Real,3>& dx = WarpX::CellSize(std::max(lev, 0));
            const amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
            const amrex::MultiFab & Ex = warpx.getEfield(lev,0);
            const amrex::MultiFab & Ey = warpx.getEfield(lev,1);
            const amrex::MultiFab & Ez = warpx.getEfield(lev,2);
            const amrex::MultiFab & Bx = warpx.getBfield(lev,0);
            const amrex::MultiFab & By = warpx.getBfield(lev,1);
            const amrex::MultiFab & Bz = warpx.getBfield(lev,2);

            // declare reduce_op
            ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpSum> reduce_op;
            ReduceData<Real, Real, Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

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
                const amrex::Dim3 lo = amrex::lbound(box);
                const std::array<amrex::Real, 3>& xyzmin = WarpX::LowerCorner(box, lev, 0._rt);
                const amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};
                const amrex::Array4<const amrex::Real> & ex_arr = Ex[pti].array();
                const amrex::Array4<const amrex::Real> & ey_arr = Ey[pti].array();
                const amrex::Array4<const amrex::Real> & ez_arr = Ez[pti].array();
                const amrex::Array4<const amrex::Real> & bx_arr = Bx[pti].array();
                const amrex::Array4<const amrex::Real> & by_arr = By[pti].array();
                const amrex::Array4<const amrex::Real> & bz_arr = Bz[pti].array();
                const amrex::IndexType ex_type = Ex[pti].box().ixType();
                const amrex::IndexType ey_type = Ey[pti].box().ixType();
                const amrex::IndexType ez_type = Ez[pti].box().ixType();
                const amrex::IndexType bx_type = Bx[pti].box().ixType();
                const amrex::IndexType by_type = By[pti].box().ixType();
                const amrex::IndexType bz_type = Bz[pti].box().ixType();

                // evaluate reduce_op
                reduce_op.eval(pti.numParticles(), reduce_data,
                [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                {
                    // get external fields
                    amrex::ParticleReal xp, yp, zp;
                    GetPosition(i, xp, yp, zp);
                    amrex::ParticleReal ex = 0._rt, ey = 0._rt, ez = 0._rt;
                    amrex::ParticleReal bx = 0._rt, by = 0._rt, bz = 0._rt;
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
                    amrex::Real chi = 0.0_rt;
                    if (is_photon) {
                        chi = QedUtils::chi_photon(ux[i]*m, uy[i]*m, uz[i]*m,
                                            ex, ey, ez, bx, by, bz);
                    } else {
                        chi = QedUtils::chi_ele_pos(ux[i]*m, uy[i]*m, uz[i]*m,
                                            ex, ey, ez, bx, by, bz);
                    }
                    return {chi, chi, chi*w[i]};
                });
            }
            auto val = reduce_data.value();
            chimin_f = get<0>(val);
            chimax_f = get<1>(val);
            chiave_f = get<2>(val);
            amrex::ParallelDescriptor::ReduceRealMin(chimin_f, amrex::ParallelDescriptor::IOProcessorNumber());
            amrex::ParallelDescriptor::ReduceRealMax(chimax_f, amrex::ParallelDescriptor::IOProcessorNumber());
            amrex::ParallelDescriptor::ReduceRealSum(chiave_f, amrex::ParallelDescriptor::IOProcessorNumber());

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
    amrex::Real const n1_dot_n2 = amrex::MultiFab::Dot(mf_dst1, 0, mf_dst2, 0, 1, 0);
    amrex::Real const lumi = 2. * PhysConst::c * n1_dot_n2 * dV;
    m_data[get_idx("dL_dt")] = lumi;
#endif // not RZ
}
