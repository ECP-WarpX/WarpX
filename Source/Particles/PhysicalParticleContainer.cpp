/* Copyright 2019-2020 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Glenn Richardson, Jean-Luc Vay
 * Ligia Diana Amorim, Luca Fedeli, Maxence Thevenet
 * Michael Rowan, Remi Lehe, Revathi Jambunathan
 * Weiqun Zhang, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PhysicalParticleContainer.H"

#include "MultiParticleContainer.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "Python/WarpXWrappers.h"
#include "Utils/IonizationEnergiesTable.H"
#include "Particles/Gather/FieldGather.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Pusher/CopyParticleAttribs.H"
#include "Particles/Pusher/PushSelector.H"
#include "Particles/Gather/GetExternalFields.H"
#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_Print.H>
#include <AMReX.H>

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif

#include <cmath>
#include <limits>
#include <sstream>
#include <string>

using namespace amrex;

namespace
{
    // Since the user provides the density distribution
    // at t_lab=0 and in the lab-frame coordinates,
    // we need to find the lab-frame position of this
    // particle at t_lab=0, from its boosted-frame coordinates
    // Assuming ballistic motion, this is given by:
    // z0_lab = gamma*( z_boost*(1-beta*betaz_lab) - ct_boost*(betaz_lab-beta) )
    // where betaz_lab is the speed of the particle in the lab frame
    //
    // In order for this equation to be solvable, betaz_lab
    // is explicitly assumed to have no dependency on z0_lab
    //
    // Note that we use the bulk momentum to perform the ballistic correction
    // Assume no z0_lab dependency
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real applyBallisticCorrection(const XDim3& pos, const InjectorMomentum* inj_mom,
                                  Real gamma_boost, Real beta_boost, Real t) noexcept
    {
        const XDim3 u_bulk = inj_mom->getBulkMomentum(pos.x, pos.y, pos.z);
        const Real gamma_bulk = std::sqrt(1._rt +
                  (u_bulk.x*u_bulk.x+u_bulk.y*u_bulk.y+u_bulk.z*u_bulk.z));
        const Real betaz_bulk = u_bulk.z/gamma_bulk;
        const Real z0 = gamma_boost * ( pos.z*(1.0_rt-beta_boost*betaz_bulk)
                             - PhysConst::c*t*(betaz_bulk-beta_boost) );
        return z0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    XDim3 getCellCoords (const GpuArray<Real, AMREX_SPACEDIM>& lo_corner,
                         const GpuArray<Real, AMREX_SPACEDIM>& dx,
                         const XDim3& r, const IntVect& iv) noexcept
    {
        XDim3 pos;
#if (AMREX_SPACEDIM == 3)
        pos.x = lo_corner[0] + (iv[0]+r.x)*dx[0];
        pos.y = lo_corner[1] + (iv[1]+r.y)*dx[1];
        pos.z = lo_corner[2] + (iv[2]+r.z)*dx[2];
#else
        pos.x = lo_corner[0] + (iv[0]+r.x)*dx[0];
        pos.y = 0.0_rt;
#if   defined WARPX_DIM_XZ
        pos.z = lo_corner[1] + (iv[1]+r.y)*dx[1];
#elif defined WARPX_DIM_RZ
        // Note that for RZ, r.y will be theta
        pos.z = lo_corner[1] + (iv[1]+r.z)*dx[1];
#endif
#endif
        return pos;
    }
}

PhysicalParticleContainer::PhysicalParticleContainer (AmrCore* amr_core, int ispecies,
                                                      const std::string& name)
    : WarpXParticleContainer(amr_core, ispecies),
      species_name(name)
{
    BackwardCompatibility();

    plasma_injector = std::make_unique<PlasmaInjector>(species_id, species_name);
    physical_species = plasma_injector->getPhysicalSpecies();
    charge = plasma_injector->getCharge();
    mass = plasma_injector->getMass();

    ParmParse pp(species_name);

    pp.query("boost_adjust_transverse_positions", boost_adjust_transverse_positions);
    pp.query("do_backward_propagation", do_backward_propagation);

    // Initialize splitting
    pp.query("do_splitting", do_splitting);
    pp.query("split_type", split_type);
    pp.query("do_not_deposit", do_not_deposit);
    pp.query("do_not_gather", do_not_gather);
    pp.query("do_not_push", do_not_push);

    pp.query("do_continuous_injection", do_continuous_injection);
    pp.query("initialize_self_fields", initialize_self_fields);
    pp.query("self_fields_required_precision", self_fields_required_precision);
    pp.query("self_fields_max_iters", self_fields_max_iters);
    // Whether to plot back-transformed (lab-frame) diagnostics
    // for this species.
    pp.query("do_back_transformed_diagnostics", do_back_transformed_diagnostics);

    pp.query("do_field_ionization", do_field_ionization);

    pp.query("do_resampling", do_resampling);
    if (do_resampling) m_resampler = Resampling(species_name);

    //check if Radiation Reaction is enabled and do consistency checks
    pp.query("do_classical_radiation_reaction", do_classical_radiation_reaction);
    //if the species is not a lepton, do_classical_radiation_reaction
    //should be false
    WarpXUtilMsg::AlwaysAssert(
        !(do_classical_radiation_reaction &&
        !(AmIA<PhysicalSpecies::electron>() ||
        AmIA<PhysicalSpecies::positron>() )),
        "ERROR: can't enable classical radiation reaction for non lepton species '"
        + species_name + "'."
    );

    //Only Boris pusher is compatible with radiation reaction
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !(do_classical_radiation_reaction &&
        WarpX::particle_pusher_algo != ParticlePusherAlgo::Boris),
        "Radiation reaction can be enabled only if Boris pusher is used");
    //_____________________________

#ifdef WARPX_QED
    pp.query("do_qed", m_do_qed);
    if(m_do_qed){
        //If do_qed is enabled, find out if Quantum Synchrotron process is enabled
        pp.query("do_qed_quantum_sync", m_do_qed_quantum_sync);
        if (m_do_qed_quantum_sync)
            AddRealComp("optical_depth_QSR");
        pp.query("do_qed_breit_wheeler", m_do_qed_breit_wheeler);
        if (m_do_qed_breit_wheeler)
            AddRealComp("optical_depth_BW");
    }

    if(m_do_qed_quantum_sync){
        pp.get("qed_quantum_sync_phot_product_species",
            m_qed_quantum_sync_phot_product_name);
    }


#endif

    // Parse galilean velocity
    ParmParse ppsatd("psatd");
    ppsatd.query("v_galilean", m_v_galilean);
    // Scale the velocity by the speed of light
    for (int i=0; i<3; i++) m_v_galilean[i] *= PhysConst::c;

    // build filter functors
    m_do_random_filter  = pp.query("random_fraction", m_random_fraction);
    m_do_uniform_filter = pp.query("uniform_stride",  m_uniform_stride);
    std::string buf;
    m_do_parser_filter  = pp.query("plot_filter_function(t,x,y,z,ux,uy,uz)", buf);
    if (m_do_parser_filter) {
        std::string function_string = "";
        Store_parserString(pp,"plot_filter_function(t,x,y,z,ux,uy,uz)",
                           function_string);
        m_particle_filter_parser = std::make_unique<ParserWrapper<7>>(
            makeParser(function_string,{"t","x","y","z","ux","uy","uz"}));
    }

}

PhysicalParticleContainer::PhysicalParticleContainer (AmrCore* amr_core)
    : WarpXParticleContainer(amr_core, 0)
{
    plasma_injector = std::make_unique<PlasmaInjector>();
}

void
PhysicalParticleContainer::BackwardCompatibility ()
{
    ParmParse pps(species_name);
    std::vector<std::string> backward_strings;
    if (pps.queryarr("plot_vars", backward_strings)){
        amrex::Abort("<species>.plot_vars is not supported anymore. "
                     "Please use the new syntax for diagnostics, see documentation.");
    }

    int backward_int;
    if (pps.query("plot_species", backward_int)){
        amrex::Abort("<species>.plot_species is not supported anymore. "
                     "Please use the new syntax for diagnostics, see documentation.");
    }
}

void PhysicalParticleContainer::InitData ()
{
    // Init ionization module here instead of in the PhysicalParticleContainer
    // constructor because dt is required
    if (do_field_ionization) {InitIonizationModule();}
    AddParticles(0); // Note - add on level 0
    Redistribute();  // We then redistribute
}

void PhysicalParticleContainer::MapParticletoBoostedFrame (
    Real& x, Real& y, Real& z, std::array<Real, 3>& u)
{
    // Map the particles from the lab frame to the boosted frame.
    // This boosts the particle to the lab frame and calculates
    // the particle time in the boosted frame. It then maps
    // the position to the time in the boosted frame.

    // For now, start with the assumption that this will only happen
    // at the start of the simulation.
    const Real t_lab = 0.;

    const Real uz_boost = WarpX::gamma_boost*WarpX::beta_boost*PhysConst::c;

    // tpr is the particle's time in the boosted frame
    Real tpr = WarpX::gamma_boost*t_lab - uz_boost*z/(PhysConst::c*PhysConst::c);

    // The particle's transformed location in the boosted frame
    Real xpr = x;
    Real ypr = y;
    Real zpr = WarpX::gamma_boost*z - uz_boost*t_lab;

    // transform u and gamma to the boosted frame
    Real gamma_lab = std::sqrt(1._rt + (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/(PhysConst::c*PhysConst::c));
    // u[0] = u[0];
    // u[1] = u[1];
    u[2] = WarpX::gamma_boost*u[2] - uz_boost*gamma_lab;
    Real gammapr = std::sqrt(1._rt + (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/(PhysConst::c*PhysConst::c));

    Real vxpr = u[0]/gammapr;
    Real vypr = u[1]/gammapr;
    Real vzpr = u[2]/gammapr;

    if (do_backward_propagation){
        u[2] = -u[2];
    }

    // Move the particles to where they will be at t = 0 in the boosted frame
    if (boost_adjust_transverse_positions) {
        x = xpr - tpr*vxpr;
        y = ypr - tpr*vypr;
    }

    z = zpr - tpr*vzpr;

}

void
PhysicalParticleContainer::AddGaussianBeam (
    const Real x_m, const Real y_m, const Real z_m,
    const Real x_rms, const Real y_rms, const Real z_rms,
    const Real x_cut, const Real y_cut, const Real z_cut,
    const Real q_tot, long npart,
    const int do_symmetrize) {

    std::mt19937_64 mt(0451);
    std::normal_distribution<double> distx(x_m, x_rms);
    std::normal_distribution<double> disty(y_m, y_rms);
    std::normal_distribution<double> distz(z_m, z_rms);

    // Declare temporary vectors on the CPU
    Gpu::HostVector<ParticleReal> particle_x;
    Gpu::HostVector<ParticleReal> particle_y;
    Gpu::HostVector<ParticleReal> particle_z;
    Gpu::HostVector<ParticleReal> particle_ux;
    Gpu::HostVector<ParticleReal> particle_uy;
    Gpu::HostVector<ParticleReal> particle_uz;
    Gpu::HostVector<ParticleReal> particle_w;
    int np = 0;

    if (ParallelDescriptor::IOProcessor()) {
        // If do_symmetrize, create 4x fewer particles, and
        // Replicate each particle 4 times (x,y) (-x,y) (x,-y) (-x,-y)
        if (do_symmetrize){
            npart /= 4;
        }
        for (long i = 0; i < npart; ++i) {
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
            const Real weight = q_tot/(npart*charge);
            const Real x = distx(mt);
            const Real y = disty(mt);
            const Real z = distz(mt);
#elif (defined WARPX_DIM_XZ)
            const Real weight = q_tot/(npart*charge*y_rms);
            const Real x = distx(mt);
            constexpr Real y = 0._prt;
            const Real z = distz(mt);
#endif
            if (plasma_injector->insideBounds(x, y, z)  &&
                std::abs( x - x_m ) < x_cut * x_rms     &&
                std::abs( y - y_m ) < y_cut * y_rms     &&
                std::abs( z - z_m ) < z_cut * z_rms   ) {
                XDim3 u = plasma_injector->getMomentum(x, y, z);
                u.x *= PhysConst::c;
                u.y *= PhysConst::c;
                u.z *= PhysConst::c;
                if (do_symmetrize){
                    // Add four particles to the beam:
                    CheckAndAddParticle(x, y, z, { u.x, u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                    CheckAndAddParticle(x, -y, z, { u.x, -u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                    CheckAndAddParticle(-x, y, z, { -u.x, u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                    CheckAndAddParticle(-x, -y, z, { -u.x, -u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                } else {
                    CheckAndAddParticle(x, y, z, { u.x, u.y, u.z}, weight,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                }
            }
        }
    }
    // Add the temporary CPU vectors to the particle structure
    np = particle_z.size();
    AddNParticles(0,np,
                  particle_x.dataPtr(),  particle_y.dataPtr(),  particle_z.dataPtr(),
                  particle_ux.dataPtr(), particle_uy.dataPtr(), particle_uz.dataPtr(),
                  1, particle_w.dataPtr(),1);
}

void
PhysicalParticleContainer::AddPlasmaFromFile(ParticleReal q_tot,
                                             ParticleReal z_shift)
{
    // Declare temporary vectors on the CPU
    Gpu::HostVector<ParticleReal> particle_x;
    Gpu::HostVector<ParticleReal> particle_z;
    Gpu::HostVector<ParticleReal> particle_ux;
    Gpu::HostVector<ParticleReal> particle_uz;
    Gpu::HostVector<ParticleReal> particle_w;
    Gpu::HostVector<ParticleReal> particle_y;
    Gpu::HostVector<ParticleReal> particle_uy;

#ifdef WARPX_USE_OPENPMD
    //TODO: Make changes for read/write in multiple MPI ranks
    if (ParallelDescriptor::IOProcessor()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(plasma_injector,
                                         "AddPlasmaFromFile: plasma injector not initialized.\n");
        // take ownership of the series and close it when done
        auto series = std::move(plasma_injector->m_openpmd_input_series);

        // assumption asserts: see PlasmaInjector
        openPMD::Iteration it = series->iterations.begin()->second;
        std::string const ps_name = it.particles.begin()->first;
        openPMD::ParticleSpecies ps = it.particles.begin()->second;

        auto const npart = ps["position"]["x"].getExtent()[0];
        std::shared_ptr<ParticleReal> ptr_x = ps["position"]["x"].loadChunk<ParticleReal>();
        double const position_unit_x = ps["position"]["x"].unitSI();
        std::shared_ptr<ParticleReal> ptr_z = ps["position"]["z"].loadChunk<ParticleReal>();
        double const position_unit_z = ps["position"]["z"].unitSI();
        std::shared_ptr<ParticleReal> ptr_ux = ps["momentum"]["x"].loadChunk<ParticleReal>();
        double const momentum_unit_x = ps["momentum"]["x"].unitSI();
        std::shared_ptr<ParticleReal> ptr_uz = ps["momentum"]["z"].loadChunk<ParticleReal>();
        double const momentum_unit_z = ps["momentum"]["z"].unitSI();
#   ifndef WARPX_DIM_XZ
        std::shared_ptr<ParticleReal> ptr_y = ps["position"]["y"].loadChunk<ParticleReal>();
        double const position_unit_y = ps["position"]["y"].unitSI();
#   endif
        std::shared_ptr<ParticleReal> ptr_uy = nullptr;
        double momentum_unit_y = 1.0;
        if (ps["momentum"].contains("y")) {
            ptr_uy = ps["momentum"]["y"].loadChunk<ParticleReal>();
             momentum_unit_y = ps["momentum"]["y"].unitSI();
        }
        series->flush();  // shared_ptr data can be read now

        ParticleReal weight = 1.0_prt;  // base standard: no info means "real" particles
        if (q_tot != 0.0) {
            weight = std::abs(q_tot) / ( std::abs(charge) * ParticleReal(npart) );
            if (ps.contains("weighting")) {
                Print() << "WARNING: Both '" << ps_name << ".q_tot' and '"
                        << ps_name << ".injection_file' specify a total charge.\n'"
                        << ps_name << ".q_tot' will take precedence.\n";
            }
        }
        // ED-PIC extension?
        else if (ps.contains("weighting")) {
            // TODO: Add ASSERT_WITH_MESSAGE to test if weighting is a constant record
            // TODO: Add ASSERT_WITH_MESSAGE for macroWeighted value in ED-PIC
            ParticleReal w = ps["weighting"][openPMD::RecordComponent::SCALAR].loadChunk<ParticleReal>().get()[0];
            double const w_unit = ps["weighting"][openPMD::RecordComponent::SCALAR].unitSI();
            weight = w * w_unit;
        }

        for (auto i = decltype(npart){0}; i<npart; ++i){
            ParticleReal const x = ptr_x.get()[i]*position_unit_x;
            ParticleReal const z = ptr_z.get()[i]*position_unit_z+z_shift;
#   if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
            ParticleReal const y = ptr_y.get()[i]*position_unit_y;
#   else
            ParticleReal const y = 0.0_prt;
#   endif
            if (plasma_injector->insideBounds(x, y, z)) {
                ParticleReal const ux = ptr_ux.get()[i]*momentum_unit_x/PhysConst::m_e;
                ParticleReal const uz = ptr_uz.get()[i]*momentum_unit_z/PhysConst::m_e;
                ParticleReal uy = 0.0_prt;
                if (ps["momentum"].contains("y")) {
                    uy = ptr_uy.get()[i]*momentum_unit_y/PhysConst::m_e;
                }
                CheckAndAddParticle(x, y, z, { ux, uy, uz}, weight,
                                    particle_x,  particle_y,  particle_z,
                                    particle_ux, particle_uy, particle_uz,
                                    particle_w);
            }
        }
        auto const np = particle_z.size();
        if (np < npart) {
            Print() << "WARNING: Simulation box doesn't cover all particles\n";
        }
    } // IO Processor
    auto const np = particle_z.size();
    AddNParticles(0, np,
                  particle_x.dataPtr(),  particle_y.dataPtr(),  particle_z.dataPtr(),
                  particle_ux.dataPtr(), particle_uy.dataPtr(), particle_uz.dataPtr(),
                  1, particle_w.dataPtr(),1);
#endif // WARPX_USE_OPENPMD

    ignore_unused(q_tot, z_shift);

    return;
}

void
PhysicalParticleContainer::CheckAndAddParticle (
    Real x, Real y, Real z,
    std::array<Real, 3> u,
    Real weight,
    Gpu::HostVector<ParticleReal>& particle_x,
    Gpu::HostVector<ParticleReal>& particle_y,
    Gpu::HostVector<ParticleReal>& particle_z,
    Gpu::HostVector<ParticleReal>& particle_ux,
    Gpu::HostVector<ParticleReal>& particle_uy,
    Gpu::HostVector<ParticleReal>& particle_uz,
    Gpu::HostVector<ParticleReal>& particle_w)
{
    if (WarpX::gamma_boost > 1.) {
        MapParticletoBoostedFrame(x, y, z, u);
    }
    particle_x.push_back(x);
    particle_y.push_back(y);
    particle_z.push_back(z);
    particle_ux.push_back(u[0]);
    particle_uy.push_back(u[1]);
    particle_uz.push_back(u[2]);
    particle_w.push_back(weight);
}

void
PhysicalParticleContainer::AddParticles (int lev)
{
    WARPX_PROFILE("PhysicalParticleContainer::AddParticles()");

    if (plasma_injector->add_single_particle) {
        AddNParticles(lev, 1,
                      &(plasma_injector->single_particle_pos[0]),
                      &(plasma_injector->single_particle_pos[1]),
                      &(plasma_injector->single_particle_pos[2]),
                      &(plasma_injector->single_particle_vel[0]),
                      &(plasma_injector->single_particle_vel[1]),
                      &(plasma_injector->single_particle_vel[2]),
                      1, &(plasma_injector->single_particle_weight), 0);
        return;
    }

    if (plasma_injector->gaussian_beam) {
        AddGaussianBeam(plasma_injector->x_m,
                        plasma_injector->y_m,
                        plasma_injector->z_m,
                        plasma_injector->x_rms,
                        plasma_injector->y_rms,
                        plasma_injector->z_rms,
                        plasma_injector->x_cut,
                        plasma_injector->y_cut,
                        plasma_injector->z_cut,
                        plasma_injector->q_tot,
                        plasma_injector->npart,
                        plasma_injector->do_symmetrize);


        return;
    }

    if (plasma_injector->external_file) {
        AddPlasmaFromFile(plasma_injector->q_tot,
                          plasma_injector->z_shift);
        return;
    }

    if ( plasma_injector->doInjection() ) {
        AddPlasma( lev );
    }
}

void
PhysicalParticleContainer::AddPlasma (int lev, RealBox part_realbox)
{
    WARPX_PROFILE("PhysicalParticleContainer::AddPlasma()");

    // If no part_realbox is provided, initialize particles in the whole domain
    const Geometry& geom = Geom(lev);
    if (!part_realbox.ok()) part_realbox = geom.ProbDomain();

    int num_ppc = plasma_injector->num_particles_per_cell;
#ifdef WARPX_DIM_RZ
    Real rmax = std::min(plasma_injector->xmax, part_realbox.hi(0));
#endif

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    Real scale_fac;
#if AMREX_SPACEDIM==3
    scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
#elif AMREX_SPACEDIM==2
    scale_fac = dx[0]*dx[1]/num_ppc;
#endif

    defineAllParticleTiles();

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    const int nlevs = numLevels();
    static bool refine_injection = false;
    static Box fine_injection_box;
    static int rrfac = 1;
    // This does not work if the mesh is dynamic.  But in that case, we should
    // not use refined injected either.  We also assume there is only one fine level.
    if (WarpX::do_moving_window and WarpX::refine_plasma
        and do_continuous_injection and nlevs == 2)
    {
        refine_injection = true;
        fine_injection_box = ParticleBoxArray(1).minimalBox();
        fine_injection_box.setSmall(WarpX::moving_window_dir, std::numeric_limits<int>::lowest());
        fine_injection_box.setBig(WarpX::moving_window_dir, std::numeric_limits<int>::max());
        rrfac = m_gdb->refRatio(0)[0];
        fine_injection_box.coarsen(rrfac);
    }

    InjectorPosition* inj_pos = plasma_injector->getInjectorPosition();
    InjectorDensity*  inj_rho = plasma_injector->getInjectorDensity();
    InjectorMomentum* inj_mom = plasma_injector->getInjectorMomentum();
    Real gamma_boost = WarpX::gamma_boost;
    Real beta_boost = WarpX::beta_boost;
    Real t = WarpX::GetInstance().gett_new(lev);
    Real density_min = plasma_injector->density_min;
    Real density_max = plasma_injector->density_max;

#ifdef WARPX_DIM_RZ
    const long nmodes = WarpX::n_rz_azimuthal_modes;
    bool radially_weighted = plasma_injector->radially_weighted;
#endif

    MFItInfo info;
    if (do_tiling && Gpu::notInLaunchRegion()) {
        info.EnableTiling(tile_size);
    }
#ifdef _OPENMP
    info.SetDynamic(true);
#pragma omp parallel if (not WarpX::serialize_ics)
#endif
    for (MFIter mfi = MakeMFIter(lev, info); mfi.isValid(); ++mfi)
    {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        Real wt = amrex::second();

        const Box& tile_box = mfi.tilebox();
        const RealBox tile_realbox = WarpX::getRealBox(tile_box, lev);

        // Find the cells of part_box that overlap with tile_realbox
        // If there is no overlap, just go to the next tile in the loop
        RealBox overlap_realbox;
        Box overlap_box;
        IntVect shifted;
        bool no_overlap = false;

        for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
            if ( tile_realbox.lo(dir) <= part_realbox.hi(dir) ) {
                Real ncells_adjust = std::floor( (tile_realbox.lo(dir) - part_realbox.lo(dir))/dx[dir] );
                overlap_realbox.setLo( dir, part_realbox.lo(dir) + std::max(ncells_adjust, 0._rt) * dx[dir]);
            } else {
                no_overlap = true; break;
            }
            if ( tile_realbox.hi(dir) >= part_realbox.lo(dir) ) {
                Real ncells_adjust = std::floor( (part_realbox.hi(dir) - tile_realbox.hi(dir))/dx[dir] );
                overlap_realbox.setHi( dir, part_realbox.hi(dir) - std::max(ncells_adjust, 0._rt) * dx[dir]);
            } else {
                no_overlap = true; break;
            }
            // Count the number of cells in this direction in overlap_realbox
            overlap_box.setSmall( dir, 0 );
            overlap_box.setBig( dir,
                int( std::round((overlap_realbox.hi(dir)-overlap_realbox.lo(dir))
                                /dx[dir] )) - 1);
            shifted[dir] =
                static_cast<int>(std::round((overlap_realbox.lo(dir)-problo[dir])/dx[dir]));
            // shifted is exact in non-moving-window direction.  That's all we care.
        }
        if (no_overlap == 1) {
            continue; // Go to the next tile
        }

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        const GpuArray<Real,AMREX_SPACEDIM> overlap_corner
            {AMREX_D_DECL(overlap_realbox.lo(0),
                          overlap_realbox.lo(1),
                          overlap_realbox.lo(2))};

        // count the number of particles that each cell in overlap_box could add
        Gpu::DeviceVector<int> counts(overlap_box.numPts(), 0);
        Gpu::DeviceVector<int> offset(overlap_box.numPts());
        auto pcounts = counts.data();
        int lrrfac = rrfac;
        int lrefine_injection = refine_injection;
        Box lfine_box = fine_injection_box;
        amrex::ParallelFor(overlap_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv(AMREX_D_DECL(i, j, k));
            auto lo = getCellCoords(overlap_corner, dx, {0._rt, 0._rt, 0._rt}, iv);
            auto hi = getCellCoords(overlap_corner, dx, {1._rt, 1._rt, 1._rt}, iv);

            lo.z = applyBallisticCorrection(lo, inj_mom, gamma_boost, beta_boost, t);
            hi.z = applyBallisticCorrection(hi, inj_mom, gamma_boost, beta_boost, t);

            if (inj_pos->overlapsWith(lo, hi))
            {
                auto index = overlap_box.index(iv);
                if (lrefine_injection) {
                    Box fine_overlap_box = overlap_box & amrex::shift(lfine_box, shifted);
                    if (fine_overlap_box.ok()) {
                        int r = (fine_overlap_box.contains(iv)) ?
                            AMREX_D_TERM(lrrfac,*lrrfac,*lrrfac) : 1;
                        pcounts[index] = num_ppc*r;
                    }
                } else {
                    pcounts[index] = num_ppc;
                }
            }
#if (AMREX_SPACEDIM != 3)
            amrex::ignore_unused(k);
#endif
        });

        // Max number of new particles. All of them are created,
        // and invalid ones are then discarded
        int max_new_particles = Scan::ExclusiveSum(counts.size(), counts.data(), offset.data());

        // Update NextID to include particles created in this function
        Long pid;
#ifdef _OPENMP
#pragma omp critical (add_plasma_nextid)
#endif
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+max_new_particles);
        }
        WarpXUtilMsg::AlwaysAssert(static_cast<Long>(pid + max_new_particles) < LastParticleID,
                                   "ERROR: overflow on particle id numbers");

        const int cpuid = ParallelDescriptor::MyProc();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ) {
            DefineAndReturnParticleTile(lev, grid_id, tile_id);
        }

        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + max_new_particles;
        particle_tile.resize(new_size);

        ParticleType* pp = particle_tile.GetArrayOfStructs()().data() + old_size;
        auto& soa = particle_tile.GetStructOfArrays();
        GpuArray<ParticleReal*,PIdx::nattribs> pa;
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            pa[ia] = soa.GetRealData(ia).data() + old_size;
        }

        int* pi = nullptr;
        if (do_field_ionization) {
            pi = soa.GetIntData(particle_icomps["ionization_level"]).data() + old_size;
        }

#ifdef WARPX_QED
        //Pointer to the optical depth component
        amrex::Real* p_optical_depth_QSR = nullptr;
        amrex::Real* p_optical_depth_BW  = nullptr;

        // If a QED effect is enabled, the corresponding optical depth
        // has to be initialized
        bool loc_has_quantum_sync = has_quantum_sync();
        bool loc_has_breit_wheeler = has_breit_wheeler();
        if (loc_has_quantum_sync)
            p_optical_depth_QSR = soa.GetRealData(
                particle_comps["optical_depth_QSR"]).data() + old_size;
        if(loc_has_breit_wheeler)
            p_optical_depth_BW = soa.GetRealData(
                particle_comps["optical_depth_BW"]).data() + old_size;

        //If needed, get the appropriate functors from the engines
        QuantumSynchrotronGetOpticalDepth quantum_sync_get_opt;
        BreitWheelerGetOpticalDepth breit_wheeler_get_opt;
        if(loc_has_quantum_sync){
            quantum_sync_get_opt =
                m_shr_p_qs_engine->build_optical_depth_functor();
        }
        if(loc_has_breit_wheeler){
            breit_wheeler_get_opt =
                m_shr_p_bw_engine->build_optical_depth_functor();
        }
#endif

        bool loc_do_field_ionization = do_field_ionization;
        int loc_ionization_initial_level = ionization_initial_level;

        // Loop over all new particles and inject them (creates too many
        // particles, in particular does not consider xmin, xmax etc.).
        // The invalid ones are given negative ID and are deleted during the
        // next redistribute.
        const auto poffset = offset.data();
        amrex::ParallelForRNG(overlap_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            IntVect iv = IntVect(AMREX_D_DECL(i, j, k));
            const auto index = overlap_box.index(iv);
            for (int i_part = 0; i_part < pcounts[index]; ++i_part)
            {
                long ip = poffset[index] + i_part;
                ParticleType& p = pp[ip];
                p.id() = pid+ip;
                p.cpu() = cpuid;

                const XDim3 r =
                    inj_pos->getPositionUnitBox(i_part, lrrfac, engine);
                auto pos = getCellCoords(overlap_corner, dx, r, iv);

#if (AMREX_SPACEDIM == 3)
                if (!tile_realbox.contains(XDim3{pos.x,pos.y,pos.z})) {
                    p.id() = -1;
                    continue;
                }
#else
                amrex::ignore_unused(k);
                if (!tile_realbox.contains(XDim3{pos.x,pos.z,0.0_rt})) {
                    p.id() = -1;
                    continue;
                }
#endif

                // Save the x and y values to use in the insideBounds checks.
                // This is needed with WARPX_DIM_RZ since x and y are modified.
                Real xb = pos.x;
                Real yb = pos.y;

#ifdef WARPX_DIM_RZ
                // Replace the x and y, setting an angle theta.
                // These x and y are used to get the momentum and density
                Real theta;
                if (nmodes == 1) {
                    // With only 1 mode, the angle doesn't matter so
                    // choose it randomly.
                    theta = 2._rt*MathConst::pi*amrex::Random(engine);
                } else {
                    theta = 2._rt*MathConst::pi*r.y;
                }
                pos.x = xb*std::cos(theta);
                pos.y = xb*std::sin(theta);
#endif

                Real dens;
                XDim3 u;
                if (gamma_boost == 1._rt) {
                    // Lab-frame simulation
                    // If the particle is not within the species's
                    // xmin, xmax, ymin, ymax, zmin, zmax, go to
                    // the next generated particle.

                    // include ballistic correction for plasma species with bulk motion
                    const Real z0 = applyBallisticCorrection(pos, inj_mom, gamma_boost,
                                                             beta_boost, t);
                    if (!inj_pos->insideBounds(xb, yb, z0)) {
                        p.id() = -1;
                        continue;
                    }

                    u = inj_mom->getMomentum(pos.x, pos.y, z0, engine);
                    dens = inj_rho->getDensity(pos.x, pos.y, z0);

                    // Remove particle if density below threshold
                    if ( dens < density_min ){
                        p.id() = -1;
                        continue;
                    }
                    // Cut density if above threshold
                    dens = amrex::min(dens, density_max);
                } else {
                    // Boosted-frame simulation
                    const Real z0_lab = applyBallisticCorrection(pos, inj_mom, gamma_boost,
                                                                 beta_boost, t);

                    // If the particle is not within the lab-frame zmin, zmax, etc.
                    // go to the next generated particle.
                    if (!inj_pos->insideBounds(xb, yb, z0_lab)) {
                        p.id() = -1;
                        continue;
                    }
                    // call `getDensity` with lab-frame parameters
                    dens = inj_rho->getDensity(pos.x, pos.y, z0_lab);
                    // Remove particle if density below threshold
                    if ( dens < density_min ){
                        p.id() = -1;
                        continue;
                    }
                    // Cut density if above threshold
                    dens = amrex::min(dens, density_max);

                    // get the full momentum, including thermal motion
                    u = inj_mom->getMomentum(pos.x, pos.y, 0._rt, engine);
                    const Real gamma_lab = std::sqrt( 1._rt+(u.x*u.x+u.y*u.y+u.z*u.z) );
                    const Real betaz_lab = u.z/(gamma_lab);

                    // At this point u and dens are the lab-frame quantities
                    // => Perform Lorentz transform
                    dens = gamma_boost * dens * ( 1.0_rt - beta_boost*betaz_lab );
                    u.z = gamma_boost * ( u.z -beta_boost*gamma_lab );
                }

                if (loc_do_field_ionization) {
                    pi[ip] = loc_ionization_initial_level;
                }

#ifdef WARPX_QED
                if(loc_has_quantum_sync){
                    p_optical_depth_QSR[ip] = quantum_sync_get_opt();
                }

                if(loc_has_breit_wheeler){
                    p_optical_depth_BW[ip] = breit_wheeler_get_opt(engine);
                }
#endif

                u.x *= PhysConst::c;
                u.y *= PhysConst::c;
                u.z *= PhysConst::c;

                // Real weight = dens * scale_fac / (AMREX_D_TERM(fac, *fac, *fac));
                Real weight = dens * scale_fac;
#ifdef WARPX_DIM_RZ
                if (radially_weighted) {
                    weight *= 2._rt*MathConst::pi*xb;
                } else {
                    // This is not correct since it might shift the particle
                    // out of the local grid
                    pos.x = std::sqrt(xb*rmax);
                    weight *= dx[0];
                }
#endif
                pa[PIdx::w ][ip] = weight;
                pa[PIdx::ux][ip] = u.x;
                pa[PIdx::uy][ip] = u.y;
                pa[PIdx::uz][ip] = u.z;

#if (AMREX_SPACEDIM == 3)
                p.pos(0) = pos.x;
                p.pos(1) = pos.y;
                p.pos(2) = pos.z;
#elif (AMREX_SPACEDIM == 2)
#ifdef WARPX_DIM_RZ
                pa[PIdx::theta][ip] = theta;
#endif
                p.pos(0) = xb;
                p.pos(1) = pos.z;
#endif
            }
        });

        amrex::Gpu::synchronize();

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }

    // The function that calls this is responsible for redistributing particles.
}

void
PhysicalParticleContainer::Evolve (int lev,
                                   const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                   const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                   const MultiFab& Ex_avg, const MultiFab& Ey_avg, const MultiFab& Ez_avg,
                                   const MultiFab& Bx_avg, const MultiFab& By_avg, const MultiFab& Bz_avg,
                                   MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                   MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                   MultiFab* rho, MultiFab* crho,
                                   const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                   const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                   Real /*t*/, Real dt, DtType a_dt_type)
{

    WARPX_PROFILE("PhysicalParticleContainer::Evolve()");
    WARPX_PROFILE_VAR_NS("PhysicalParticleContainer::Evolve::GatherAndPush", blp_fg);

    BL_ASSERT(OnSameGrids(lev,jx));

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    const iMultiFab* current_masks = WarpX::CurrentBufferMasks(lev);
    const iMultiFab* gather_masks = WarpX::GatherBufferMasks(lev);

    bool has_buffer = cEx || cjx;

    if (WarpX::do_back_transformed_diagnostics && do_back_transformed_diagnostics)
    {
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const auto np = pti.numParticles();
            const auto t_lev = pti.GetLevel();
            const auto index = pti.GetPairIndex();
            tmp_particle_data.resize(finestLevel()+1);
            for (int i = 0; i < TmpIdx::nattribs; ++i)
                tmp_particle_data[t_lev][index][i].resize(np);
        }
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif

        FArrayBox filtered_Ex, filtered_Ey, filtered_Ez;
        FArrayBox filtered_Bx, filtered_By, filtered_Bz;

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            Real wt = amrex::second();

            const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto&  wp = attribs[PIdx::w];
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];

            const long np = pti.numParticles();

            // Data on the grid
            FArrayBox const* exfab = WarpX::fft_do_time_averaging ? &(Ex_avg[pti]) : &(Ex[pti]);
            FArrayBox const* eyfab = WarpX::fft_do_time_averaging ? &(Ey_avg[pti]) : &(Ey[pti]);
            FArrayBox const* ezfab = WarpX::fft_do_time_averaging ? &(Ez_avg[pti]) : &(Ez[pti]);
            FArrayBox const* bxfab = WarpX::fft_do_time_averaging ? &(Bx_avg[pti]) : &(Bx[pti]);
            FArrayBox const* byfab = WarpX::fft_do_time_averaging ? &(By_avg[pti]) : &(By[pti]);
            FArrayBox const* bzfab = WarpX::fft_do_time_averaging ? &(Bz_avg[pti]) : &(Bz[pti]);

            Elixir exeli, eyeli, ezeli, bxeli, byeli, bzeli;

            if (WarpX::use_fdtd_nci_corr)
            {
                // Filter arrays Ex[pti], store the result in
                // filtered_Ex and update pointer exfab so that it
                // points to filtered_Ex (and do the same for all
                // components of E and B).
                applyNCIFilter(lev, pti.tilebox(), exeli, eyeli, ezeli, bxeli, byeli, bzeli,
                               filtered_Ex, filtered_Ey, filtered_Ez,
                               filtered_Bx, filtered_By, filtered_Bz,
                               Ex[pti], Ey[pti], Ez[pti], Bx[pti], By[pti], Bz[pti],
                               exfab, eyfab, ezfab, bxfab, byfab, bzfab);
            }

            // Determine which particles deposit/gather in the buffer, and
            // which particles deposit/gather in the fine patch
            long nfine_current = np;
            long nfine_gather = np;
            if (has_buffer && !do_not_push) {
                // - Modify `nfine_current` and `nfine_gather` (in place)
                //    so that they correspond to the number of particles
                //    that deposit/gather in the fine patch respectively.
                // - Reorder the particle arrays,
                //    so that the `nfine_current`/`nfine_gather` first particles
                //    deposit/gather in the fine patch
                //    and (thus) the `np-nfine_current`/`np-nfine_gather` last particles
                //    deposit/gather in the buffer
                PartitionParticlesInBuffers( nfine_current, nfine_gather, np,
                    pti, lev, current_masks, gather_masks, uxp, uyp, uzp, wp );
            }

            const long np_current = (cjx) ? nfine_current : np;

            if (rho) {
                // Deposit charge before particle push, in component 0 of MultiFab rho.
                int* AMREX_RESTRICT ion_lev;
                if (do_field_ionization){
                    ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
                } else {
                    ion_lev = nullptr;
                }
                DepositCharge(pti, wp, ion_lev, rho, 0, 0,
                              np_current, thread_num, lev, lev);
                if (has_buffer){
                    DepositCharge(pti, wp, ion_lev, crho, 0, np_current,
                                  np-np_current, thread_num, lev, lev-1);
                }
            }

            if (! do_not_push)
            {
                const long np_gather = (cEx) ? nfine_gather : np;

                int e_is_nodal = Ex.is_nodal() and Ey.is_nodal() and Ez.is_nodal();

                //
                // Gather and push for particles not in the buffer
                //
                WARPX_PROFILE_VAR_START(blp_fg);
                PushPX(pti, exfab, eyfab, ezfab,
                       bxfab, byfab, bzfab,
                       Ex.nGrow(), e_is_nodal,
                       0, np_gather, lev, lev, dt, ScaleFields(false), a_dt_type);

                if (np_gather < np)
                {
                    const IntVect& ref_ratio = WarpX::RefRatio(lev-1);
                    const Box& cbox = amrex::coarsen(box,ref_ratio);

                    // Data on the grid
                    FArrayBox const* cexfab = &(*cEx)[pti];
                    FArrayBox const* ceyfab = &(*cEy)[pti];
                    FArrayBox const* cezfab = &(*cEz)[pti];
                    FArrayBox const* cbxfab = &(*cBx)[pti];
                    FArrayBox const* cbyfab = &(*cBy)[pti];
                    FArrayBox const* cbzfab = &(*cBz)[pti];

                    if (WarpX::use_fdtd_nci_corr)
                    {
                        // Filter arrays (*cEx)[pti], store the result in
                        // filtered_Ex and update pointer cexfab so that it
                        // points to filtered_Ex (and do the same for all
                        // components of E and B)
                        applyNCIFilter(lev-1, cbox, exeli, eyeli, ezeli, bxeli, byeli, bzeli,
                                       filtered_Ex, filtered_Ey, filtered_Ez,
                                       filtered_Bx, filtered_By, filtered_Bz,
                                       (*cEx)[pti], (*cEy)[pti], (*cEz)[pti],
                                       (*cBx)[pti], (*cBy)[pti], (*cBz)[pti],
                                       cexfab, ceyfab, cezfab, cbxfab, cbyfab, cbzfab);
                    }

                    // Field gather and push for particles in gather buffers
                    e_is_nodal = cEx->is_nodal() and cEy->is_nodal() and cEz->is_nodal();
                    PushPX(pti, cexfab, ceyfab, cezfab,
                           cbxfab, cbyfab, cbzfab,
                           cEx->nGrow(), e_is_nodal,
                           nfine_gather, np-nfine_gather,
                           lev, lev-1, dt, ScaleFields(false), a_dt_type);
                }

                WARPX_PROFILE_VAR_STOP(blp_fg);

                //
                // Current Deposition (only needed for electromagnetic solver)
                //
                if (!WarpX::do_electrostatic) {
                    int* AMREX_RESTRICT ion_lev;
                    if (do_field_ionization){
                        ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
                    } else {
                        ion_lev = nullptr;
                    }
                    // Deposit inside domains
                    DepositCurrent(pti, wp, uxp, uyp, uzp, ion_lev, &jx, &jy, &jz,
                                   0, np_current, thread_num,
                                   lev, lev, dt);
                    if (has_buffer){
                        // Deposit in buffers
                        DepositCurrent(pti, wp, uxp, uyp, uzp, ion_lev, cjx, cjy, cjz,
                                       np_current, np-np_current, thread_num,
                                       lev, lev-1, dt);
                    }
                } // end of "if !do_electrostatic"
            } // end of "if do_not_push"

            if (rho) {
                // Deposit charge after particle push, in component 1 of MultiFab rho.
                // (Skipped for electrostatic solver, as this may lead to out-of-bounds)
                if (!WarpX::do_electrostatic) {
                    int* AMREX_RESTRICT ion_lev;
                    if (do_field_ionization){
                        ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
                    } else {
                        ion_lev = nullptr;
                    }
                    DepositCharge(pti, wp, ion_lev, rho, 1, 0,
                                  np_current, thread_num, lev, lev);
                    if (has_buffer){
                        DepositCharge(pti, wp, ion_lev, crho, 1, np_current,
                                      np-np_current, thread_num, lev, lev-1);
                    }
                }
            }

            amrex::Gpu::synchronize();

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add( &(*cost)[pti.index()], wt);
            }
        }
    }
    // Split particles at the end of the timestep.
    // When subcycling is ON, the splitting is done on the last call to
    // PhysicalParticleContainer::Evolve on the finest level, i.e., at the
    // end of the large timestep. Otherwise, the pushes on different levels
    // are not consistent, and the call to Redistribute (inside
    // SplitParticles) may result in split particles to deposit twice on the
    // coarse level.
    if (do_splitting && (a_dt_type == DtType::SecondHalf || a_dt_type == DtType::Full) ){
        SplitParticles(lev);
    }
}

void
PhysicalParticleContainer::applyNCIFilter (
    int lev, const Box& box,
    Elixir& exeli, Elixir& eyeli, Elixir& ezeli,
    Elixir& bxeli, Elixir& byeli, Elixir& bzeli,
    FArrayBox& filtered_Ex, FArrayBox& filtered_Ey, FArrayBox& filtered_Ez,
    FArrayBox& filtered_Bx, FArrayBox& filtered_By, FArrayBox& filtered_Bz,
    const FArrayBox& Ex, const FArrayBox& Ey, const FArrayBox& Ez,
    const FArrayBox& Bx, const FArrayBox& By, const FArrayBox& Bz,
    FArrayBox const * & ex_ptr, FArrayBox const * & ey_ptr,
    FArrayBox const * & ez_ptr, FArrayBox const * & bx_ptr,
    FArrayBox const * & by_ptr, FArrayBox const * & bz_ptr)
{

    // Get instances of NCI Godfrey filters
    const auto& nci_godfrey_filter_exeybz = WarpX::GetInstance().nci_godfrey_filter_exeybz;
    const auto& nci_godfrey_filter_bxbyez = WarpX::GetInstance().nci_godfrey_filter_bxbyez;

#if (AMREX_SPACEDIM == 2)
    const Box& tbox = amrex::grow(box,{static_cast<int>(WarpX::nox),
                static_cast<int>(WarpX::noz)});
#else
    const Box& tbox = amrex::grow(box,{static_cast<int>(WarpX::nox),
                static_cast<int>(WarpX::noy),
                static_cast<int>(WarpX::noz)});
#endif

    // Filter Ex (Both 2D and 3D)
    filtered_Ex.resize(amrex::convert(tbox,Ex.box().ixType()));
    // Safeguard for GPU
    exeli = filtered_Ex.elixir();
    // Apply filter on Ex, result stored in filtered_Ex

    nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Ex, Ex, filtered_Ex.box());
    // Update ex_ptr reference
    ex_ptr = &filtered_Ex;

    // Filter Ez
    filtered_Ez.resize(amrex::convert(tbox,Ez.box().ixType()));
    ezeli = filtered_Ez.elixir();
    nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_Ez, Ez, filtered_Ez.box());
    ez_ptr = &filtered_Ez;

    // Filter By
    filtered_By.resize(amrex::convert(tbox,By.box().ixType()));
    byeli = filtered_By.elixir();
    nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_By, By, filtered_By.box());
    by_ptr = &filtered_By;
#if (AMREX_SPACEDIM == 3)
    // Filter Ey
    filtered_Ey.resize(amrex::convert(tbox,Ey.box().ixType()));
    eyeli = filtered_Ey.elixir();
    nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Ey, Ey, filtered_Ey.box());
    ey_ptr = &filtered_Ey;

    // Filter Bx
    filtered_Bx.resize(amrex::convert(tbox,Bx.box().ixType()));
    bxeli = filtered_Bx.elixir();
    nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_Bx, Bx, filtered_Bx.box());
    bx_ptr = &filtered_Bx;

    // Filter Bz
    filtered_Bz.resize(amrex::convert(tbox,Bz.box().ixType()));
    bzeli = filtered_Bz.elixir();
    nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Bz, Bz, filtered_Bz.box());
    bz_ptr = &filtered_Bz;
#else
    amrex::ignore_unused(eyeli, bxeli, bzeli,
        filtered_Ey, filtered_Bx, filtered_Bz,
        Ey, Bx, Bz, ey_ptr, bx_ptr, bz_ptr);
#endif
}

// Loop over all particles in the particle container and
// split particles tagged with p.id()=DoSplitParticleID
void
PhysicalParticleContainer::SplitParticles (int lev)
{
    auto& mypc = WarpX::GetInstance().GetPartContainer();
    auto& pctmp_split = mypc.GetPCtmp();
    RealVector psplit_x, psplit_y, psplit_z, psplit_w;
    RealVector psplit_ux, psplit_uy, psplit_uz;
    long np_split_to_add = 0;
    long np_split;
    if(split_type==0)
    {
        np_split = (AMREX_SPACEDIM == 3) ? 8 : 4;
    } else {
        np_split = 2*AMREX_SPACEDIM;
    }

    // Loop over particle interator
    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const auto GetPosition = GetParticlePosition(pti);

        const amrex::Vector<int> ppc_nd = plasma_injector->num_particles_per_cell_each_dim;
        const std::array<Real,3>& dx = WarpX::CellSize(lev);
        amrex::Vector<Real> split_offset = {dx[0]/2._rt,
                                            dx[1]/2._rt,
                                            dx[2]/2._rt};
        if (ppc_nd[0] > 0){
            // offset for split particles is computed as a function of cell size
            // and number of particles per cell, so that a uniform distribution
            // before splitting results in a uniform distribution after splitting
            split_offset[0] /= ppc_nd[0];
            split_offset[1] /= ppc_nd[1];
            split_offset[2] /= ppc_nd[2];
        }
        // particle Array Of Structs data
        auto& particles = pti.GetArrayOfStructs();
        // particle Struct Of Arrays data
        auto& attribs = pti.GetAttribs();
        auto& wp  = attribs[PIdx::w ];
        auto& uxp = attribs[PIdx::ux];
        auto& uyp = attribs[PIdx::uy];
        auto& uzp = attribs[PIdx::uz];
        const long np = pti.numParticles();
        for(int i=0; i<np; i++){
            ParticleReal xp, yp, zp;
            GetPosition(i, xp, yp, zp);
            auto& p = particles[i];
            if (p.id() == DoSplitParticleID){
                // If particle is tagged, split it and put the
                // split particles in local arrays psplit_x etc.
                np_split_to_add += np_split;
#if (AMREX_SPACEDIM==2)
                if (split_type==0){
                    // Split particle in two along each diagonals
                    // 4 particles in 2d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        for (int kshift = -1; kshift < 2; kshift +=2 ){
                            // Add one particle with offset in x and z
                            psplit_x.push_back( xp + ishift*split_offset[0] );
                            psplit_y.push_back( yp );
                            psplit_z.push_back( zp + kshift*split_offset[2] );
                            psplit_ux.push_back( uxp[i] );
                            psplit_uy.push_back( uyp[i] );
                            psplit_uz.push_back( uzp[i] );
                            psplit_w.push_back( wp[i]/np_split );
                        }
                    }
                } else {
                    // Split particle in two along each axis
                    // 4 particles in 2d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        // Add one particle with offset in x
                        psplit_x.push_back( xp + ishift*split_offset[0] );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in z
                        psplit_x.push_back( xp );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp + ishift*split_offset[2] );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                    }
                }
#elif (AMREX_SPACEDIM==3)
                if (split_type==0){
                    // Split particle in two along each diagonals
                    // 8 particles in 3d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        for (int jshift = -1; jshift < 2; jshift +=2 ){
                            for (int kshift = -1; kshift < 2; kshift +=2 ){
                                // Add one particle with offset in x, y and z
                                psplit_x.push_back( xp + ishift*split_offset[0] );
                                psplit_y.push_back( yp + jshift*split_offset[1] );
                                psplit_z.push_back( zp + kshift*split_offset[2] );
                                psplit_ux.push_back( uxp[i] );
                                psplit_uy.push_back( uyp[i] );
                                psplit_uz.push_back( uzp[i] );
                                psplit_w.push_back( wp[i]/np_split );
                            }
                        }
                    }
                } else {
                    // Split particle in two along each axis
                    // 6 particles in 3d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        // Add one particle with offset in x
                        psplit_x.push_back( xp + ishift*split_offset[0] );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in y
                        psplit_x.push_back( xp );
                        psplit_y.push_back( yp + ishift*split_offset[1] );
                        psplit_z.push_back( zp );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in z
                        psplit_x.push_back( xp );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp + ishift*split_offset[2] );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                    }
                }
#endif
                // invalidate the particle
                p.id() = -p.id();
            }
        }
    }
    // Add local arrays psplit_x etc. to the temporary
    // particle container pctmp_split. Split particles
    // are tagged with p.id()=NoSplitParticleID so that
    // they are not re-split when entering a higher level
    // AddNParticles calls Redistribute, so that particles
    // in pctmp_split are in the proper grids and tiles
    pctmp_split.AddNParticles(lev,
                              np_split_to_add,
                              psplit_x.dataPtr(),
                              psplit_y.dataPtr(),
                              psplit_z.dataPtr(),
                              psplit_ux.dataPtr(),
                              psplit_uy.dataPtr(),
                              psplit_uz.dataPtr(),
                              1,
                              psplit_w.dataPtr(),
                              1, NoSplitParticleID);
    // Copy particles from tmp to current particle container
    addParticles(pctmp_split,1);
    // Clear tmp container
    pctmp_split.clearParticles();
}

void
PhysicalParticleContainer::PushP (int lev, Real dt,
                                  const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                  const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    WARPX_PROFILE("PhysicalParticleContainer::PushP()");

    if (do_not_push) return;

    const std::array<amrex::Real,3>& dx = WarpX::CellSize(std::max(lev,0));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            amrex::Box box = pti.tilebox();
            box.grow(Ex.nGrow());

            const long np = pti.numParticles();

            // Data on the grid
            const FArrayBox& exfab = Ex[pti];
            const FArrayBox& eyfab = Ey[pti];
            const FArrayBox& ezfab = Ez[pti];
            const FArrayBox& bxfab = Bx[pti];
            const FArrayBox& byfab = By[pti];
            const FArrayBox& bzfab = Bz[pti];

            const auto getPosition = GetParticlePosition(pti);

            const auto getExternalE = GetExternalEField(pti);
            const auto getExternalB = GetExternalBField(pti);

            const auto& xyzmin = WarpX::GetInstance().LowerCornerWithGalilean(box,m_v_galilean,lev);

            const Dim3 lo = lbound(box);

            bool galerkin_interpolation = WarpX::galerkin_interpolation;
            int nox = WarpX::nox;
            int n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;

            amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
            amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};

            amrex::Array4<const amrex::Real> const& ex_arr = exfab.array();
            amrex::Array4<const amrex::Real> const& ey_arr = eyfab.array();
            amrex::Array4<const amrex::Real> const& ez_arr = ezfab.array();
            amrex::Array4<const amrex::Real> const& bx_arr = bxfab.array();
            amrex::Array4<const amrex::Real> const& by_arr = byfab.array();
            amrex::Array4<const amrex::Real> const& bz_arr = bzfab.array();

            amrex::IndexType const ex_type = exfab.box().ixType();
            amrex::IndexType const ey_type = eyfab.box().ixType();
            amrex::IndexType const ez_type = ezfab.box().ixType();
            amrex::IndexType const bx_type = bxfab.box().ixType();
            amrex::IndexType const by_type = byfab.box().ixType();
            amrex::IndexType const bz_type = bzfab.box().ixType();

            auto& attribs = pti.GetAttribs();
            ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

            int* AMREX_RESTRICT ion_lev = nullptr;
            if (do_field_ionization) {
                ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
            }

            // Loop over the particles and update their momentum
            const amrex::Real q = this->charge;
            const amrex::Real m = this-> mass;

            const auto pusher_algo = WarpX::particle_pusher_algo;
            const auto do_crr = do_classical_radiation_reaction;

            const auto t_do_not_gather = do_not_gather;

            amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
            {
                amrex::ParticleReal xp, yp, zp;
                getPosition(ip, xp, yp, zp);

                amrex::ParticleReal Exp = 0._rt, Eyp = 0._rt, Ezp = 0._rt;
                amrex::ParticleReal Bxp = 0._rt, Byp = 0._rt, Bzp = 0._rt;

                if (!t_do_not_gather){
                    // first gather E and B to the particle positions
                    doGatherShapeN(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                   ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                   ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                   dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes,
                                   nox, galerkin_interpolation);
                }
                // Externally applied E-field in Cartesian co-ordinates
                getExternalE(ip, Exp, Eyp, Ezp);
                // Externally applied B-field in Cartesian co-ordinates
                getExternalB(ip, Bxp, Byp, Bzp);

                if (do_crr) {
                    amrex::Real qp = q;
                    if (ion_lev) { qp *= ion_lev[ip]; }
                    UpdateMomentumBorisWithRadiationReaction(ux[ip], uy[ip], uz[ip],
                                                             Exp, Eyp, Ezp, Bxp,
                                                             Byp, Bzp, qp, m, dt);
                } else if (pusher_algo == ParticlePusherAlgo::Boris) {
                    amrex::Real qp = q;
                    if (ion_lev) { qp *= ion_lev[ip]; }
                    UpdateMomentumBoris( ux[ip], uy[ip], uz[ip],
                                         Exp, Eyp, Ezp, Bxp,
                                         Byp, Bzp, qp, m, dt);
                } else if (pusher_algo == ParticlePusherAlgo::Vay) {
                    amrex::Real qp = q;
                    if (ion_lev){ qp *= ion_lev[ip]; }
                    UpdateMomentumVay( ux[ip], uy[ip], uz[ip],
                                       Exp, Eyp, Ezp, Bxp,
                                       Byp, Bzp, qp, m, dt);
                } else if (pusher_algo == ParticlePusherAlgo::HigueraCary) {
                    amrex::Real qp = q;
                    if (ion_lev){ qp *= ion_lev[ip]; }
                    UpdateMomentumHigueraCary( ux[ip], uy[ip], uz[ip],
                                               Exp, Eyp, Ezp, Bxp,
                                               Byp, Bzp, qp, m, dt);
                } else {
                    amrex::Abort("Unknown particle pusher");
                }
            });
        }
    }
}

void
PhysicalParticleContainer::GetParticleSlice (
    const int direction, const Real z_old,
    const Real z_new, const Real t_boost,
    const Real t_lab, const Real dt,
    DiagnosticParticles& diagnostic_particles)
{
    WARPX_PROFILE("PhysicalParticleContainer::GetParticleSlice()");

    // Assume that the boost in the positive z direction.
#if (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT(direction == 1);
#else
    AMREX_ALWAYS_ASSERT(direction == 2);
#endif

    // Note the the slice should always move in the negative boost direction.
    AMREX_ALWAYS_ASSERT(z_new < z_old);

    AMREX_ALWAYS_ASSERT(do_back_transformed_diagnostics == 1);

    const int nlevs = std::max(0, finestLevel()+1);

    // we figure out a box for coarse-grained rejection. If the RealBox corresponding to a
    // given tile doesn't intersect with this, there is no need to check any particles.
    const Real* base_dx = Geom(0).CellSize();
    const Real z_min = z_new - base_dx[direction];
    const Real z_max = z_old + base_dx[direction];

    RealBox slice_box = Geom(0).ProbDomain();
    slice_box.setLo(direction, z_min);
    slice_box.setHi(direction, z_max);

    diagnostic_particles.resize(finestLevel()+1);

    for (int lev = 0; lev < nlevs; ++lev) {

        const Real* dx  = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();

        // first we touch each map entry in serial
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
            diagnostic_particles[lev][index];
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            // Temporary arrays to store copy_flag and copy_index
            // for particles that cross the z-slice
            // These arrays are defined before the WarpXParIter to prevent them
            // from going out of scope after each iteration, while the kernels
            // may still need access to them.
            // Note that the destructor for WarpXParIter is synchronized.
            amrex::Gpu::DeviceVector<int> FlagForPartCopy;
            amrex::Gpu::DeviceVector<int> IndexForPartCopy;
            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
                const Box& box = pti.validbox();
                auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                const RealBox tile_real_box(box, dx, plo);

                if ( !slice_box.intersects(tile_real_box) ) continue;

                const auto GetPosition = GetParticlePosition(pti);

                auto& attribs = pti.GetAttribs();
                Real* const AMREX_RESTRICT wpnew = attribs[PIdx::w].dataPtr();
                Real* const AMREX_RESTRICT uxpnew = attribs[PIdx::ux].dataPtr();
                Real* const AMREX_RESTRICT uypnew = attribs[PIdx::uy].dataPtr();
                Real* const AMREX_RESTRICT uzpnew = attribs[PIdx::uz].dataPtr();

                Real* const AMREX_RESTRICT
                  xpold = tmp_particle_data[lev][index][TmpIdx::xold].dataPtr();
                Real* const AMREX_RESTRICT
                  ypold = tmp_particle_data[lev][index][TmpIdx::yold].dataPtr();
                Real* const AMREX_RESTRICT
                  zpold = tmp_particle_data[lev][index][TmpIdx::zold].dataPtr();
                Real* const AMREX_RESTRICT
                  uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
                Real* const AMREX_RESTRICT
                  uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
                Real* const AMREX_RESTRICT
                  uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();

                const long np = pti.numParticles();

                Real uzfrm = -WarpX::gamma_boost*WarpX::beta_boost*PhysConst::c;
                Real inv_c2 = 1.0/PhysConst::c/PhysConst::c;

                FlagForPartCopy.resize(np);
                IndexForPartCopy.resize(np);

                int* const AMREX_RESTRICT Flag = FlagForPartCopy.dataPtr();
                int* const AMREX_RESTRICT IndexLocation = IndexForPartCopy.dataPtr();

                //Flag particles that need to be copied if they cross the z_slice
                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE(int i)
                {
                    ParticleReal xp, yp, zp;
                    GetPosition(i, xp, yp, zp);
                    Flag[i] = 0;
                    if ( (((zp >= z_new) && (zpold[i] <= z_old)) ||
                          ((zp <= z_new) && (zpold[i] >= z_old))) )
                    {
                        Flag[i] = 1;
                    }
                });

                // exclusive scan to obtain location indices using flag values
                // These location indices are used to copy data from
                // src to dst when the copy-flag is set to 1.
                const int total_partdiag_size = amrex::Scan::ExclusiveSum(np,Flag,IndexLocation);

                // allocate array size for diagnostic particle array
                diagnostic_particles[lev][index].resize(total_partdiag_size);

                amrex::Real gammaboost = WarpX::gamma_boost;
                amrex::Real betaboost = WarpX::beta_boost;
                amrex::Real Phys_c = PhysConst::c;

                Real* const AMREX_RESTRICT diag_wp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::w).data();
                Real* const AMREX_RESTRICT diag_xp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::x).data();
                Real* const AMREX_RESTRICT diag_yp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::y).data();
                Real* const AMREX_RESTRICT diag_zp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::z).data();
                Real* const AMREX_RESTRICT diag_uxp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::ux).data();
                Real* const AMREX_RESTRICT diag_uyp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::uy).data();
                Real* const AMREX_RESTRICT diag_uzp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::uz).data();

                // Copy particle data to diagnostic particle array on the GPU
                //  using flag and index values
                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE(int i)
                {
                    ParticleReal xp_new, yp_new, zp_new;
                    GetPosition(i, xp_new, yp_new, zp_new);
                    if (Flag[i] == 1)
                    {
                         // Lorentz Transform particles to lab-frame
                         const Real gamma_new_p = std::sqrt(1.0 + inv_c2*
                                                  (uxpnew[i]*uxpnew[i]
                                                 + uypnew[i]*uypnew[i]
                                                 + uzpnew[i]*uzpnew[i]));
                         const Real t_new_p = gammaboost*t_boost - uzfrm*zp_new*inv_c2;
                         const Real z_new_p = gammaboost*(zp_new + betaboost*Phys_c*t_boost);
                         const Real uz_new_p = gammaboost*uzpnew[i] - gamma_new_p*uzfrm;

                         const Real gamma_old_p = std::sqrt(1.0 + inv_c2*
                                                  (uxpold[i]*uxpold[i]
                                                 + uypold[i]*uypold[i]
                                                 + uzpold[i]*uzpold[i]));
                         const Real t_old_p = gammaboost*(t_boost - dt)
                                              - uzfrm*zpold[i]*inv_c2;
                         const Real z_old_p = gammaboost*(zpold[i]
                                              + betaboost*Phys_c*(t_boost-dt));
                         const Real uz_old_p = gammaboost*uzpold[i]
                                              - gamma_old_p*uzfrm;

                         // interpolate in time to t_lab
                         const Real weight_old = (t_new_p - t_lab)
                                               / (t_new_p - t_old_p);
                         const Real weight_new = (t_lab - t_old_p)
                                               / (t_new_p - t_old_p);

                         const Real xp = xpold[i]*weight_old + xp_new*weight_new;
                         const Real yp = ypold[i]*weight_old + yp_new*weight_new;
                         const Real zp = z_old_p*weight_old  + z_new_p*weight_new;

                         const Real uxp = uxpold[i]*weight_old
                                        + uxpnew[i]*weight_new;
                         const Real uyp = uypold[i]*weight_old
                                        + uypnew[i]*weight_new;
                         const Real uzp = uz_old_p*weight_old
                                        + uz_new_p  *weight_new;

                         const int loc = IndexLocation[i];
                         diag_wp[loc] = wpnew[i];
                         diag_xp[loc] = xp;
                         diag_yp[loc] = yp;
                         diag_zp[loc] = zp;
                         diag_uxp[loc] = uxp;
                         diag_uyp[loc] = uyp;
                         diag_uzp[loc] = uzp;
                    }
                });
                Gpu::synchronize(); // because of FlagForPartCopy & IndexForPartCopy
            }
        }
    }
}

/* \brief Inject particles during the simulation
 * \param injection_box: domain where particles should be injected.
 */
void
PhysicalParticleContainer::ContinuousInjection (const RealBox& injection_box)
{
    // Inject plasma on level 0. Paticles will be redistributed.
    const int lev=0;
    AddPlasma(lev, injection_box);
}

/* \brief Perform the field gather and particle push operations in one fused kernel
 *
 */
void
PhysicalParticleContainer::PushPX (WarpXParIter& pti,
                                   amrex::FArrayBox const * exfab,
                                   amrex::FArrayBox const * eyfab,
                                   amrex::FArrayBox const * ezfab,
                                   amrex::FArrayBox const * bxfab,
                                   amrex::FArrayBox const * byfab,
                                   amrex::FArrayBox const * bzfab,
                                   const int ngE, const int /*e_is_nodal*/,
                                   const long offset,
                                   const long np_to_push,
                                   int lev, int gather_lev,
                                   amrex::Real dt, ScaleFields scaleFields,
                                   DtType a_dt_type)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((gather_lev==(lev-1)) ||
                                     (gather_lev==(lev  )),
                                     "Gather buffers only work for lev-1");
    // If no particles, do not do anything
    if (np_to_push == 0) return;

    // Get cell size on gather_lev
    const std::array<Real,3>& dx = WarpX::CellSize(std::max(gather_lev,0));

    // Get box from which field is gathered.
    // If not gathering from the finest level, the box is coarsened.
    Box box;
    if (lev == gather_lev) {
        box = pti.tilebox();
    } else {
        const IntVect& ref_ratio = WarpX::RefRatio(gather_lev);
        box = amrex::coarsen(pti.tilebox(),ref_ratio);
    }

    // Add guard cells to the box.
    box.grow(ngE);

    const auto getPosition = GetParticlePosition(pti, offset);
          auto setPosition = SetParticlePosition(pti, offset);

    const auto getExternalE = GetExternalEField(pti, offset);
    const auto getExternalB = GetExternalBField(pti, offset);

    // Lower corner of tile box physical domain (take into account Galilean shift)
    Real cur_time = WarpX::GetInstance().gett_new(lev);
    const auto& time_of_last_gal_shift = WarpX::GetInstance().time_of_last_gal_shift;
    Real time_shift = (cur_time - time_of_last_gal_shift);
    amrex::Array<amrex::Real,3> galilean_shift ={
        m_v_galilean[0]*time_shift,
        m_v_galilean[1]*time_shift,
        m_v_galilean[2]*time_shift };
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(box, galilean_shift, gather_lev);

    const Dim3 lo = lbound(box);

    bool galerkin_interpolation = WarpX::galerkin_interpolation;
    int nox = WarpX::nox;
    int n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;

    amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
    amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};

    amrex::Array4<const amrex::Real> const& ex_arr = exfab->array();
    amrex::Array4<const amrex::Real> const& ey_arr = eyfab->array();
    amrex::Array4<const amrex::Real> const& ez_arr = ezfab->array();
    amrex::Array4<const amrex::Real> const& bx_arr = bxfab->array();
    amrex::Array4<const amrex::Real> const& by_arr = byfab->array();
    amrex::Array4<const amrex::Real> const& bz_arr = bzfab->array();

    amrex::IndexType const ex_type = exfab->box().ixType();
    amrex::IndexType const ey_type = eyfab->box().ixType();
    amrex::IndexType const ez_type = ezfab->box().ixType();
    amrex::IndexType const bx_type = bxfab->box().ixType();
    amrex::IndexType const by_type = byfab->box().ixType();
    amrex::IndexType const bz_type = bzfab->box().ixType();

    auto& attribs = pti.GetAttribs();
    ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    auto copyAttribs = CopyParticleAttribs(pti, tmp_particle_data, offset);
    int do_copy = (WarpX::do_back_transformed_diagnostics &&
                          do_back_transformed_diagnostics &&
                   (a_dt_type!=DtType::SecondHalf));

    int* AMREX_RESTRICT ion_lev = nullptr;
    if (do_field_ionization) {
        ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
    }

    // Loop over the particles and update their momentum
    const amrex::Real q = this->charge;
    const amrex::Real m = this-> mass;

    const auto pusher_algo = WarpX::particle_pusher_algo;
    const auto do_crr = do_classical_radiation_reaction;
#ifdef WARPX_QED
    const auto do_sync = m_do_qed_quantum_sync;
    amrex::Real t_chi_max = 0.0;
    if (do_sync) t_chi_max = m_shr_p_qs_engine->get_minimum_chi_part();

    QuantumSynchrotronEvolveOpticalDepth evolve_opt;
    amrex::ParticleReal* AMREX_RESTRICT p_optical_depth_QSR = nullptr;
    const bool local_has_quantum_sync = has_quantum_sync();
    if (local_has_quantum_sync) {
        evolve_opt = m_shr_p_qs_engine->build_evolve_functor();
        p_optical_depth_QSR = pti.GetAttribs(particle_comps["optical_depth_QSR"]).dataPtr();
    }
#endif

    const auto t_do_not_gather = do_not_gather;

    amrex::ParallelFor( np_to_push, [=] AMREX_GPU_DEVICE (long ip)
    {
        amrex::ParticleReal xp, yp, zp;
        getPosition(ip, xp, yp, zp);

        amrex::ParticleReal Exp = 0._rt, Eyp = 0._rt, Ezp = 0._rt;
        amrex::ParticleReal Bxp = 0._rt, Byp = 0._rt, Bzp = 0._rt;

        if(!t_do_not_gather){
            // first gather E and B to the particle positions
            doGatherShapeN(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                           ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                           ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                           dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes,
                           nox, galerkin_interpolation);
        }
        // Externally applied E-field in Cartesian co-ordinates
        getExternalE(ip, Exp, Eyp, Ezp);
        // Externally applied B-field in Cartesian co-ordinates
        getExternalB(ip, Bxp, Byp, Bzp);

        scaleFields(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp);

        doParticlePush(getPosition, setPosition, copyAttribs, ip,
                       ux[ip+offset], uy[ip+offset], uz[ip+offset],
                       Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                       ion_lev ? ion_lev[ip] : 0,
                       m, q, pusher_algo, do_crr, do_copy,
#ifdef WARPX_QED
                       do_sync,
                       t_chi_max,
#endif
                       dt);

#ifdef WARPX_QED
    if (local_has_quantum_sync) {
        evolve_opt(ux[ip], uy[ip], uz[ip],
                   Exp, Eyp, Ezp,Bxp, Byp, Bzp,
                   dt, p_optical_depth_QSR[ip]);
    }
#endif

    });
}

void
PhysicalParticleContainer::InitIonizationModule ()
{
    if (!do_field_ionization) return;
    ParmParse pp(species_name);
    if (charge != PhysConst::q_e){
        amrex::Warning(
            "charge != q_e for ionizable species: overriding user value and setting charge = q_e.");
        charge = PhysConst::q_e;
    }
    pp.query("ionization_initial_level", ionization_initial_level);
    pp.get("ionization_product_species", ionization_product_name);
    pp.get("physical_element", physical_element);
    // Add runtime integer component for ionization level
    AddIntComp("ionization_level");
    // Get atomic number and ionization energies from file
    int ion_element_id = ion_map_ids[physical_element];
    ion_atomic_number = ion_atomic_numbers[ion_element_id];
    Vector<Real> h_ionization_energies(ion_atomic_number);
    int offset = ion_energy_offsets[ion_element_id];
    for(int i=0; i<ion_atomic_number; i++){
        h_ionization_energies[i] = table_ionization_energies[i+offset];
    }
    // Compute ADK prefactors (See Chen, JCP 236 (2013), equation (2))
    // For now, we assume l=0 and m=0.
    // The approximate expressions are used,
    // without Gamma function
    Real wa = std::pow(PhysConst::alpha,3) * PhysConst::c / PhysConst::r_e;
    Real Ea = PhysConst::m_e * PhysConst::c*PhysConst::c /PhysConst::q_e *
        std::pow(PhysConst::alpha,4)/PhysConst::r_e;
    Real UH = table_ionization_energies[0];
    Real l_eff = std::sqrt(UH/h_ionization_energies[0]) - 1.;

    const Real dt = WarpX::GetInstance().getdt(0);

    ionization_energies.resize(ion_atomic_number);
    adk_power.resize(ion_atomic_number);
    adk_prefactor.resize(ion_atomic_number);
    adk_exp_prefactor.resize(ion_atomic_number);

    Gpu::copyAsync(Gpu::hostToDevice,
                   h_ionization_energies.begin(), h_ionization_energies.end(),
                   ionization_energies.begin());

    Real const* AMREX_RESTRICT p_ionization_energies = ionization_energies.data();
    Real * AMREX_RESTRICT p_adk_power = adk_power.data();
    Real * AMREX_RESTRICT p_adk_prefactor = adk_prefactor.data();
    Real * AMREX_RESTRICT p_adk_exp_prefactor = adk_exp_prefactor.data();
    amrex::ParallelFor(ion_atomic_number, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        Real n_eff = (i+1) * std::sqrt(UH/p_ionization_energies[i]);
        Real C2 = std::pow(2,2*n_eff)/(n_eff*tgamma(n_eff+l_eff+1)*tgamma(n_eff-l_eff));
        p_adk_power[i] = -(2*n_eff - 1);
        Real Uion = p_ionization_energies[i];
        p_adk_prefactor[i] = dt * wa * C2 * ( Uion/(2*UH) )
            * std::pow(2*std::pow((Uion/UH),3./2)*Ea,2*n_eff - 1);
        p_adk_exp_prefactor[i] = -2./3 * std::pow( Uion/UH,3./2) * Ea;
    });

    Gpu::synchronize();
}

IonizationFilterFunc
PhysicalParticleContainer::getIonizationFunc (const WarpXParIter& pti,
                                              int lev,
                                              int ngE,
                                              const amrex::FArrayBox& Ex,
                                              const amrex::FArrayBox& Ey,
                                              const amrex::FArrayBox& Ez,
                                              const amrex::FArrayBox& Bx,
                                              const amrex::FArrayBox& By,
                                              const amrex::FArrayBox& Bz)
{
    WARPX_PROFILE("PhysicalParticleContainer::getIonizationFunc()");

    return IonizationFilterFunc(pti, lev, ngE, Ex, Ey, Ez, Bx, By, Bz,
                                m_v_galilean,
                                ionization_energies.dataPtr(),
                                adk_prefactor.dataPtr(),
                                adk_exp_prefactor.dataPtr(),
                                adk_power.dataPtr(),
                                particle_icomps["ionization_level"],
                                ion_atomic_number);
}

void PhysicalParticleContainer::resample (const int timestep)
{
    // In heavily load imbalanced simulations, MPI processes with few particles will spend most of
    // the time at the MPI synchronization in TotalNumberOfParticles(). Having two profiler entries
    // here is thus useful to avoid confusing time spent waiting for other processes with time
    // spent doing actual resampling.
    WARPX_PROFILE_VAR_NS("MultiParticleContainer::doResampling::MPI_synchronization",
                         blp_resample_synchronization);
    WARPX_PROFILE_VAR_NS("MultiParticleContainer::doResampling::ActualResampling",
                         blp_resample_actual);

    WARPX_PROFILE_VAR_START(blp_resample_synchronization);
    const amrex::Real global_numparts = TotalNumberOfParticles();
    WARPX_PROFILE_VAR_STOP(blp_resample_synchronization);

    WARPX_PROFILE_VAR_START(blp_resample_actual);
    if (m_resampler.triggered(timestep, global_numparts))
    {
        amrex::Print() << "Resampling " << species_name << ".\n";
        for (int lev = 0; lev <= maxLevel(); lev++)
        {
            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
                m_resampler(pti, lev, this);
            }
        }
    }
    WARPX_PROFILE_VAR_STOP(blp_resample_actual);

}


#ifdef WARPX_QED


bool PhysicalParticleContainer::has_quantum_sync () const
{
    return m_do_qed_quantum_sync;
}

bool PhysicalParticleContainer::has_breit_wheeler () const
{
    return m_do_qed_breit_wheeler;
}

void
PhysicalParticleContainer::
set_breit_wheeler_engine_ptr (std::shared_ptr<BreitWheelerEngine> ptr)
{
    m_shr_p_bw_engine = ptr;
}

void
PhysicalParticleContainer::
set_quantum_sync_engine_ptr (std::shared_ptr<QuantumSynchrotronEngine> ptr)
{
    m_shr_p_qs_engine = ptr;
}

PhotonEmissionFilterFunc
PhysicalParticleContainer::getPhotonEmissionFilterFunc ()
{
    WARPX_PROFILE("PhysicalParticleContainer::getPhotonEmissionFunc()");
    return PhotonEmissionFilterFunc{particle_runtime_comps["optical_depth_QSR"]};
}

PairGenerationFilterFunc
PhysicalParticleContainer::getPairGenerationFilterFunc ()
{
    WARPX_PROFILE("PhysicalParticleContainer::getPairGenerationFunc()");
    return PairGenerationFilterFunc{particle_runtime_comps["optical_depth_BW"]};
}

#endif
