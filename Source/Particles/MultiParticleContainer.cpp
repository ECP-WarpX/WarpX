/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Jean-Luc Vay, Luca Fedeli
 * Mathieu Lobet, Maxence Thevenet, Neil Zaim
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MultiParticleContainer.H"
#include "Particles/ElementaryProcess/Ionization.H"
#ifdef WARPX_QED
#   include "Particles/ElementaryProcess/QEDInternals/BreitWheelerEngineWrapper.H"
#   include "Particles/ElementaryProcess/QEDInternals/QuantumSyncEngineWrapper.H"
#   include "Particles/ElementaryProcess/QEDSchwingerProcess.H"
#   include "Particles/ElementaryProcess/QEDPairGeneration.H"
#   include "Particles/ElementaryProcess/QEDPhotonEmission.H"
#endif
#include "Particles/LaserParticleContainer.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#ifdef WARPX_QED
#   include "Particles/ParticleCreation/FilterCreateTransformFromFAB.H"
#endif
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartCreate.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/PhotonParticleContainer.H"
#include "Particles/PhysicalParticleContainer.H"
#include "Particles/RigidInjectedParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "SpeciesPhysicalProperties.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXProfilerWrapper.H"
#ifdef AMREX_USE_EB
#   include "EmbeddedBoundary/ParticleScraper.H"
#   include "EmbeddedBoundary/ParticleBoundaryProcess.H"
#endif

#include "WarpX.H"

#include <ablastr/utils/Communication.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace amrex;

namespace
{
    /** A little collection to transport six Array4 that point to the EM fields */
    struct MyFieldList
    {
        Array4< amrex::Real const > const Ex, Ey, Ez, Bx, By, Bz;
    };
}

MultiParticleContainer::MultiParticleContainer (AmrCore* amr_core)
{

    ReadParameters();

    auto const nspecies = static_cast<int>(species_names.size());
    auto const nlasers = static_cast<int>(lasers_names.size());

    allcontainers.resize(nspecies + nlasers);
    for (int i = 0; i < nspecies; ++i) {
        if (species_types[i] == PCTypes::Physical) {
            allcontainers[i] = std::make_unique<PhysicalParticleContainer>(amr_core, i, species_names[i]);
        }
        else if (species_types[i] == PCTypes::RigidInjected) {
            allcontainers[i] = std::make_unique<RigidInjectedParticleContainer>(amr_core, i, species_names[i]);
        }
        else if (species_types[i] == PCTypes::Photon) {
            allcontainers[i] = std::make_unique<PhotonParticleContainer>(amr_core, i, species_names[i]);
        }
        allcontainers[i]->m_deposit_on_main_grid = m_deposit_on_main_grid[i];
        allcontainers[i]->m_gather_from_main_grid = m_gather_from_main_grid[i];
    }

    for (int i = nspecies; i < nspecies+nlasers; ++i) {
        allcontainers[i] = std::make_unique<LaserParticleContainer>(amr_core, i, lasers_names[i-nspecies]);
    }

    pc_tmp = std::make_unique<PhysicalParticleContainer>(amr_core);

    // Compute the number of species for which lab-frame data is dumped
    // nspecies_lab_frame_diags, and map their ID to MultiParticleContainer
    // particle IDs in map_species_lab_diags.
    map_species_back_transformed_diagnostics.resize(nspecies);
    nspecies_back_transformed_diagnostics = 0;
    for (int i=0; i<nspecies; i++){
        auto& pc = allcontainers[i];
        if (pc->do_back_transformed_diagnostics){
            map_species_back_transformed_diagnostics[nspecies_back_transformed_diagnostics] = i;
            do_back_transformed_diagnostics = 1;
            nspecies_back_transformed_diagnostics += 1;
        }
    }

    // Setup particle collisions
    collisionhandler = std::make_unique<CollisionHandler>(this);

}

void
MultiParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        ParmParse pp_particles("particles");

        // allocating and initializing default values of external fields for particles
        m_E_external_particle.resize(3);
        m_B_external_particle.resize(3);
        // initialize E and B fields to 0.0
        for (int idim = 0; idim < 3; ++idim) {
            m_E_external_particle[idim] = 0.0;
            m_B_external_particle[idim] = 0.0;
        }
        // default values of E_external_particle and B_external_particle
        // are used to set the E and B field when "constant" or "parser"
        // is not explicitly used in the input
        pp_particles.query("B_ext_particle_init_style", m_B_ext_particle_s);
        std::transform(m_B_ext_particle_s.begin(),
                       m_B_ext_particle_s.end(),
                       m_B_ext_particle_s.begin(),
                       ::tolower);
        pp_particles.query("E_ext_particle_init_style", m_E_ext_particle_s);
        std::transform(m_E_ext_particle_s.begin(),
                       m_E_ext_particle_s.end(),
                       m_E_ext_particle_s.begin(),
                       ::tolower);
        // if the input string for B_external on particles is "constant"
        // then the values for the external B on particles must
        // be provided in the input file.
        if (m_B_ext_particle_s == "constant")
            getArrWithParser(pp_particles, "B_external_particle", m_B_external_particle);

        // if the input string for E_external on particles is "constant"
        // then the values for the external E on particles must
        // be provided in the input file.
        if (m_E_ext_particle_s == "constant")
            getArrWithParser(pp_particles, "E_external_particle", m_E_external_particle);

        // if the input string for B_ext_particle_s is
        // "parse_b_ext_particle_function" then the mathematical expression
        // for the Bx_, By_, Bz_external_particle_function(x,y,z)
        // must be provided in the input file.
        if (m_B_ext_particle_s == "parse_b_ext_particle_function") {
           // store the mathematical expression as string
           std::string str_Bx_ext_particle_function;
           std::string str_By_ext_particle_function;
           std::string str_Bz_ext_particle_function;
           Store_parserString(pp_particles, "Bx_external_particle_function(x,y,z,t)",
                                      str_Bx_ext_particle_function);
           Store_parserString(pp_particles, "By_external_particle_function(x,y,z,t)",
                                      str_By_ext_particle_function);
           Store_parserString(pp_particles, "Bz_external_particle_function(x,y,z,t)",
                                      str_Bz_ext_particle_function);

           // Parser for B_external on the particle
           m_Bx_particle_parser = std::make_unique<amrex::Parser>(
                                    makeParser(str_Bx_ext_particle_function,{"x","y","z","t"}));
           m_By_particle_parser = std::make_unique<amrex::Parser>(
                                    makeParser(str_By_ext_particle_function,{"x","y","z","t"}));
           m_Bz_particle_parser = std::make_unique<amrex::Parser>(
                                    makeParser(str_Bz_ext_particle_function,{"x","y","z","t"}));

        }

        // if the input string for E_ext_particle_s is
        // "parse_e_ext_particle_function" then the mathematical expression
        // for the Ex_, Ey_, Ez_external_particle_function(x,y,z)
        // must be provided in the input file.
        if (m_E_ext_particle_s == "parse_e_ext_particle_function") {
           // store the mathematical expression as string
           std::string str_Ex_ext_particle_function;
           std::string str_Ey_ext_particle_function;
           std::string str_Ez_ext_particle_function;
           Store_parserString(pp_particles, "Ex_external_particle_function(x,y,z,t)",
                                      str_Ex_ext_particle_function);
           Store_parserString(pp_particles, "Ey_external_particle_function(x,y,z,t)",
                                      str_Ey_ext_particle_function);
           Store_parserString(pp_particles, "Ez_external_particle_function(x,y,z,t)",
                                      str_Ez_ext_particle_function);
           // Parser for E_external on the particle
           m_Ex_particle_parser = std::make_unique<amrex::Parser>(
                                    makeParser(str_Ex_ext_particle_function,{"x","y","z","t"}));
           m_Ey_particle_parser = std::make_unique<amrex::Parser>(
                                    makeParser(str_Ey_ext_particle_function,{"x","y","z","t"}));
           m_Ez_particle_parser = std::make_unique<amrex::Parser>(
                                    makeParser(str_Ez_ext_particle_function,{"x","y","z","t"}));

        }

        // if the input string for E_ext_particle_s or B_ext_particle_s is
        // "repeated_plasma_lens" then the plasma lens properties
        // must be provided in the input file.
        if (m_E_ext_particle_s == "repeated_plasma_lens" ||
            m_B_ext_particle_s == "repeated_plasma_lens") {
            getWithParser(pp_particles, "repeated_plasma_lens_period", m_repeated_plasma_lens_period);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_repeated_plasma_lens_period > 0._rt,
                                             "The period of the repeated plasma lens must be greater than zero");
            getArrWithParser(pp_particles, "repeated_plasma_lens_starts", h_repeated_plasma_lens_starts);
            getArrWithParser(pp_particles, "repeated_plasma_lens_lengths", h_repeated_plasma_lens_lengths);

            int n_lenses = static_cast<int>(h_repeated_plasma_lens_starts.size());
            d_repeated_plasma_lens_starts.resize(n_lenses);
            d_repeated_plasma_lens_lengths.resize(n_lenses);
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                       h_repeated_plasma_lens_starts.begin(), h_repeated_plasma_lens_starts.end(),
                       d_repeated_plasma_lens_starts.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                       h_repeated_plasma_lens_lengths.begin(), h_repeated_plasma_lens_lengths.end(),
                       d_repeated_plasma_lens_lengths.begin());

            h_repeated_plasma_lens_strengths_E.resize(n_lenses);
            h_repeated_plasma_lens_strengths_B.resize(n_lenses);

            if (m_E_ext_particle_s == "repeated_plasma_lens") {
                getArrWithParser(pp_particles, "repeated_plasma_lens_strengths_E", h_repeated_plasma_lens_strengths_E);
            }
            if (m_B_ext_particle_s == "repeated_plasma_lens") {
                getArrWithParser(pp_particles, "repeated_plasma_lens_strengths_B", h_repeated_plasma_lens_strengths_B);
            }

            d_repeated_plasma_lens_strengths_E.resize(n_lenses);
            d_repeated_plasma_lens_strengths_B.resize(n_lenses);
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                       h_repeated_plasma_lens_strengths_E.begin(), h_repeated_plasma_lens_strengths_E.end(),
                       d_repeated_plasma_lens_strengths_E.begin());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                       h_repeated_plasma_lens_strengths_B.begin(), h_repeated_plasma_lens_strengths_B.end(),
                       d_repeated_plasma_lens_strengths_B.begin());

            amrex::Gpu::synchronize();
        }


        // particle species
        pp_particles.queryarr("species_names", species_names);
        auto const nspecies = species_names.size();

        if (nspecies > 0) {
            // Get species to deposit on main grid
            m_deposit_on_main_grid.resize(nspecies, false);
            std::vector<std::string> tmp;
            pp_particles.queryarr("deposit_on_main_grid", tmp);
            for (auto const& name : tmp) {
                auto it = std::find(species_names.begin(), species_names.end(), name);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    it != species_names.end(),
                    "species '" + name
                    + "' in particles.deposit_on_main_grid must be part of particles.species_names");
                int i = std::distance(species_names.begin(), it);
                m_deposit_on_main_grid[i] = true;
            }

            m_gather_from_main_grid.resize(nspecies, false);
            std::vector<std::string> tmp_gather;
            pp_particles.queryarr("gather_from_main_grid", tmp_gather);
            for (auto const& name : tmp_gather) {
                auto it = std::find(species_names.begin(), species_names.end(), name);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    it != species_names.end(),
                    "species '" + name
                        + "' in particles.gather_from_main_grid must be part of particles.species_names");
                int i = std::distance(species_names.begin(), it);
                m_gather_from_main_grid.at(i) = true;
            }

            species_types.resize(nspecies, PCTypes::Physical);

            // Get rigid-injected species
            std::vector<std::string> rigid_injected_species;
            pp_particles.queryarr("rigid_injected_species", rigid_injected_species);
            if (!rigid_injected_species.empty()) {
                for (auto const& name : rigid_injected_species) {
                    auto it = std::find(species_names.begin(), species_names.end(), name);
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        it != species_names.end(),
                        "species '" + name
                        + "' in particles.rigid_injected_species must be part of particles.species_names");
                    int i = std::distance(species_names.begin(), it);
                    species_types[i] = PCTypes::RigidInjected;
                }
            }
            // Get photon species
            std::vector<std::string> photon_species;
            pp_particles.queryarr("photon_species", photon_species);
            if (!photon_species.empty()) {
                for (auto const& name : photon_species) {
                    auto it = std::find(species_names.begin(), species_names.end(), name);
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        it != species_names.end(),
                        "species '" + name
                        + "' in particles.photon_species must be part of particles.species_names");
                    int i = std::distance(species_names.begin(), it);
                    species_types[i] = PCTypes::Photon;
                }
            }

        }
        pp_particles.query("use_fdtd_nci_corr", WarpX::use_fdtd_nci_corr);
#ifdef WARPX_DIM_RZ
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::use_fdtd_nci_corr==0,
                            "ERROR: use_fdtd_nci_corr is not supported in RZ");
#endif
#ifdef WARPX_DIM_1D_Z
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::use_fdtd_nci_corr==0,
                            "ERROR: use_fdtd_nci_corr is not supported in 1D");
#endif

        ParmParse pp_lasers("lasers");
        pp_lasers.queryarr("names", lasers_names);

#ifdef WARPX_QED
        ParmParse pp_warpx("warpx");
        pp_warpx.query("do_qed_schwinger", m_do_qed_schwinger);

        if (m_do_qed_schwinger) {
            ParmParse pp_qed_schwinger("qed_schwinger");
            pp_qed_schwinger.get("ele_product_species", m_qed_schwinger_ele_product_name);
            pp_qed_schwinger.get("pos_product_species", m_qed_schwinger_pos_product_name);
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            getWithParser(pp_qed_schwinger, "y_size",m_qed_schwinger_y_size);
#endif
            queryWithParser(pp_qed_schwinger, "threshold_poisson_gaussian",
                                              m_qed_schwinger_threshold_poisson_gaussian);
            queryWithParser(pp_qed_schwinger, "xmin", m_qed_schwinger_xmin);
            queryWithParser(pp_qed_schwinger, "xmax", m_qed_schwinger_xmax);
#if defined(WARPX_DIM_3D)
            queryWithParser(pp_qed_schwinger, "ymin", m_qed_schwinger_ymin);
            queryWithParser(pp_qed_schwinger, "ymax", m_qed_schwinger_ymax);
#endif
            queryWithParser(pp_qed_schwinger, "zmin", m_qed_schwinger_zmin);
            queryWithParser(pp_qed_schwinger, "zmax", m_qed_schwinger_zmax);
        }
#endif
        initialized = true;
    }
}

WarpXParticleContainer&
MultiParticleContainer::GetParticleContainerFromName (const std::string& name) const
{
    auto it = std::find(species_names.begin(), species_names.end(), name);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        it != species_names.end(),
        "unknown species name");
    int i = std::distance(species_names.begin(), it);
    return *allcontainers[i];
}

void
MultiParticleContainer::AllocData ()
{
    for (auto& pc : allcontainers) {
        pc->AllocData();
    }
    pc_tmp->AllocData();
}

void
MultiParticleContainer::InitData ()
{
    for (auto& pc : allcontainers) {
        pc->InitData();
    }
    pc_tmp->InitData();
    // For each species, get the ID of its product species.
    // This is used for ionization and pair creation processes.
    mapSpeciesProduct();

    CheckIonizationProductSpecies();

#ifdef WARPX_QED
    CheckQEDProductSpecies();
    InitQED();
#endif

}

void
MultiParticleContainer::Evolve (int lev,
                                const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                MultiFab* cjx,  MultiFab* cjy, MultiFab* cjz,
                                MultiFab* rho, MultiFab* crho,
                                const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                Real t, Real dt, DtType a_dt_type, bool skip_deposition)
{
    if (! skip_deposition) {
        jx.setVal(0.0);
        jy.setVal(0.0);
        jz.setVal(0.0);
        if (cjx) cjx->setVal(0.0);
        if (cjy) cjy->setVal(0.0);
        if (cjz) cjz->setVal(0.0);
        if (rho) rho->setVal(0.0);
        if (crho) crho->setVal(0.0);
    }
    for (auto& pc : allcontainers) {
        pc->Evolve(lev, Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, cjx, cjy, cjz,
                   rho, crho, cEx, cEy, cEz, cBx, cBy, cBz, t, dt, a_dt_type, skip_deposition);
    }
}

void
MultiParticleContainer::PushX (Real dt)
{
    for (auto& pc : allcontainers) {
        pc->PushX(dt);
    }
}

void
MultiParticleContainer::PushP (int lev, Real dt,
                               const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                               const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    for (auto& pc : allcontainers) {
        pc->PushP(lev, dt, Ex, Ey, Ez, Bx, By, Bz);
    }
}

std::unique_ptr<MultiFab>
MultiParticleContainer::GetZeroChargeDensity (const int lev)
{
    WarpX& warpx = WarpX::GetInstance();

    BoxArray ba = warpx.boxArray(lev);
    DistributionMapping dmap = warpx.DistributionMap(lev);
    const int ng_rho = warpx.get_ng_depos_rho().max();

    auto zero_rho = std::make_unique<MultiFab>(amrex::convert(ba,IntVect::TheNodeVector()),
                                               dmap,WarpX::ncomps,ng_rho);
    zero_rho->setVal(amrex::Real(0.0));
    return zero_rho;
}

void
MultiParticleContainer::DepositCurrent (
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& J,
    const amrex::Real dt, const amrex::Real relative_time)
{
    // Reset the J arrays
    for (int lev = 0; lev < J.size(); ++lev)
    {
        J[lev][0]->setVal(0.0, J[lev][0]->nGrowVect());
        J[lev][1]->setVal(0.0, J[lev][1]->nGrowVect());
        J[lev][2]->setVal(0.0, J[lev][2]->nGrowVect());
    }

    // Call the deposition kernel for each species
    for (auto& pc : allcontainers)
    {
        pc->DepositCurrent(J, dt, relative_time);
    }

#ifdef WARPX_DIM_RZ
    for (int lev = 0; lev < J.size(); ++lev)
    {
        WarpX::GetInstance().ApplyInverseVolumeScalingToCurrentDensity(J[lev][0].get(), J[lev][1].get(), J[lev][2].get(), lev);
    }
#endif
}

void
MultiParticleContainer::DepositCharge (
    amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
    const amrex::Real relative_time)
{
    // Reset the rho array
    for (int lev = 0; lev < rho.size(); ++lev)
    {
        rho[lev]->setVal(0.0, 0, WarpX::ncomps, rho[lev]->nGrowVect());
    }

    // Push the particles in time, if needed
    if (relative_time != 0.) PushX(relative_time);

    // Call the deposition kernel for each species
    for (auto& pc : allcontainers)
    {
        bool const local = true;
        bool const reset = false;
        bool const do_rz_volume_scaling = false;
        bool const interpolate_across_levels = false;
        pc->DepositCharge(rho, local, reset, do_rz_volume_scaling,
                              interpolate_across_levels);
    }

    // Push the particles back in time
    if (relative_time != 0.) PushX(-relative_time);

#ifdef WARPX_DIM_RZ
    for (int lev = 0; lev < rho.size(); ++lev)
    {
        WarpX::GetInstance().ApplyInverseVolumeScalingToChargeDensity(rho[lev].get(), lev);
    }
#endif
}

std::unique_ptr<MultiFab>
MultiParticleContainer::GetChargeDensity (int lev, bool local)
{
    if (allcontainers.empty())
    {
        std::unique_ptr<MultiFab> rho = GetZeroChargeDensity(lev);
        return rho;
    }
    else
    {
        std::unique_ptr<MultiFab> rho = allcontainers[0]->GetChargeDensity(lev, true);
        for (unsigned i = 1, n = allcontainers.size(); i < n; ++i) {
            std::unique_ptr<MultiFab> rhoi = allcontainers[i]->GetChargeDensity(lev, true);
            MultiFab::Add(*rho, *rhoi, 0, 0, rho->nComp(), rho->nGrowVect());
        }
        if (!local) {
            const Geometry& gm = allcontainers[0]->Geom(lev);
            ablastr::utils::communication::SumBoundary(*rho, WarpX::do_single_precision_comms, gm.periodicity());
        }
        return rho;
    }
}

void
MultiParticleContainer::SortParticlesByBin (amrex::IntVect bin_size)
{
    for (auto& pc : allcontainers) {
        pc->SortParticlesByBin(bin_size);
    }
}

void
MultiParticleContainer::Redistribute ()
{
    for (auto& pc : allcontainers) {
        pc->Redistribute();
    }
}

void
MultiParticleContainer::defineAllParticleTiles ()
{
    for (auto& pc : allcontainers) {
        pc->defineAllParticleTiles();
    }
}

void
MultiParticleContainer::RedistributeLocal (const int num_ghost)
{
    for (auto& pc : allcontainers) {
        pc->Redistribute(0, 0, 0, num_ghost);
    }
}

void
MultiParticleContainer::ApplyBoundaryConditions ()
{
    for (auto& pc : allcontainers) {
        pc->ApplyBoundaryConditions();
    }
}

Vector<Long>
MultiParticleContainer::GetZeroParticlesInGrid (const int lev) const
{
    WarpX& warpx = WarpX::GetInstance();
    const int num_boxes = warpx.boxArray(lev).size();
    const Vector<Long> r(num_boxes, 0);
    return r;
}

Vector<Long>
MultiParticleContainer::NumberOfParticlesInGrid (int lev) const
{
    if (allcontainers.empty())
    {
        const Vector<Long> r = GetZeroParticlesInGrid(lev);
        return r;
    }
    else
    {
        const bool only_valid=true, only_local=true;
        Vector<Long> r = allcontainers[0]->NumberOfParticlesInGrid(lev,only_valid,only_local);
        for (unsigned i = 1, n = allcontainers.size(); i < n; ++i) {
            const auto& ri = allcontainers[i]->NumberOfParticlesInGrid(lev,only_valid,only_local);
            for (unsigned j=0, m=ri.size(); j<m; ++j) {
                r[j] += ri[j];
            }
        }
        ParallelDescriptor::ReduceLongSum(r.data(),r.size());
        return r;
    }
}

void
MultiParticleContainer::Increment (MultiFab& mf, int lev)
{
    for (auto& pc : allcontainers) {
        pc->Increment(mf,lev);
    }
}

void
MultiParticleContainer::SetParticleBoxArray (int lev, BoxArray& new_ba)
{
    for (auto& pc : allcontainers) {
        pc->SetParticleBoxArray(lev,new_ba);
    }
}

void
MultiParticleContainer::SetParticleDistributionMap (int lev, DistributionMapping& new_dm)
{
    for (auto& pc : allcontainers) {
        pc->SetParticleDistributionMap(lev,new_dm);
    }
}

void
MultiParticleContainer::PostRestart ()
{
    for (auto& pc : allcontainers) {
        pc->PostRestart();
    }
    pc_tmp->PostRestart();
}

void
MultiParticleContainer
::GetLabFrameData (const std::string& /*snapshot_name*/,
                   const int /*i_lab*/, const int direction,
                   const Real z_old, const Real z_new,
                   const Real t_boost, const Real t_lab, const Real dt,
                   Vector<WarpXParticleContainer::DiagnosticParticleData>& parts) const
{

    WARPX_PROFILE("MultiParticleContainer::GetLabFrameData()");

    // Loop over particle species
    for (int i = 0; i < nspecies_back_transformed_diagnostics; ++i){
        int isp = map_species_back_transformed_diagnostics[i];
        WarpXParticleContainer* pc = allcontainers[isp].get();
        WarpXParticleContainer::DiagnosticParticles diagnostic_particles;
        pc->GetParticleSlice(direction, z_old, z_new, t_boost, t_lab, dt, diagnostic_particles);
        // Here, diagnostic_particles[lev][index] is a WarpXParticleContainer::DiagnosticParticleData
        // where "lev" is the AMR level and "index" is a [grid index][tile index] pair.

        // Loop over AMR levels
        for (int lev = 0; lev <= pc->finestLevel(); ++lev){
            // Loop over [grid index][tile index] pairs
            // and Fills parts[species number i] with particle data from all grids and
            // tiles in diagnostic_particles. parts contains particles from all
            // AMR levels indistinctly.
            for (const auto& dp : diagnostic_particles[lev]){
                // it->first is the [grid index][tile index] key
                // it->second is the corresponding
                // WarpXParticleContainer::DiagnosticParticleData value
                parts[i].GetRealData(DiagIdx::w).insert(  parts[i].GetRealData(DiagIdx::w  ).end(),
                                                          dp.second.GetRealData(DiagIdx::w  ).begin(),
                                                          dp.second.GetRealData(DiagIdx::w  ).end());

                parts[i].GetRealData(DiagIdx::x).insert(  parts[i].GetRealData(DiagIdx::x  ).end(),
                                                          dp.second.GetRealData(DiagIdx::x  ).begin(),
                                                          dp.second.GetRealData(DiagIdx::x  ).end());

                parts[i].GetRealData(DiagIdx::y).insert(  parts[i].GetRealData(DiagIdx::y  ).end(),
                                                          dp.second.GetRealData(DiagIdx::y  ).begin(),
                                                          dp.second.GetRealData(DiagIdx::y  ).end());

                parts[i].GetRealData(DiagIdx::z).insert(  parts[i].GetRealData(DiagIdx::z  ).end(),
                                                          dp.second.GetRealData(DiagIdx::z  ).begin(),
                                                          dp.second.GetRealData(DiagIdx::z  ).end());

                parts[i].GetRealData(DiagIdx::ux).insert(  parts[i].GetRealData(DiagIdx::ux).end(),
                                                           dp.second.GetRealData(DiagIdx::ux).begin(),
                                                           dp.second.GetRealData(DiagIdx::ux).end());

                parts[i].GetRealData(DiagIdx::uy).insert(  parts[i].GetRealData(DiagIdx::uy).end(),
                                                           dp.second.GetRealData(DiagIdx::uy).begin(),
                                                           dp.second.GetRealData(DiagIdx::uy).end());

                parts[i].GetRealData(DiagIdx::uz).insert(  parts[i].GetRealData(DiagIdx::uz).end(),
                                                           dp.second.GetRealData(DiagIdx::uz).begin(),
                                                           dp.second.GetRealData(DiagIdx::uz).end());
            }
        }
    }
}

/* \brief Continuous injection for particles initially outside of the domain.
 * \param injection_box: Domain where new particles should be injected.
 * Loop over all WarpXParticleContainer in MultiParticleContainer and
 * calls virtual function ContinuousInjection.
 */
void
MultiParticleContainer::ContinuousInjection (const RealBox& injection_box) const
{
    for (auto& pc : allcontainers){
        if (pc->do_continuous_injection){
            pc->ContinuousInjection(injection_box);
        }
    }
}

/* \brief Update position of continuous injection parameters.
 * \param dt: simulation time step (level 0)
 * All classes inherited from WarpXParticleContainer do not have
 * a position to update (PhysicalParticleContainer does not do anything).
 */
void
MultiParticleContainer::UpdateContinuousInjectionPosition (Real dt) const
{
    for (auto& pc : allcontainers){
        if (pc->do_continuous_injection){
            pc->UpdateContinuousInjectionPosition(dt);
        }
    }
}

int
MultiParticleContainer::doContinuousInjection () const
{
    int warpx_do_continuous_injection = 0;
    for (auto& pc : allcontainers){
        if (pc->do_continuous_injection){
            warpx_do_continuous_injection = 1;
        }
    }
    return warpx_do_continuous_injection;
}

/* \brief Continuous injection of a flux of particles
 * Loop over all WarpXParticleContainer in MultiParticleContainer and
 * calls virtual function ContinuousFluxInjection.
 */
void
MultiParticleContainer::ContinuousFluxInjection (amrex::Real t, amrex::Real dt) const
{
    for (auto& pc : allcontainers){
        pc->ContinuousFluxInjection(t, dt);
    }
}

/* \brief Get ID of product species of each species.
 * The users specifies the name of the product species,
 * this routine get its ID.
 */
void
MultiParticleContainer::mapSpeciesProduct ()
{
    for (int i=0; i < static_cast<int>(species_names.size()); i++){
        auto& pc = allcontainers[i];
        // If species pc has ionization on, find species with name
        // pc->ionization_product_name and store its ID into
        // pc->ionization_product.
        if (pc->do_field_ionization){
            const int i_product = getSpeciesID(pc->ionization_product_name);
            pc->ionization_product = i_product;
        }

#ifdef WARPX_QED
        if (pc->has_breit_wheeler()){
            const int i_product_ele = getSpeciesID(
                pc->m_qed_breit_wheeler_ele_product_name);
            pc->m_qed_breit_wheeler_ele_product = i_product_ele;

            const int i_product_pos = getSpeciesID(
                pc->m_qed_breit_wheeler_pos_product_name);
            pc->m_qed_breit_wheeler_pos_product = i_product_pos;
        }

        if(pc->has_quantum_sync()){
            const int i_product_phot = getSpeciesID(
                pc->m_qed_quantum_sync_phot_product_name);
            pc->m_qed_quantum_sync_phot_product = i_product_phot;
        }
#endif

    }

#ifdef WARPX_QED
    if (m_do_qed_schwinger) {
    m_qed_schwinger_ele_product =
        getSpeciesID(m_qed_schwinger_ele_product_name);
    m_qed_schwinger_pos_product =
        getSpeciesID(m_qed_schwinger_pos_product_name);
    }
#endif
}

/* \brief Given a species name, return its ID.
 */
int
MultiParticleContainer::getSpeciesID (std::string product_str) const
{
    auto species_and_lasers_names = GetSpeciesAndLasersNames();
    int i_product = 0;
    bool found = 0;
    // Loop over species
    for (int i=0; i < static_cast<int>(species_and_lasers_names.size()); i++){
        // If species name matches, store its ID
        // into i_product
        if (species_and_lasers_names[i] == product_str){
            found = 1;
            i_product = i;
        }
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        found != 0,
        "could not find the ID of product species '"
        + product_str + "'" + ". Wrong name?");

    return i_product;
}

void
MultiParticleContainer::SetDoBackTransformedParticles (const bool do_back_transformed_particles) {
    m_do_back_transformed_particles = do_back_transformed_particles;
}

void
MultiParticleContainer::SetDoBackTransformedParticles (std::string species_name, const bool do_back_transformed_particles) {
    auto species_names_list = GetSpeciesNames();
    bool found = 0;
    // Loop over species
    for (int i = 0; i < static_cast<int>(species_names.size()); ++i) {
        // If species name matches, set back-transformed particles parameters
        if (species_names_list[i] == species_name) {
           found = 1;
           auto& pc = allcontainers[i];
           pc->SetDoBackTransformedParticles(do_back_transformed_particles);
        }
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        found != 0,
        "ERROR: could not find the ID of product species '"
        + species_name + "'" + ". Wrong name?");
}

void
MultiParticleContainer::doFieldIonization (int lev,
                                           const MultiFab& Ex,
                                           const MultiFab& Ey,
                                           const MultiFab& Ez,
                                           const MultiFab& Bx,
                                           const MultiFab& By,
                                           const MultiFab& Bz)
{
    WARPX_PROFILE("MultiParticleContainer::doFieldIonization()");

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop over all species.
    // Ionized particles in pc_source create particles in pc_product
    for (auto& pc_source : allcontainers)
    {
        if (!pc_source->do_field_ionization){ continue; }

        auto& pc_product = allcontainers[pc_source->ionization_product];

        SmartCopyFactory copy_factory(*pc_source, *pc_product);
        auto phys_pc_ptr = static_cast<PhysicalParticleContainer*>(pc_source.get());

        auto Copy      = copy_factory.getSmartCopy();
        auto Transform = IonizationTransformFunc();

        pc_source ->defineAllParticleTiles();
        pc_product->defineAllParticleTiles();

        auto info = getMFItInfo(*pc_source, *pc_product);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*pc_source, lev, info); pti.isValid(); ++pti)
        {
            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            Real wt = amrex::second();

            auto& src_tile = pc_source ->ParticlesAt(lev, pti);
            auto& dst_tile = pc_product->ParticlesAt(lev, pti);

            auto Filter = phys_pc_ptr->getIonizationFunc(pti, lev, Ex.nGrowVect(),
                                                         Ex[pti], Ey[pti], Ez[pti],
                                                         Bx[pti], By[pti], Bz[pti]);

            const auto np_dst = dst_tile.numParticles();
            const auto num_added = filterCopyTransformParticles<1>(dst_tile, src_tile, np_dst,
                                                                   Filter, Copy, Transform);

            setNewParticleIDs(dst_tile, np_dst, num_added);

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add( &(*cost)[pti.index()], wt);
            }
        }
    }
}

void
MultiParticleContainer::doCollisions ( Real cur_time, amrex::Real dt )
{
    WARPX_PROFILE("MultiParticleContainer::doCollisions()");
    collisionhandler->doCollisions(cur_time, dt, this);
}

void MultiParticleContainer::doResampling (const int timestep)
{
    for (auto& pc : allcontainers)
    {
        // do_resampling can only be true for PhysicalParticleContainers
        if (!pc->do_resampling){ continue; }

        pc->resample(timestep);
    }
}

void MultiParticleContainer::CheckIonizationProductSpecies()
{
    for (int i=0; i < static_cast<int>(species_names.size()); i++){
        if (allcontainers[i]->do_field_ionization){
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != allcontainers[i]->ionization_product,
                "ERROR: ionization product cannot be the same species");
        }
    }
}

void MultiParticleContainer::ScrapeParticles (const amrex::Vector<const amrex::MultiFab*>& distance_to_eb)
{
#ifdef AMREX_USE_EB
    for (auto& pc : allcontainers) {
        scrapeParticles(*pc, distance_to_eb, ParticleBoundaryProcess::Absorb());
    }
#else
    amrex::ignore_unused(distance_to_eb);
#endif
}

#ifdef WARPX_QED
void MultiParticleContainer::InitQED ()
{
    m_shr_p_qs_engine = std::make_shared<QuantumSynchrotronEngine>();
    m_shr_p_bw_engine = std::make_shared<BreitWheelerEngine>();

    m_nspecies_quantum_sync = 0;
    m_nspecies_breit_wheeler = 0;

    for (auto& pc : allcontainers) {
        if(pc->has_quantum_sync()){
            pc->set_quantum_sync_engine_ptr
                (m_shr_p_qs_engine);
            m_nspecies_quantum_sync++;
        }
        if(pc->has_breit_wheeler()){
            pc->set_breit_wheeler_engine_ptr
                (m_shr_p_bw_engine);
            m_nspecies_breit_wheeler++;
        }
    }

    if(m_nspecies_quantum_sync != 0)
        InitQuantumSync();

    if(m_nspecies_breit_wheeler !=0)
        InitBreitWheeler();

}

void MultiParticleContainer::InitQuantumSync ()
{
    std::string lookup_table_mode;
    ParmParse pp_qed_qs("qed_qs");

    //If specified, use a user-defined energy threshold for photon creation
    ParticleReal temp;
    constexpr auto mec2 = PhysConst::c * PhysConst::c * PhysConst::m_e;
    if(queryWithParser(pp_qed_qs, "photon_creation_energy_threshold", temp)){
        temp *= mec2;
        m_quantum_sync_photon_creation_energy_threshold = temp;
    }
    else{
        ablastr::warn_manager::WMRecordWarning("QED",
            "Using default value (2*me*c^2) for photon energy creation threshold",
            ablastr::warn_manager::WarnPriority::low);
    }

    // qs_minimum_chi_part is the minimum chi parameter to be
    // considered for Synchrotron emission. If a lepton has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real qs_minimum_chi_part;
    getWithParser(pp_qed_qs, "chi_min", qs_minimum_chi_part);


    pp_qed_qs.query("lookup_table_mode", lookup_table_mode);
    if(lookup_table_mode.empty()){
        amrex::Abort("Quantum Synchrotron table mode should be provided");
    }

    if(lookup_table_mode == "generate"){
        ablastr::warn_manager::WMRecordWarning("QED",
            "A new Quantum Synchrotron table will be generated.",
            ablastr::warn_manager::WarnPriority::low);
#ifndef WARPX_QED_TABLE_GEN
        amrex::Error("Error: Compile with QED_TABLE_GEN=TRUE to enable table generation!\n");
#else
        QuantumSyncGenerateTable();
#endif
    }
    else if(lookup_table_mode == "load"){
        std::string load_table_name;
        pp_qed_qs.query("load_table_from", load_table_name);
        ablastr::warn_manager::WMRecordWarning("QED",
            "The Quantum Synchrotron table will be read from the file: " + load_table_name,
            ablastr::warn_manager::WarnPriority::low);
        if(load_table_name.empty()){
            amrex::Abort("Quantum Synchrotron table name should be provided");
        }
        Vector<char> table_data;
        ParallelDescriptor::ReadAndBcastFile(load_table_name, table_data);
        ParallelDescriptor::Barrier();
        m_shr_p_qs_engine->init_lookup_tables_from_raw_data(table_data,
            qs_minimum_chi_part);
    }
    else if(lookup_table_mode == "builtin"){
        ablastr::warn_manager::WMRecordWarning("QED",
            "The built-in Quantum Synchrotron table will be used. "
            "This low resolution table is intended for testing purposes only.",
            ablastr::warn_manager::WarnPriority::medium);
        m_shr_p_qs_engine->init_builtin_tables(qs_minimum_chi_part);
    }
    else{
        amrex::Abort("Unknown Quantum Synchrotron table mode");
    }

    if(!m_shr_p_qs_engine->are_lookup_tables_initialized()){
        amrex::Abort("Table initialization has failed!");
    }
}

void MultiParticleContainer::InitBreitWheeler ()
{
    std::string lookup_table_mode;
    ParmParse pp_qed_bw("qed_bw");

    // bw_minimum_chi_phot is the minimum chi parameter to be
    // considered for pair production. If a photon has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real bw_minimum_chi_part;
    if(!queryWithParser(pp_qed_bw, "chi_min", bw_minimum_chi_part))
        amrex::Abort("qed_bw.chi_min should be provided!");

    pp_qed_bw.query("lookup_table_mode", lookup_table_mode);
    if(lookup_table_mode.empty()){
        amrex::Abort("Breit Wheeler table mode should be provided");
    }

    if(lookup_table_mode == "generate"){
        ablastr::warn_manager::WMRecordWarning("QED",
            "A new Breit Wheeler table will be generated.",
            ablastr::warn_manager::WarnPriority::low);
#ifndef WARPX_QED_TABLE_GEN
        amrex::Error("Error: Compile with QED_TABLE_GEN=TRUE to enable table generation!\n");
#else
        BreitWheelerGenerateTable();
#endif
    }
    else if(lookup_table_mode == "load"){
        std::string load_table_name;
        pp_qed_bw.query("load_table_from", load_table_name);
        ablastr::warn_manager::WMRecordWarning("QED",
            "The Breit Wheeler table will be read from the file:" + load_table_name,
            ablastr::warn_manager::WarnPriority::low);
        if(load_table_name.empty()){
            amrex::Abort("Breit Wheeler table name should be provided");
        }
        Vector<char> table_data;
        ParallelDescriptor::ReadAndBcastFile(load_table_name, table_data);
        ParallelDescriptor::Barrier();
        m_shr_p_bw_engine->init_lookup_tables_from_raw_data(
            table_data, bw_minimum_chi_part);
    }
    else if(lookup_table_mode == "builtin"){
        ablastr::warn_manager::WMRecordWarning("QED",
            "The built-in Breit Wheeler table will be used. "
            "This low resolution table is intended for testing purposes only.",
            ablastr::warn_manager::WarnPriority::medium);
        m_shr_p_bw_engine->init_builtin_tables(bw_minimum_chi_part);
    }
    else{
        amrex::Abort("Unknown Breit Wheeler table mode");
    }

    if(!m_shr_p_bw_engine->are_lookup_tables_initialized()){
        amrex::Abort("Table initialization has failed!");
    }
}

void
MultiParticleContainer::QuantumSyncGenerateTable ()
{
    ParmParse pp_qed_qs("qed_qs");
    std::string table_name;
    pp_qed_qs.query("save_table_in", table_name);
    if(table_name.empty())
        amrex::Abort("qed_qs.save_table_in should be provided!");

    // qs_minimum_chi_part is the minimum chi parameter to be
    // considered for Synchrotron emission. If a lepton has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real qs_minimum_chi_part;
    getWithParser(pp_qed_qs, "chi_min", qs_minimum_chi_part);

    if(ParallelDescriptor::IOProcessor()){
        PicsarQuantumSyncCtrl ctrl;

        //==Table parameters==

        //--- sub-table 1 (1D)
        //These parameters are used to pre-compute a function
        //which appears in the evolution of the optical depth

        //Minimun chi for the table. If a lepton has chi < tab_dndt_chi_min,
        //chi is considered as if it were equal to tab_dndt_chi_min
        getWithParser(pp_qed_qs, "tab_dndt_chi_min", ctrl.dndt_params.chi_part_min);

        //Maximum chi for the table. If a lepton has chi > tab_dndt_chi_max,
        //chi is considered as if it were equal to tab_dndt_chi_max
        getWithParser(pp_qed_qs, "tab_dndt_chi_max", ctrl.dndt_params.chi_part_max);

        //How many points should be used for chi in the table
        getWithParser(pp_qed_qs, "tab_dndt_how_many", ctrl.dndt_params.chi_part_how_many);
        //------

        //--- sub-table 2 (2D)
        //These parameters are used to pre-compute a function
        //which is used to extract the properties of the generated
        //photons.

        //Minimun chi for the table. If a lepton has chi < tab_em_chi_min,
        //chi is considered as if it were equal to tab_em_chi_min
        getWithParser(pp_qed_qs, "tab_em_chi_min", ctrl.phot_em_params.chi_part_min);

        //Maximum chi for the table. If a lepton has chi > tab_em_chi_max,
        //chi is considered as if it were equal to tab_em_chi_max
        getWithParser(pp_qed_qs, "tab_em_chi_max", ctrl.phot_em_params.chi_part_max);

        //How many points should be used for chi in the table
        getWithParser(pp_qed_qs, "tab_em_chi_how_many", ctrl.phot_em_params.chi_part_how_many);

        //The other axis of the table is the ratio between the quantum
        //parameter of the emitted photon and the quantum parameter of the
        //lepton. This parameter is the minimum ratio to consider for the table.
        getWithParser(pp_qed_qs, "tab_em_frac_min", ctrl.phot_em_params.frac_min);

        //This parameter is the number of different points to consider for the second
        //axis
        getWithParser(pp_qed_qs, "tab_em_frac_how_many", ctrl.phot_em_params.frac_how_many);
        //====================

        m_shr_p_qs_engine->compute_lookup_tables(ctrl, qs_minimum_chi_part);
        const auto data = m_shr_p_qs_engine->export_lookup_tables_data();
        WarpXUtilIO::WriteBinaryDataOnFile(table_name,
            Vector<char>{data.begin(), data.end()});
    }

    ParallelDescriptor::Barrier();
    Vector<char> table_data;
    ParallelDescriptor::ReadAndBcastFile(table_name, table_data);
    ParallelDescriptor::Barrier();

    //No need to initialize from raw data for the processor that
    //has just generated the table
    if(!ParallelDescriptor::IOProcessor()){
        m_shr_p_qs_engine->init_lookup_tables_from_raw_data(
            table_data, qs_minimum_chi_part);
    }
}

void
MultiParticleContainer::BreitWheelerGenerateTable ()
{
    ParmParse pp_qed_bw("qed_bw");
    std::string table_name;
    pp_qed_bw.query("save_table_in", table_name);
    if(table_name.empty())
        amrex::Abort("qed_bw.save_table_in should be provided!");

    // bw_minimum_chi_phot is the minimum chi parameter to be
    // considered for pair production. If a photon has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real bw_minimum_chi_part;
    getWithParser(pp_qed_bw, "chi_min", bw_minimum_chi_part);

    if(ParallelDescriptor::IOProcessor()){
        PicsarBreitWheelerCtrl ctrl;

        //==Table parameters==

        //--- sub-table 1 (1D)
        //These parameters are used to pre-compute a function
        //which appears in the evolution of the optical depth

        //Minimun chi for the table. If a photon has chi < tab_dndt_chi_min,
        //an analytical approximation is used.
        getWithParser(pp_qed_bw, "tab_dndt_chi_min", ctrl.dndt_params.chi_phot_min);

        //Maximum chi for the table. If a photon has chi > tab_dndt_chi_max,
        //an analytical approximation is used.
        getWithParser(pp_qed_bw, "tab_dndt_chi_max", ctrl.dndt_params.chi_phot_max);

        //How many points should be used for chi in the table
        getWithParser(pp_qed_bw, "tab_dndt_how_many", ctrl.dndt_params.chi_phot_how_many);
        //------

        //--- sub-table 2 (2D)
        //These parameters are used to pre-compute a function
        //which is used to extract the properties of the generated
        //particles.

        //Minimun chi for the table. If a photon has chi < tab_pair_chi_min
        //chi is considered as it were equal to chi_phot_tpair_min
        getWithParser(pp_qed_bw, "tab_pair_chi_min", ctrl.pair_prod_params.chi_phot_min);

        //Maximum chi for the table. If a photon has chi > tab_pair_chi_max
        //chi is considered as it were equal to chi_phot_tpair_max
        getWithParser(pp_qed_bw, "tab_pair_chi_max", ctrl.pair_prod_params.chi_phot_max);

        //How many points should be used for chi in the table
        getWithParser(pp_qed_bw, "tab_pair_chi_how_many", ctrl.pair_prod_params.chi_phot_how_many);

        //The other axis of the table is the fraction of the initial energy
        //'taken away' by the most energetic particle of the pair.
        //This parameter is the number of different fractions to consider
        getWithParser(pp_qed_bw, "tab_pair_frac_how_many", ctrl.pair_prod_params.frac_how_many);
        //====================

        m_shr_p_bw_engine->compute_lookup_tables(ctrl, bw_minimum_chi_part);
        const auto data = m_shr_p_bw_engine->export_lookup_tables_data();
        WarpXUtilIO::WriteBinaryDataOnFile(table_name,
            Vector<char>{data.begin(), data.end()});
    }

    ParallelDescriptor::Barrier();
    Vector<char> table_data;
    ParallelDescriptor::ReadAndBcastFile(table_name, table_data);
    ParallelDescriptor::Barrier();

    //No need to initialize from raw data for the processor that
    //has just generated the table
    if(!ParallelDescriptor::IOProcessor()){
        m_shr_p_bw_engine->init_lookup_tables_from_raw_data(
            table_data, bw_minimum_chi_part);
    }
}

void
MultiParticleContainer::doQEDSchwinger ()
{
    WARPX_PROFILE("MultiParticleContainer::doQEDSchwinger()");

    if (!m_do_qed_schwinger) {return;}

    auto & warpx = WarpX::GetInstance();

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(warpx.do_nodal ||
       warpx.field_gathering_algo == GatheringAlgo::MomentumConserving,
          "ERROR: Schwinger process only implemented for warpx.do_nodal = 1"
                                 "or algo.field_gathering = momentum-conserving");

    constexpr int level_0 = 0;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(warpx.maxLevel() == level_0,
        "ERROR: Schwinger process not implemented with mesh refinement");

#ifdef WARPX_DIM_RZ
    amrex::Abort("Schwinger process not implemented in rz geometry");
#endif
#ifdef WARPX_DIM_1D_Z
    amrex::Abort("Schwinger process not implemented in 1D geometry");
#endif

// Get cell volume. In 2D the transverse size is
// chosen by the user in the input file.
    amrex::Geometry const & geom = warpx.Geom(level_0);
#if defined(WARPX_DIM_1D_Z)
    const auto dV = geom.CellSize(0); // TODO: scale properly
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const auto dV = geom.CellSize(0) * geom.CellSize(1)
        * m_qed_schwinger_y_size;
#elif defined(WARPX_DIM_3D)
    const auto dV = geom.CellSize(0) * geom.CellSize(1)
        * geom.CellSize(2);
#endif

   // Get the temporal step
   const auto dt =  warpx.getdt(level_0);

    auto& pc_product_ele =
            allcontainers[m_qed_schwinger_ele_product];
    auto& pc_product_pos =
            allcontainers[m_qed_schwinger_pos_product];

    pc_product_ele->defineAllParticleTiles();
    pc_product_pos->defineAllParticleTiles();

    const MultiFab & Ex = warpx.getEfield(level_0,0);
    const MultiFab & Ey = warpx.getEfield(level_0,1);
    const MultiFab & Ez = warpx.getEfield(level_0,2);
    const MultiFab & Bx = warpx.getBfield(level_0,0);
    const MultiFab & By = warpx.getBfield(level_0,1);
    const MultiFab & Bz = warpx.getBfield(level_0,2);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
     {
        // Make the box cell centered to avoid creating particles twice on the tile edges
        amrex::Box box = enclosedCells(mfi.nodaltilebox());

        // Get the box representing global Schwinger boundaries
        const amrex::Box global_schwinger_box = ComputeSchwingerGlobalBox();

        // If Schwinger process is not activated anywhere in the current box, we move to the next
        // one. Otherwise we use the intersection of current box with global Schwinger box.
        if (!box.intersects(global_schwinger_box)) {continue;}
        box &= global_schwinger_box;

        const MyFieldList fieldsEB = {
            Ex[mfi].array(), Ey[mfi].array(), Ez[mfi].array(),
            Bx[mfi].array(), By[mfi].array(), Bz[mfi].array()};

        auto& dst_ele_tile = pc_product_ele->ParticlesAt(level_0, mfi);
        auto& dst_pos_tile = pc_product_pos->ParticlesAt(level_0, mfi);

        const auto np_ele_dst = dst_ele_tile.numParticles();
        const auto np_pos_dst = dst_pos_tile.numParticles();

        const auto Filter  = SchwingerFilterFunc{
                              m_qed_schwinger_threshold_poisson_gaussian,
                              dV, dt};

        const SmartCreateFactory create_factory_ele(*pc_product_ele);
        const SmartCreateFactory create_factory_pos(*pc_product_pos);
        const auto CreateEle = create_factory_ele.getSmartCreate();
        const auto CreatePos = create_factory_pos.getSmartCreate();

        const auto Transform = SchwingerTransformFunc{m_qed_schwinger_y_size,
                            ParticleStringNames::to_index.find("w")->second};

        const auto num_added = filterCreateTransformFromFAB<1>( dst_ele_tile,
                               dst_pos_tile, box, fieldsEB, np_ele_dst,
                               np_pos_dst,Filter, CreateEle, CreatePos,
                               Transform);

        setNewParticleIDs(dst_ele_tile, np_ele_dst, num_added);
        setNewParticleIDs(dst_pos_tile, np_pos_dst, num_added);

    }
}

amrex::Box
MultiParticleContainer::ComputeSchwingerGlobalBox () const
{
    auto & warpx = WarpX::GetInstance();
    constexpr int level_0 = 0;
    amrex::Geometry const & geom = warpx.Geom(level_0);

#if defined(WARPX_DIM_3D)
    const amrex::Array<amrex::Real,3> schwinger_min{m_qed_schwinger_xmin,
                                                    m_qed_schwinger_ymin,
                                                    m_qed_schwinger_zmin};
    const amrex::Array<amrex::Real,3> schwinger_max{m_qed_schwinger_xmax,
                                                    m_qed_schwinger_ymax,
                                                    m_qed_schwinger_zmax};
#else
    const amrex::Array<amrex::Real,2> schwinger_min{m_qed_schwinger_xmin,
                                                    m_qed_schwinger_zmin};
    const amrex::Array<amrex::Real,2> schwinger_max{m_qed_schwinger_xmax,
                                                    m_qed_schwinger_zmax};
#endif

    // Box inside which Schwinger is activated
    amrex::Box schwinger_global_box;

    for (int dir=0; dir<AMREX_SPACEDIM; dir++)
    {
        // Dealing with these corner cases should ensure that we don't overflow on the integers
        if (schwinger_min[dir] < geom.ProbLo(dir))
        {
            schwinger_global_box.setSmall(dir, std::numeric_limits<int>::lowest());
        }
        else if (schwinger_min[dir] > geom.ProbHi(dir))
        {
            schwinger_global_box.setSmall(dir, std::numeric_limits<int>::max());
        }
        else
        {
            // Schwinger pairs are currently created on the lower nodes of a cell. Using ceil here
            // excludes all cells whose lower node is strictly lower than schwinger_min[dir].
            schwinger_global_box.setSmall(dir, static_cast<int>(std::ceil(
                               (schwinger_min[dir] - geom.ProbLo(dir)) / geom.CellSize(dir))));
        }

        if (schwinger_max[dir] < geom.ProbLo(dir))
        {
            schwinger_global_box.setBig(dir, std::numeric_limits<int>::lowest());
        }
        else if (schwinger_max[dir] > geom.ProbHi(dir))
        {
            schwinger_global_box.setBig(dir, std::numeric_limits<int>::max());
        }
        else
        {
            // Schwinger pairs are currently created on the lower nodes of a cell. Using floor here
            // excludes all cells whose lower node is strictly higher than schwinger_max[dir].
            schwinger_global_box.setBig(dir, static_cast<int>(std::floor(
                               (schwinger_max[dir] - geom.ProbLo(dir)) / geom.CellSize(dir))));
        }
    }

    return schwinger_global_box;
}

void MultiParticleContainer::doQedEvents (int lev,
                                          const MultiFab& Ex,
                                          const MultiFab& Ey,
                                          const MultiFab& Ez,
                                          const MultiFab& Bx,
                                          const MultiFab& By,
                                          const MultiFab& Bz)
{
    WARPX_PROFILE("MultiParticleContainer::doQedEvents()");

    doQedBreitWheeler(lev, Ex, Ey, Ez, Bx, By, Bz);
    doQedQuantumSync(lev, Ex, Ey, Ez, Bx, By, Bz);
}

void MultiParticleContainer::doQedBreitWheeler (int lev,
                                                const MultiFab& Ex,
                                                const MultiFab& Ey,
                                                const MultiFab& Ez,
                                                const MultiFab& Bx,
                                                const MultiFab& By,
                                                const MultiFab& Bz)
{
    WARPX_PROFILE("MultiParticleContainer::doQedBreitWheeler()");

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop over all species.
    // Photons undergoing Breit Wheeler process create electrons
    // in pc_product_ele and positrons in pc_product_pos

    for (auto& pc_source : allcontainers){
        if(!pc_source->has_breit_wheeler()) continue;

        // Get product species
        auto& pc_product_ele =
            allcontainers[pc_source->m_qed_breit_wheeler_ele_product];
        auto& pc_product_pos =
            allcontainers[pc_source->m_qed_breit_wheeler_pos_product];

        SmartCopyFactory copy_factory_ele(*pc_source, *pc_product_ele);
        SmartCopyFactory copy_factory_pos(*pc_source, *pc_product_pos);
        auto phys_pc_ptr = static_cast<PhysicalParticleContainer*>(pc_source.get());

        const auto Filter  = phys_pc_ptr->getPairGenerationFilterFunc();
        const auto CopyEle = copy_factory_ele.getSmartCopy();
        const auto CopyPos = copy_factory_pos.getSmartCopy();

        const auto pair_gen_functor = m_shr_p_bw_engine->build_pair_functor();

        pc_source ->defineAllParticleTiles();
        pc_product_pos->defineAllParticleTiles();
        pc_product_ele->defineAllParticleTiles();

        auto info = getMFItInfo(*pc_source, *pc_product_ele, *pc_product_pos);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*pc_source, lev, info); pti.isValid(); ++pti)
        {
            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            Real wt = amrex::second();

            auto Transform = PairGenerationTransformFunc(pair_gen_functor,
                                                         pti, lev, Ex.nGrowVect(),
                                                         Ex[pti], Ey[pti], Ez[pti],
                                                         Bx[pti], By[pti], Bz[pti]);

            auto& src_tile = pc_source->ParticlesAt(lev, pti);
            auto& dst_ele_tile = pc_product_ele->ParticlesAt(lev, pti);
            auto& dst_pos_tile = pc_product_pos->ParticlesAt(lev, pti);

            const auto np_dst_ele = dst_ele_tile.numParticles();
            const auto np_dst_pos = dst_pos_tile.numParticles();
            const auto num_added = filterCopyTransformParticles<1>(
                                                      dst_ele_tile, dst_pos_tile,
                                                      src_tile, np_dst_ele, np_dst_pos,
                                                      Filter, CopyEle, CopyPos, Transform);

            setNewParticleIDs(dst_ele_tile, np_dst_ele, num_added);
            setNewParticleIDs(dst_pos_tile, np_dst_pos, num_added);

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add( &(*cost)[pti.index()], wt);
            }
        }
    }
}

void MultiParticleContainer::doQedQuantumSync (int lev,
                                               const MultiFab& Ex,
                                               const MultiFab& Ey,
                                               const MultiFab& Ez,
                                               const MultiFab& Bx,
                                               const MultiFab& By,
                                               const MultiFab& Bz)
{
    WARPX_PROFILE("MultiParticleContainer::doQedQuantumSync()");

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

    // Loop over all species.
    // Electrons or positrons undergoing Quantum photon emission process
    // create photons in pc_product_phot

    for (auto& pc_source : allcontainers){
        if(!pc_source->has_quantum_sync()){ continue; }

        // Get product species
        auto& pc_product_phot =
            allcontainers[pc_source->m_qed_quantum_sync_phot_product];

        SmartCopyFactory copy_factory_phot(*pc_source, *pc_product_phot);
        auto phys_pc_ptr =
            static_cast<PhysicalParticleContainer*>(pc_source.get());

        const auto Filter   = phys_pc_ptr->getPhotonEmissionFilterFunc();
        const auto CopyPhot = copy_factory_phot.getSmartCopy();

        pc_source ->defineAllParticleTiles();
        pc_product_phot->defineAllParticleTiles();

        auto info = getMFItInfo(*pc_source, *pc_product_phot);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*pc_source, lev, info); pti.isValid(); ++pti)
        {
            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            Real wt = amrex::second();

            auto Transform = PhotonEmissionTransformFunc(
                  m_shr_p_qs_engine->build_optical_depth_functor(),
                  pc_source->particle_runtime_comps["opticalDepthQSR"],
                  m_shr_p_qs_engine->build_phot_em_functor(),
                  pti, lev, Ex.nGrowVect(),
                  Ex[pti], Ey[pti], Ez[pti],
                  Bx[pti], By[pti], Bz[pti]);

            auto& src_tile = pc_source->ParticlesAt(lev, pti);
            auto& dst_tile = pc_product_phot->ParticlesAt(lev, pti);

            const auto np_dst = dst_tile.numParticles();

            const auto num_added =
                filterCopyTransformParticles<1>(dst_tile, src_tile, np_dst,
                                                Filter, CopyPhot, Transform);

            setNewParticleIDs(dst_tile, np_dst, num_added);

            cleanLowEnergyPhotons(
                                  dst_tile, np_dst, num_added,
                                  m_quantum_sync_photon_creation_energy_threshold);

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add( &(*cost)[pti.index()], wt);
            }
        }
    }
}

void MultiParticleContainer::CheckQEDProductSpecies()
{
    auto const nspecies = static_cast<int>(species_names.size());
    for (int i=0; i<nspecies; i++){
        const auto& pc = allcontainers[i];
        if (pc->has_breit_wheeler()){
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != pc->m_qed_breit_wheeler_ele_product,
                "ERROR: Breit Wheeler product cannot be the same species");

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != pc->m_qed_breit_wheeler_pos_product,
                "ERROR: Breit Wheeler product cannot be the same species");

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                allcontainers[pc->m_qed_breit_wheeler_ele_product]->
                    AmIA<PhysicalSpecies::electron>()
                &&
                allcontainers[pc->m_qed_breit_wheeler_pos_product]->
                    AmIA<PhysicalSpecies::positron>(),
                "ERROR: Breit Wheeler product species are of wrong type");
        }

        if(pc->has_quantum_sync()){
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != pc->m_qed_quantum_sync_phot_product,
                "ERROR: Quantum Synchrotron product cannot be the same species");

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                allcontainers[pc->m_qed_quantum_sync_phot_product]->
                    AmIA<PhysicalSpecies::photon>(),
                "ERROR: Quantum Synchrotron product species is of wrong type");
        }
    }

    if (m_do_qed_schwinger) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                allcontainers[m_qed_schwinger_ele_product]->
                    AmIA<PhysicalSpecies::electron>()
                &&
                allcontainers[m_qed_schwinger_pos_product]->
                    AmIA<PhysicalSpecies::positron>(),
                "ERROR: Schwinger process product species are of wrong type");
    }

}

#endif
