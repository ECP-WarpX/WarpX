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
#include "SpeciesPhysicalProperties.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"
#ifdef WARPX_QED
    #include "Particles/ElementaryProcess/QEDInternals/SchwingerProcessWrapper.H"
    #include "Particles/ElementaryProcess/QEDSchwingerProcess.H"
    #include "Particles/ParticleCreation/FilterCreateTransformFromFAB.H"
#endif

#include <AMReX_Vector.H>

#include <limits>
#include <algorithm>
#include <string>
#include <vector>

using namespace amrex;

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

    // collision
    auto const ncollisions = collision_names.size();
    allcollisions.resize(ncollisions);
    for (int i = 0; i < static_cast<int>(ncollisions); ++i) {
        allcollisions[i] =
            std::make_unique<CollisionType>(species_names, collision_names[i]);
    }

}

void
MultiParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        ParmParse pp("particles");

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
        pp.query("B_ext_particle_init_style", m_B_ext_particle_s);
        std::transform(m_B_ext_particle_s.begin(),
                       m_B_ext_particle_s.end(),
                       m_B_ext_particle_s.begin(),
                       ::tolower);
        pp.query("E_ext_particle_init_style", m_E_ext_particle_s);
        std::transform(m_E_ext_particle_s.begin(),
                       m_E_ext_particle_s.end(),
                       m_E_ext_particle_s.begin(),
                       ::tolower);
        // if the input string for B_external on particles is "constant"
        // then the values for the external B on particles must
        // be provided in the input file.
        if (m_B_ext_particle_s == "constant")
            pp.getarr("B_external_particle", m_B_external_particle);

        // if the input string for E_external on particles is "constant"
        // then the values for the external E on particles must
        // be provided in the input file.
        if (m_E_ext_particle_s == "constant")
            pp.getarr("E_external_particle", m_E_external_particle);

        // if the input string for B_ext_particle_s is
        // "parse_b_ext_particle_function" then the mathematical expression
        // for the Bx_, By_, Bz_external_particle_function(x,y,z)
        // must be provided in the input file.
        if (m_B_ext_particle_s == "parse_b_ext_particle_function") {
           // store the mathematical expression as string
           std::string str_Bx_ext_particle_function;
           std::string str_By_ext_particle_function;
           std::string str_Bz_ext_particle_function;
           Store_parserString(pp, "Bx_external_particle_function(x,y,z,t)",
                                      str_Bx_ext_particle_function);
           Store_parserString(pp, "By_external_particle_function(x,y,z,t)",
                                      str_By_ext_particle_function);
           Store_parserString(pp, "Bz_external_particle_function(x,y,z,t)",
                                      str_Bz_ext_particle_function);

           // Parser for B_external on the particle
           m_Bx_particle_parser = std::make_unique<ParserWrapper<4>>(
                                    makeParser(str_Bx_ext_particle_function,{"x","y","z","t"}));
           m_By_particle_parser = std::make_unique<ParserWrapper<4>>(
                                    makeParser(str_By_ext_particle_function,{"x","y","z","t"}));
           m_Bz_particle_parser = std::make_unique<ParserWrapper<4>>(
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
           Store_parserString(pp, "Ex_external_particle_function(x,y,z,t)",
                                      str_Ex_ext_particle_function);
           Store_parserString(pp, "Ey_external_particle_function(x,y,z,t)",
                                      str_Ey_ext_particle_function);
           Store_parserString(pp, "Ez_external_particle_function(x,y,z,t)",
                                      str_Ez_ext_particle_function);
           // Parser for E_external on the particle
           m_Ex_particle_parser = std::make_unique<ParserWrapper<4>>(
                                    makeParser(str_Ex_ext_particle_function,{"x","y","z","t"}));
           m_Ey_particle_parser = std::make_unique<ParserWrapper<4>>(
                                    makeParser(str_Ey_ext_particle_function,{"x","y","z","t"}));
           m_Ez_particle_parser = std::make_unique<ParserWrapper<4>>(
                                    makeParser(str_Ez_ext_particle_function,{"x","y","z","t"}));

        }


        // particle species
        pp.queryarr("species_names", species_names);
        auto const nspecies = species_names.size();

        if (nspecies > 0) {
            // Get species to deposit on main grid
            m_deposit_on_main_grid.resize(nspecies, false);
            std::vector<std::string> tmp;
            pp.queryarr("deposit_on_main_grid", tmp);
            for (auto const& name : tmp) {
                auto it = std::find(species_names.begin(), species_names.end(), name);
                WarpXUtilMsg::AlwaysAssert(
                    it != species_names.end(),
                    "ERROR: species '" + name
                    + "' in particles.deposit_on_main_grid must be part of particles.species_names"
                );
                int i = std::distance(species_names.begin(), it);
                m_deposit_on_main_grid[i] = true;
            }

            m_gather_from_main_grid.resize(nspecies, false);
            std::vector<std::string> tmp_gather;
            pp.queryarr("gather_from_main_grid", tmp_gather);
            for (auto const& name : tmp_gather) {
                auto it = std::find(species_names.begin(), species_names.end(), name);
                WarpXUtilMsg::AlwaysAssert(
                    it != species_names.end(),
                    "ERROR: species '" + name
                    + "' in particles.gather_from_main_grid must be part of particles.species_names"
                );
                int i = std::distance(species_names.begin(), it);
                m_gather_from_main_grid.at(i) = true;
            }

            species_types.resize(nspecies, PCTypes::Physical);

            // Get rigid-injected species
            std::vector<std::string> rigid_injected_species;
            pp.queryarr("rigid_injected_species", rigid_injected_species);
            if (!rigid_injected_species.empty()) {
                for (auto const& name : rigid_injected_species) {
                    auto it = std::find(species_names.begin(), species_names.end(), name);
                    WarpXUtilMsg::AlwaysAssert(
                        it != species_names.end(),
                        "ERROR: species '" + name
                        + "' in particles.rigid_injected_species must be part of particles.species_names"
                    );
                    int i = std::distance(species_names.begin(), it);
                    species_types[i] = PCTypes::RigidInjected;
                }
            }
            // Get photon species
            std::vector<std::string> photon_species;
            pp.queryarr("photon_species", photon_species);
            if (!photon_species.empty()) {
                for (auto const& name : photon_species) {
                    auto it = std::find(species_names.begin(), species_names.end(), name);
                    WarpXUtilMsg::AlwaysAssert(
                        it != species_names.end(),
                        "ERROR: species '" + name
                        + "' in particles.rigid_injected_species must be part of particles.species_names"
                    );
                    int i = std::distance(species_names.begin(), it);
                    species_types[i] = PCTypes::Photon;
                }
            }

            // binary collisions
            ParmParse pc("collisions");
            pc.queryarr("collision_names", collision_names);

        }
        pp.query("use_fdtd_nci_corr", WarpX::use_fdtd_nci_corr);
#ifdef WARPX_DIM_RZ
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::use_fdtd_nci_corr==0,
                            "ERROR: use_fdtd_nci_corr is not supported in RZ");
#endif
        pp.query("galerkin_interpolation", WarpX::galerkin_interpolation);

        std::string boundary_conditions = "none";
        pp.query("boundary_conditions", boundary_conditions);
        if        (boundary_conditions == "none"){
            m_boundary_conditions = ParticleBC::none;
        } else if (boundary_conditions == "absorbing"){
            m_boundary_conditions = ParticleBC::absorbing;
        } else {
            amrex::Abort("unknown particle BC type");
        }

        ParmParse ppl("lasers");
        ppl.queryarr("names", lasers_names);

#ifdef WARPX_QED
        ParmParse ppw("warpx");
        ppw.query("do_qed_schwinger", m_do_qed_schwinger);

        if (m_do_qed_schwinger) {
            ParmParse ppq("qed_schwinger");
            ppq.get("ele_product_species", m_qed_schwinger_ele_product_name);
            ppq.get("pos_product_species", m_qed_schwinger_pos_product_name);
#if (AMREX_SPACEDIM == 2)
            ppq.get("y_size",m_qed_schwinger_y_size);
#endif
            ppq.query("threshold_poisson_gaussian",
                      m_qed_schwinger_threshold_poisson_gaussian);
        }
#endif
        initialized = true;
    }
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
                                const MultiFab& Ex_avg, const MultiFab& Ey_avg, const MultiFab& Ez_avg,
                                const MultiFab& Bx_avg, const MultiFab& By_avg, const MultiFab& Bz_avg,
                                MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                MultiFab* cjx,  MultiFab* cjy, MultiFab* cjz,
                                MultiFab* rho, MultiFab* crho,
                                const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                Real t, Real dt, DtType a_dt_type)
{
    jx.setVal(0.0);
    jy.setVal(0.0);
    jz.setVal(0.0);
    if (cjx) cjx->setVal(0.0);
    if (cjy) cjy->setVal(0.0);
    if (cjz) cjz->setVal(0.0);
    if (rho) rho->setVal(0.0);
    if (crho) crho->setVal(0.0);
    for (auto& pc : allcontainers) {
        pc->Evolve(lev, Ex, Ey, Ez, Bx, By, Bz, Ex_avg, Ey_avg, Ez_avg, Bx_avg, By_avg, Bz_avg, jx, jy, jz, cjx, cjy, cjz,
                   rho, crho, cEx, cEy, cEz, cBx, cBy, cBz, t, dt, a_dt_type);
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

std::unique_ptr<MultiFab>
MultiParticleContainer::GetChargeDensity (int lev, bool local)
{
    if (allcontainers.size() == 0)
    {
        std::unique_ptr<MultiFab> rho = GetZeroChargeDensity(lev);
        return rho;
    }
    else
    {
        std::unique_ptr<MultiFab> rho = allcontainers[0]->GetChargeDensity(lev, true);
        for (unsigned i = 1, n = allcontainers.size(); i < n; ++i) {
            std::unique_ptr<MultiFab> rhoi = allcontainers[i]->GetChargeDensity(lev, true);
            MultiFab::Add(*rho, *rhoi, 0, 0, rho->nComp(), rho->nGrow());
        }
        if (!local) {
            const Geometry& gm = allcontainers[0]->Geom(lev);
            rho->SumBoundary(gm.periodicity());
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
        pc->ApplyBoundaryConditions(m_boundary_conditions);
    }
}

Vector<long>
MultiParticleContainer::GetZeroParticlesInGrid (const int lev) const
{
    WarpX& warpx = WarpX::GetInstance();
    const int num_boxes = warpx.boxArray(lev).size();
    const Vector<Long> r(num_boxes, 0);
    return r;
}

Vector<long>
MultiParticleContainer::NumberOfParticlesInGrid (int lev) const
{
    if (allcontainers.size() == 0)
    {
        const Vector<long> r = GetZeroParticlesInGrid(lev);
        return r;
    }
    else
    {
        const bool only_valid=true, only_local=true;
        Vector<long> r = allcontainers[0]->NumberOfParticlesInGrid(lev,only_valid,only_local);
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
            for (auto it = diagnostic_particles[lev].begin(); it != diagnostic_particles[lev].end(); ++it){
                // it->first is the [grid index][tile index] key
                // it->second is the corresponding
                // WarpXParticleContainer::DiagnosticParticleData value
                parts[i].GetRealData(DiagIdx::w).insert(  parts[i].GetRealData(DiagIdx::w  ).end(),
                                                          it->second.GetRealData(DiagIdx::w  ).begin(),
                                                          it->second.GetRealData(DiagIdx::w  ).end());

                parts[i].GetRealData(DiagIdx::x).insert(  parts[i].GetRealData(DiagIdx::x  ).end(),
                                                          it->second.GetRealData(DiagIdx::x  ).begin(),
                                                          it->second.GetRealData(DiagIdx::x  ).end());

                parts[i].GetRealData(DiagIdx::y).insert(  parts[i].GetRealData(DiagIdx::y  ).end(),
                                                          it->second.GetRealData(DiagIdx::y  ).begin(),
                                                          it->second.GetRealData(DiagIdx::y  ).end());

                parts[i].GetRealData(DiagIdx::z).insert(  parts[i].GetRealData(DiagIdx::z  ).end(),
                                                          it->second.GetRealData(DiagIdx::z  ).begin(),
                                                          it->second.GetRealData(DiagIdx::z  ).end());

                parts[i].GetRealData(DiagIdx::ux).insert(  parts[i].GetRealData(DiagIdx::ux).end(),
                                                           it->second.GetRealData(DiagIdx::ux).begin(),
                                                           it->second.GetRealData(DiagIdx::ux).end());

                parts[i].GetRealData(DiagIdx::uy).insert(  parts[i].GetRealData(DiagIdx::uy).end(),
                                                           it->second.GetRealData(DiagIdx::uy).begin(),
                                                           it->second.GetRealData(DiagIdx::uy).end());

                parts[i].GetRealData(DiagIdx::uz).insert(  parts[i].GetRealData(DiagIdx::uz).end(),
                                                           it->second.GetRealData(DiagIdx::uz).begin(),
                                                           it->second.GetRealData(DiagIdx::uz).end());
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
    int i_product = 0;
    bool found = 0;
    // Loop over species
    for (int i=0; i < static_cast<int>(species_names.size()); i++){
        // If species name matches, store its ID
        // into i_product
        if (species_names[i] == product_str){
            found = 1;
            i_product = i;
        }
    }

    WarpXUtilMsg::AlwaysAssert(
        found != 0,
        "ERROR: could not find the ID of product species '"
        + product_str + "'" + ". Wrong name?"
    );

    return i_product;
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*pc_source, lev, info); pti.isValid(); ++pti)
        {
            auto& src_tile = pc_source ->ParticlesAt(lev, pti);
            auto& dst_tile = pc_product->ParticlesAt(lev, pti);

            auto Filter = phys_pc_ptr->getIonizationFunc(pti, lev, Ex.nGrow(),
                                                         Ex[pti], Ey[pti], Ez[pti],
                                                         Bx[pti], By[pti], Bz[pti]);

            const auto np_dst = dst_tile.numParticles();
            const auto num_added = filterCopyTransformParticles<1>(dst_tile, src_tile, np_dst,
                                                                   Filter, Copy, Transform);

            setNewParticleIDs(dst_tile, np_dst, num_added);
        }
    }
}

void
MultiParticleContainer::doCoulombCollisions ( Real cur_time )
{
    WARPX_PROFILE("MultiParticleContainer::doCoulombCollisions()");

    for( auto const& collision : allcollisions )
    {

        const Real dt = WarpX::GetInstance().getdt(0);
        if ( int(std::floor(cur_time/dt)) % collision->m_ndt != 0 ) continue;

        auto& species1 = allcontainers[ collision->m_species1_index ];
        auto& species2 = allcontainers[ collision->m_species2_index ];

        // Enable tiling
        MFItInfo info;
        if (Gpu::notInLaunchRegion()) info.EnableTiling(species1->tile_size);

        // Loop over refinement levels
        for (int lev = 0; lev <= species1->finestLevel(); ++lev){

            // Loop over all grids/tiles at this level
#ifdef _OPENMP
            info.SetDynamic(true);
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi = species1->MakeMFIter(lev, info); mfi.isValid(); ++mfi){

                CollisionType::doCoulombCollisionsWithinTile
                    ( lev, mfi, species1, species2,
                      collision->m_isSameSpecies,
                      collision->m_CoulombLog,
                      collision->m_ndt );

            }
        }
    }
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
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != allcontainers[i]->ionization_product,
                "ERROR: ionization product cannot be the same species");
        }
    }
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
    ParmParse pp("qed_qs");

    //If specified, use a user-defined energy threshold for photon creaction
    ParticleReal temp;
    constexpr auto mec2 = PhysConst::c * PhysConst::c * PhysConst::m_e;
    if(pp.query("photon_creation_energy_threshold", temp)){
        temp *= mec2;
        m_quantum_sync_photon_creation_energy_threshold = temp;
    }
    else{
        amrex::Print() << "Using default value (2*me*c^2)" <<
            " for photon energy creaction threshold \n" ;
    }

    // qs_minimum_chi_part is the minimum chi parameter to be
    // considered for Synchrotron emission. If a lepton has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real qs_minimum_chi_part;
    pp.get("chi_min", qs_minimum_chi_part);


    pp.query("lookup_table_mode", lookup_table_mode);
    if(lookup_table_mode.empty()){
        amrex::Abort("Quantum Synchrotron table mode should be provided");
    }

    if(lookup_table_mode == "generate"){
        amrex::Print() << "Quantum Synchrotron table will be generated. \n" ;
#ifndef WARPX_QED_TABLE_GEN
        amrex::Error("Error: Compile with QED_TABLE_GEN=TRUE to enable table generation!\n");
#else
        QuantumSyncGenerateTable();
#endif
    }
    else if(lookup_table_mode == "load"){
        amrex::Print() << "Quantum Synchrotron table will be read from file. \n" ;
        std::string load_table_name;
        pp.query("load_table_from", load_table_name);
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
        amrex::Print() << "Built-in Quantum Synchrotron table will be used. \n" ;
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
    ParmParse pp("qed_bw");

    // bw_minimum_chi_phot is the minimum chi parameter to be
    // considered for pair production. If a photon has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real bw_minimum_chi_part;
    if(!pp.query("chi_min", bw_minimum_chi_part))
        amrex::Abort("qed_bw.chi_min should be provided!");

    pp.query("lookup_table_mode", lookup_table_mode);
    if(lookup_table_mode.empty()){
        amrex::Abort("Breit Wheeler table mode should be provided");
    }

    if(lookup_table_mode == "generate"){
        amrex::Print() << "Breit Wheeler table will be generated. \n" ;
#ifndef WARPX_QED_TABLE_GEN
        amrex::Error("Error: Compile with QED_TABLE_GEN=TRUE to enable table generation!\n");
#else
        BreitWheelerGenerateTable();
#endif
    }
    else if(lookup_table_mode == "load"){
        amrex::Print() << "Breit Wheeler table will be read from file. \n" ;
        std::string load_table_name;
        pp.query("load_table_from", load_table_name);
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
        amrex::Print() << "Built-in Breit Wheeler table will be used. \n" ;
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
    ParmParse pp("qed_qs");
    std::string table_name;
    pp.query("save_table_in", table_name);
    if(table_name.empty())
        amrex::Abort("qed_qs.save_table_in should be provided!");

    // qs_minimum_chi_part is the minimum chi parameter to be
    // considered for Synchrotron emission. If a lepton has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real qs_minimum_chi_part;
    pp.get("chi_min", qs_minimum_chi_part);

    if(ParallelDescriptor::IOProcessor()){
        PicsarQuantumSyncCtrl ctrl;

        //==Table parameters==

        //--- sub-table 1 (1D)
        //These parameters are used to pre-compute a function
        //which appears in the evolution of the optical depth

        //Minimun chi for the table. If a lepton has chi < tab_dndt_chi_min,
        //chi is considered as if it were equal to tab_dndt_chi_min
        pp.get("tab_dndt_chi_min", ctrl.dndt_params.chi_part_min);

        //Maximum chi for the table. If a lepton has chi > tab_dndt_chi_max,
        //chi is considered as if it were equal to tab_dndt_chi_max
        pp.get("tab_dndt_chi_max", ctrl.dndt_params.chi_part_max);

        //How many points should be used for chi in the table
        pp.get("tab_dndt_how_many", ctrl.dndt_params.chi_part_how_many);
        //------

        //--- sub-table 2 (2D)
        //These parameters are used to pre-compute a function
        //which is used to extract the properties of the generated
        //photons.

        //Minimun chi for the table. If a lepton has chi < tab_em_chi_min,
        //chi is considered as if it were equal to tab_em_chi_min
        pp.get("tab_em_chi_min", ctrl.phot_em_params.chi_part_min);

        //Maximum chi for the table. If a lepton has chi > tab_em_chi_max,
        //chi is considered as if it were equal to tab_em_chi_max
        pp.get("tab_em_chi_max", ctrl.phot_em_params.chi_part_max);

        //How many points should be used for chi in the table
        pp.get("tab_em_chi_how_many", ctrl.phot_em_params.chi_part_how_many);

        //The other axis of the table is the ratio between the quantum
        //parameter of the emitted photon and the quantum parameter of the
        //lepton. This parameter is the minimum ratio to consider for the table.
        pp.get("tab_em_frac_min", ctrl.phot_em_params.frac_min);

        //This parameter is the number of different points to consider for the second
        //axis
        pp.get("tab_em_frac_how_many", ctrl.phot_em_params.frac_how_many);
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
    ParmParse pp("qed_bw");
    std::string table_name;
    pp.query("save_table_in", table_name);
    if(table_name.empty())
        amrex::Abort("qed_bw.save_table_in should be provided!");

    // bw_minimum_chi_phot is the minimum chi parameter to be
    // considered for pair production. If a photon has chi < chi_min,
    // the optical depth is not evolved and photon generation is ignored
    amrex::Real bw_minimum_chi_part;
    pp.get("chi_min", bw_minimum_chi_part);

    if(ParallelDescriptor::IOProcessor()){
        PicsarBreitWheelerCtrl ctrl;

        //==Table parameters==

        //--- sub-table 1 (1D)
        //These parameters are used to pre-compute a function
        //which appears in the evolution of the optical depth

        //Minimun chi for the table. If a photon has chi < tab_dndt_chi_min,
        //an analytical approximation is used.
        pp.get("tab_dndt_chi_min", ctrl.dndt_params.chi_phot_min);

        //Maximum chi for the table. If a photon has chi > tab_dndt_chi_max,
        //an analytical approximation is used.
        pp.get("tab_dndt_chi_max", ctrl.dndt_params.chi_phot_max);

        //How many points should be used for chi in the table
        pp.get("tab_dndt_how_many", ctrl.dndt_params.chi_phot_how_many);
        //------

        //--- sub-table 2 (2D)
        //These parameters are used to pre-compute a function
        //which is used to extract the properties of the generated
        //particles.

        //Minimun chi for the table. If a photon has chi < tab_pair_chi_min
        //chi is considered as it were equal to chi_phot_tpair_min
        pp.get("tab_pair_chi_min", ctrl.pair_prod_params.chi_phot_min);

        //Maximum chi for the table. If a photon has chi > tab_pair_chi_max
        //chi is considered as it were equal to chi_phot_tpair_max
        pp.get("tab_pair_chi_max", ctrl.pair_prod_params.chi_phot_max);

        //How many points should be used for chi in the table
        pp.get("tab_pair_chi_how_many", ctrl.pair_prod_params.chi_phot_how_many);

        //The other axis of the table is the fraction of the initial energy
        //'taken away' by the most energetic particle of the pair.
        //This parameter is the number of different fractions to consider
        pp.get("tab_pair_frac_how_many", ctrl.pair_prod_params.frac_how_many);
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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(warpx.do_nodal ||
       warpx.field_gathering_algo == GatheringAlgo::MomentumConserving,
          "ERROR: Schwinger process only implemented for warpx.do_nodal = 1"
                                 "or algo.field_gathering = momentum-conserving");

    const int level_0 = 0;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(warpx.maxLevel() == level_0,
        "ERROR: Schwinger process not implemented with mesh refinement");

#ifdef WARPX_DIM_RZ
    amrex::Abort("Schwinger process not implemented in rz geometry");
#endif

// Get cell volume. In 2D the transverse size is
// chosen by the user in the input file.
    amrex::Geometry const & geom = warpx.Geom(level_0);
#if (AMREX_SPACEDIM == 2)
    const auto dV = geom.CellSize(0) * geom.CellSize(1)
        * m_qed_schwinger_y_size;
#elif (AMREX_SPACEDIM == 3)
    const auto dV = geom.CellSize(0) * geom.CellSize(1)
        * geom.CellSize(2);
#endif

   // Get the temporal step
   const auto dt =  warpx.getdt(level_0);

    auto& pc_product_ele =
            allcontainers[m_qed_schwinger_ele_product];
    auto& pc_product_pos =
            allcontainers[m_qed_schwinger_pos_product];

    const MultiFab & Ex = warpx.getEfield(level_0,0);
    const MultiFab & Ey = warpx.getEfield(level_0,1);
    const MultiFab & Ez = warpx.getEfield(level_0,2);
    const MultiFab & Bx = warpx.getBfield(level_0,0);
    const MultiFab & By = warpx.getBfield(level_0,1);
    const MultiFab & Bz = warpx.getBfield(level_0,2);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
     {
        // Make the box cell centered to avoid creating particles twice on the tile edges
        const Box& box = enclosedCells(mfi.nodaltilebox());

        const auto& arrEx = Ex[mfi].array();
        const auto& arrEy = Ey[mfi].array();
        const auto& arrEz = Ez[mfi].array();
        const auto& arrBx = Bx[mfi].array();
        const auto& arrBy = By[mfi].array();
        const auto& arrBz = Bz[mfi].array();

        const Array4<const amrex::Real> array_EMFAB [] = {arrEx,arrEy,arrEz,
                                           arrBx,arrBy,arrBz};

        pc_product_ele->defineAllParticleTiles();
        pc_product_pos->defineAllParticleTiles();

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
                              dst_pos_tile, box, array_EMFAB, np_ele_dst,
                               np_pos_dst,Filter, CreateEle, CreatePos,
                                Transform);

        setNewParticleIDs(dst_ele_tile, np_ele_dst, num_added);
        setNewParticleIDs(dst_pos_tile, np_pos_dst, num_added);

    }
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*pc_source, lev, info); pti.isValid(); ++pti)
        {
            auto Transform = PairGenerationTransformFunc(pair_gen_functor,
                                                         pti, lev, Ex.nGrow(),
                                                         Ex[pti], Ey[pti], Ez[pti],
                                                         Bx[pti], By[pti], Bz[pti],
                                                         pc_source->get_v_galilean());

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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*pc_source, lev, info); pti.isValid(); ++pti)
        {
            auto Transform = PhotonEmissionTransformFunc(
                  m_shr_p_qs_engine->build_optical_depth_functor(),
                  pc_source->particle_runtime_comps["optical_depth_QSR"],
                  m_shr_p_qs_engine->build_phot_em_functor(),
                  pti, lev, Ex.nGrow(),
                  Ex[pti], Ey[pti], Ez[pti],
                  Bx[pti], By[pti], Bz[pti],
                  pc_source->get_v_galilean());

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
        }
    }
}

void MultiParticleContainer::CheckQEDProductSpecies()
{
    auto const nspecies = static_cast<int>(species_names.size());
    for (int i=0; i<nspecies; i++){
        const auto& pc = allcontainers[i];
        if (pc->has_breit_wheeler()){
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != pc->m_qed_breit_wheeler_ele_product,
                "ERROR: Breit Wheeler product cannot be the same species");

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != pc->m_qed_breit_wheeler_pos_product,
                "ERROR: Breit Wheeler product cannot be the same species");

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                allcontainers[pc->m_qed_breit_wheeler_ele_product]->
                    AmIA<PhysicalSpecies::electron>()
                &&
                allcontainers[pc->m_qed_breit_wheeler_pos_product]->
                    AmIA<PhysicalSpecies::positron>(),
                "ERROR: Breit Wheeler product species are of wrong type");
        }

        if(pc->has_quantum_sync()){
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != pc->m_qed_quantum_sync_phot_product,
                "ERROR: Quantum Synchrotron product cannot be the same species");

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                allcontainers[pc->m_qed_quantum_sync_phot_product]->
                    AmIA<PhysicalSpecies::photon>(),
                "ERROR: Quantum Synchrotron product species is of wrong type");
        }
    }

    if (m_do_qed_schwinger) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                allcontainers[m_qed_schwinger_ele_product]->
                    AmIA<PhysicalSpecies::electron>()
                &&
                allcontainers[m_qed_schwinger_pos_product]->
                    AmIA<PhysicalSpecies::positron>(),
                "ERROR: Schwinger process product species are of wrong type");
    }

}

#endif
