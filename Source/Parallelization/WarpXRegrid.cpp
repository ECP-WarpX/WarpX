/* Copyright 2019 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Maxence Thevenet, Michael Rowan
 * Remi Lehe, Weiqun Zhang, levinem
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Diagnostics/MultiDiagnostics.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_IndexType.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MakeType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_iMultiFab.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

using namespace amrex;

void
WarpX::LoadBalance ()
{
    WARPX_PROFILE_REGION("LoadBalance");
    WARPX_PROFILE("WarpX::LoadBalance()");

    AMREX_ALWAYS_ASSERT(costs[0] != nullptr);

#ifdef AMREX_USE_MPI
    if (load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Heuristic)
    {
        // compute the costs on a per-rank basis
        ComputeCostsHeuristic(costs);
    }

    // By default, do not do a redistribute; this toggles to true if RemakeLevel
    // is called for any level
    int loadBalancedAnyLevel = false;

    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev)
    {
        int doLoadBalance = false;

        // Compute the new distribution mapping
        DistributionMapping newdm;
        const amrex::Real nboxes = costs[lev]->size();
        const amrex::Real nprocs = ParallelContext::NProcsSub();
        const int nmax = static_cast<int>(std::ceil(nboxes/nprocs*load_balance_knapsack_factor));
        // These store efficiency (meaning, the  average 'cost' over all ranks,
        // normalized to max cost) for current and proposed distribution mappings
        amrex::Real currentEfficiency = 0.0;
        amrex::Real proposedEfficiency = 0.0;

        newdm = (load_balance_with_sfc)
            ? DistributionMapping::makeSFC(*costs[lev],
                                           currentEfficiency, proposedEfficiency,
                                           false,
                                           ParallelDescriptor::IOProcessorNumber())
            : DistributionMapping::makeKnapSack(*costs[lev],
                                                currentEfficiency, proposedEfficiency,
                                                nmax,
                                                false,
                                                ParallelDescriptor::IOProcessorNumber());
        // As specified in the above calls to makeSFC and makeKnapSack, the new
        // distribution mapping is NOT communicated to all ranks; the loadbalanced
        // dm is up-to-date only on root, and we can decide whether to broadcast
        if ((load_balance_efficiency_ratio_threshold > 0.0)
            && (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()))
        {
            doLoadBalance = (proposedEfficiency > load_balance_efficiency_ratio_threshold*currentEfficiency);
        }

        ParallelDescriptor::Bcast(&doLoadBalance, 1,
                                  ParallelDescriptor::IOProcessorNumber());

        if (doLoadBalance)
        {
            Vector<int> pmap;
            if (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber())
            {
                pmap = newdm.ProcessorMap();
            } else
            {
                pmap.resize(static_cast<std::size_t>(nboxes));
            }
            ParallelDescriptor::Bcast(&pmap[0], pmap.size(), ParallelDescriptor::IOProcessorNumber());

            if (ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber())
            {
                newdm = DistributionMapping(pmap);
            }

            RemakeLevel(lev, t_new[lev], boxArray(lev), newdm);

            // Record the load balance efficiency
            setLoadBalanceEfficiency(lev, proposedEfficiency);
        }

        loadBalancedAnyLevel = loadBalancedAnyLevel || doLoadBalance;
    }
    if (loadBalancedAnyLevel)
    {
        mypc->Redistribute();
        mypc->defineAllParticleTiles();
    }
#endif
}


void
WarpX::RemakeLevel (int lev, Real /*time*/, const BoxArray& ba, const DistributionMapping& dm)
{
    if (ba == boxArray(lev))
    {
        if (ParallelDescriptor::NProcs() == 1) return;

#ifdef AMREX_USE_EB
        m_field_factory[lev] = amrex::makeEBFabFactory(Geom(lev), ba, dm,
                                                       {1,1,1}, // Not clear how many ghost cells we need yet
                                                       amrex::EBSupport::full);
#else
        m_field_factory[lev] = std::make_unique<FArrayBoxFactory>();
#endif

        // Fine patch
        for (int idim=0; idim < 3; ++idim)
        {
            {
                const IntVect& ng = Bfield_fp[lev][idim]->nGrowVect();
                auto pmf = std::make_unique<MultiFab>(Bfield_fp[lev][idim]->boxArray(),
                                                                  dm, Bfield_fp[lev][idim]->nComp(), ng);
                pmf->Redistribute(*Bfield_fp[lev][idim], 0, 0, Bfield_fp[lev][idim]->nComp(), ng);
                Bfield_fp[lev][idim] = std::move(pmf);
            }
            {
                const IntVect& ng = Efield_fp[lev][idim]->nGrowVect();
                auto pmf = std::make_unique<MultiFab>(Efield_fp[lev][idim]->boxArray(),
                                                      dm, Efield_fp[lev][idim]->nComp(), ng);
                pmf->Redistribute(*Efield_fp[lev][idim], 0, 0, Efield_fp[lev][idim]->nComp(), ng);
                Efield_fp[lev][idim] = std::move(pmf);
            }
            {
                const IntVect& ng = current_fp[lev][idim]->nGrowVect();
                auto pmf = std::make_unique<MultiFab>(current_fp[lev][idim]->boxArray(),
                                                      dm, current_fp[lev][idim]->nComp(), ng);
                current_fp[lev][idim] = std::move(pmf);
            }
            if (current_store[lev][idim])
            {
                const IntVect& ng = current_store[lev][idim]->nGrowVect();
                auto pmf = std::make_unique<MultiFab>(current_store[lev][idim]->boxArray(),
                                                      dm, current_store[lev][idim]->nComp(), ng);
                // no need to redistribute
                current_store[lev][idim] = std::move(pmf);
            }
        }

        if (F_fp[lev] != nullptr) {
            const IntVect& ng = F_fp[lev]->nGrowVect();
            auto pmf = std::make_unique<MultiFab>(F_fp[lev]->boxArray(),
                                                              dm, F_fp[lev]->nComp(), ng);
            pmf->Redistribute(*F_fp[lev], 0, 0, F_fp[lev]->nComp(), ng);
            F_fp[lev] = std::move(pmf);
        }

        if (rho_fp[lev] != nullptr) {
            const int nc = rho_fp[lev]->nComp();
            const IntVect& ng = rho_fp[lev]->nGrowVect();
            auto pmf = std::make_unique<MultiFab>(rho_fp[lev]->boxArray(),
                                                              dm, nc, ng);
            rho_fp[lev] = std::move(pmf);
        }

#ifdef WARPX_USE_PSATD
        if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            if (spectral_solver_fp[lev] != nullptr) {
                // Get the cell-centered box
                BoxArray realspace_ba = ba;   // Copy box
                realspace_ba.enclosedCells(); // Make it cell-centered
                auto ngE = getngE();
                auto dx = CellSize(lev);

#   ifdef WARPX_DIM_RZ
                if ( fft_periodic_single_box == false ) {
                    realspace_ba.grow(1, ngE[1]); // add guard cells only in z
                }
                AllocLevelSpectralSolverRZ(spectral_solver_fp,
                                           lev,
                                           realspace_ba,
                                           dm,
                                           dx);
#   else
                if ( fft_periodic_single_box == false ) {
                    realspace_ba.grow(ngE);   // add guard cells
                }
                bool const pml_flag_false = false;
                AllocLevelSpectralSolver(spectral_solver_fp,
                                         lev,
                                         realspace_ba,
                                         dm,
                                         dx,
                                         pml_flag_false);
#   endif
            }
        }
#endif

        // Aux patch
        if (lev == 0 && Bfield_aux[0][0]->ixType() == Bfield_fp[0][0]->ixType())
        {
            for (int idim = 0; idim < 3; ++idim) {
                Bfield_aux[lev][idim] = std::make_unique<MultiFab>(*Bfield_fp[lev][idim], amrex::make_alias, 0, Bfield_aux[lev][idim]->nComp());
                Efield_aux[lev][idim] = std::make_unique<MultiFab>(*Efield_fp[lev][idim], amrex::make_alias, 0, Efield_aux[lev][idim]->nComp());
            }
        } else {
            for (int idim=0; idim < 3; ++idim)
            {
                {
                    const IntVect& ng = Bfield_aux[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(Bfield_aux[lev][idim]->boxArray(),
                                                                      dm, Bfield_aux[lev][idim]->nComp(), ng);
                    // pmf->Redistribute(*Bfield_aux[lev][idim], 0, 0, Bfield_aux[lev][idim]->nComp(), ng);
                    Bfield_aux[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = Efield_aux[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(Efield_aux[lev][idim]->boxArray(),
                                                                      dm, Efield_aux[lev][idim]->nComp(), ng);
                    // pmf->Redistribute(*Efield_aux[lev][idim], 0, 0, Efield_aux[lev][idim]->nComp(), ng);
                    Efield_aux[lev][idim] = std::move(pmf);
                }
            }
        }

        // Coarse patch
        if (lev > 0) {
            for (int idim=0; idim < 3; ++idim)
            {
                {
                    const IntVect& ng = Bfield_cp[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(Bfield_cp[lev][idim]->boxArray(),
                                                                      dm, Bfield_cp[lev][idim]->nComp(), ng);
                    pmf->Redistribute(*Bfield_cp[lev][idim], 0, 0, Bfield_cp[lev][idim]->nComp(), ng);
                    Bfield_cp[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = Efield_cp[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(Efield_cp[lev][idim]->boxArray(),
                                                                      dm, Efield_cp[lev][idim]->nComp(), ng);
                    pmf->Redistribute(*Efield_cp[lev][idim], 0, 0, Efield_cp[lev][idim]->nComp(), ng);
                    Efield_cp[lev][idim] = std::move(pmf);
                }
                {
                    const IntVect& ng = current_cp[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(current_cp[lev][idim]->boxArray(),
                                                                       dm, current_cp[lev][idim]->nComp(), ng);
                    current_cp[lev][idim] = std::move(pmf);
                }
            }

            if (F_cp[lev] != nullptr) {
                const IntVect& ng = F_cp[lev]->nGrowVect();
                auto pmf = std::make_unique<MultiFab>(F_cp[lev]->boxArray(),
                                                                  dm, F_cp[lev]->nComp(), ng);
                pmf->Redistribute(*F_cp[lev], 0, 0, F_cp[lev]->nComp(), ng);
                F_cp[lev] = std::move(pmf);
            }

            if (rho_cp[lev] != nullptr) {
                const int nc = rho_cp[lev]->nComp();
                const IntVect& ng = rho_cp[lev]->nGrowVect();
                auto pmf  = std::make_unique<MultiFab>(rho_cp[lev]->boxArray(),
                                                                  dm, nc, ng);
                rho_cp[lev] = std::move(pmf);
            }

#ifdef WARPX_USE_PSATD
            if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
                if (spectral_solver_cp[lev] != nullptr) {
                    BoxArray cba = ba;
                    cba.coarsen(refRatio(lev-1));
                    std::array<Real,3> cdx = CellSize(lev-1);

                    // Get the cell-centered box
                    BoxArray c_realspace_ba = cba;  // Copy box
                    c_realspace_ba.enclosedCells(); // Make it cell-centered

                    auto ngE = getngE();

#   ifdef WARPX_DIM_RZ
                    c_realspace_ba.grow(1, ngE[1]); // add guard cells only in z
                    AllocLevelSpectralSolverRZ(spectral_solver_cp,
                                               lev,
                                               c_realspace_ba,
                                               dm,
                                               cdx);
#   else
                    c_realspace_ba.grow(ngE);
                    bool const pml_flag_false = false;
                    AllocLevelSpectralSolver(spectral_solver_cp,
                                             lev,
                                             c_realspace_ba,
                                             dm,
                                             cdx,
                                             pml_flag_false);
#   endif
                }
            }
#endif
        }

        if (lev > 0 && (n_field_gather_buffer > 0 || n_current_deposition_buffer > 0)) {
            for (int idim=0; idim < 3; ++idim)
            {
                if (Bfield_cax[lev][idim])
                {
                    const IntVect& ng = Bfield_cax[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(Bfield_cax[lev][idim]->boxArray(),
                                                                      dm, Bfield_cax[lev][idim]->nComp(), ng);
                    // pmf->ParallelCopy(*Bfield_cax[lev][idim], 0, 0, Bfield_cax[lev][idim]->nComp(), ng, ng);
                    Bfield_cax[lev][idim] = std::move(pmf);
                }
                if (Efield_cax[lev][idim])
                {
                    const IntVect& ng = Efield_cax[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(Efield_cax[lev][idim]->boxArray(),
                                                                      dm, Efield_cax[lev][idim]->nComp(), ng);
                    // pmf->ParallelCopy(*Efield_cax[lev][idim], 0, 0, Efield_cax[lev][idim]->nComp(), ng, ng);
                    Efield_cax[lev][idim] = std::move(pmf);
                }
                if (current_buf[lev][idim])
                {
                    const IntVect& ng = current_buf[lev][idim]->nGrowVect();
                    auto pmf = std::make_unique<MultiFab>(current_buf[lev][idim]->boxArray(),
                                                                      dm, current_buf[lev][idim]->nComp(), ng);
                    // pmf->ParallelCopy(*current_buf[lev][idim], 0, 0, current_buf[lev][idim]->nComp(), ng, ng);
                    current_buf[lev][idim] = std::move(pmf);
                }
            }
            if (charge_buf[lev])
            {
                const IntVect& ng = charge_buf[lev]->nGrowVect();
                auto pmf = std::make_unique<MultiFab>(charge_buf[lev]->boxArray(),
                                                                  dm, charge_buf[lev]->nComp(), ng);
                // pmf->ParallelCopy(*charge_buf[lev][idim], 0, 0, charge_buf[lev]->nComp(), ng, ng);
                charge_buf[lev] = std::move(pmf);
            }
            if (current_buffer_masks[lev])
            {
                const IntVect& ng = current_buffer_masks[lev]->nGrowVect();
                auto pmf = std::make_unique<iMultiFab>(current_buffer_masks[lev]->boxArray(),
                                                                    dm, current_buffer_masks[lev]->nComp(), ng);
                // we can avoid this since we immediately re-build the values via BuildBufferMasks()
                // pmf->Redistribute(*current_buffer_masks[lev], 0, 0, current_buffer_masks[lev]->nComp(), ng);
                current_buffer_masks[lev] = std::move(pmf);
            }
            if (gather_buffer_masks[lev])
            {
                const IntVect& ng = gather_buffer_masks[lev]->nGrowVect();
                auto pmf = std::make_unique<iMultiFab>(gather_buffer_masks[lev]->boxArray(),
                                                                    dm, gather_buffer_masks[lev]->nComp(), ng);
                // we can avoid this since we immediately re-build the values via BuildBufferMasks()
                // pmf->Redistribute(*gather_buffer_masks[lev], 0, 0, gather_buffer_masks[lev]->nComp(), ng);
                gather_buffer_masks[lev] = std::move(pmf);
            }
            if (current_buffer_masks[lev] || gather_buffer_masks[lev])
                BuildBufferMasks();
        }

        if (costs[lev] != nullptr)
        {
            costs[lev] = std::make_unique<LayoutData<Real>>(ba, dm);
            for (int i : costs[lev]->IndexArray())
            {
                (*costs[lev])[i] = 0.0;
                setLoadBalanceEfficiency(lev, -1);
            }
        }

        SetDistributionMap(lev, dm);

    } else
    {
        amrex::Abort("RemakeLevel: to be implemented");
    }
    // Re-initialize diagnostic functors that stores pointers to the user-requested fields at level, lev.
    multi_diags->InitializeFieldFunctors( lev );
}

void
WarpX::ComputeCostsHeuristic (amrex::Vector<std::unique_ptr<amrex::LayoutData<amrex::Real> > >& a_costs)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto & mypc_ref = GetInstance().GetPartContainer();
        const auto nSpecies = mypc_ref.nSpecies();

        // Species loop
        for (int i_s = 0; i_s < nSpecies; ++i_s)
        {
            auto & myspc = mypc_ref.GetParticleContainer(i_s);

            // Particle loop
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
            {
                (*a_costs[lev])[pti.index()] += costs_heuristic_particles_wt*pti.numParticles();
            }
        }

        //Cell loop
        MultiFab* Ex = Efield_fp[lev][0].get();
        for (MFIter mfi(*Ex, false); mfi.isValid(); ++mfi)
        {
            const Box& gbx = mfi.growntilebox();
            (*a_costs[lev])[mfi.index()] += costs_heuristic_cells_wt*gbx.numPts();
        }
    }
}

void
WarpX::ResetCosts ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int i : costs[lev]->IndexArray())
        {
            // Reset costs
            (*costs[lev])[i] = 0.0;
        }
    }
}
