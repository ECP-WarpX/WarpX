/* Copyright 2019-2021 Axel Huebl, Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_DEPOSIT_CHARGE_H_
#define ABLASTR_DEPOSIT_CHARGE_H_

#include "ablastr/profiler/ProfilerWrapper.H"
#include "Parallelization/KernelTimer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"
#include "Particles/Deposition/ChargeDeposition.H"
#ifdef WARPX_DIM_RZ
#   include "Utils/WarpX_Complex.H"
#endif

#include <AMReX.H>


namespace ablastr {
namespace particles {

/** Perform charge deposition for the particles on a tile.
 *
 * \tparam T_PC a type of amrex::ParticleContainer
 *
 * \param pti an amrex::ParIter pointing to the tile to operate on
 * \param wp vector of the particle weights for those particles.
 * \param ion_lev pointer to array of particle ionization level. This is
                  required to have the charge of each macroparticle
                  since q is a scalar. For non-ionizable species,
                  ion_lev is a null pointer.
 * \param rho MultiFab of the charge density
 * \param icomp component in MultiFab to start depositing to
 * \param nc number of components to deposit
 * \param offset index to start at when looping over particles to depose
 * \param np_to_depose number of particles to depose
 * \param local_rho temporary FArrayBox for deposition with OpenMP
 * \param lev the level of the particles we are on
 * \param depos_lev the level to deposit the particles to
 * \param charge charge of the particle species
 * \param nox shape factor in the x direction
 * \param noy shape factor in the y direction
 * \param noz shape factor in the z direction
 * \param ng_rho number of ghost cells to use for rho
 * \param dx cell spacing at level lev
 * \param xyzmin lo corner of the current tile in physical coordinates.
 * \param ref_ratio mesh refinement ratio between lev and depos_lev
 * \param cost pointer to (load balancing) cost corresponding to box where present
               particles deposit current. If nullptr, costs are not updated.
 * \param n_rz_azimuthal_modes number of azimuthal modes in use, irrelevant outside RZ geometry.
 * \param load_balance_costs_update_algo selected method for updating load balance costs.
 * \param do_device_synchronize call amrex::Gpu::synchronize() for tiny profiler regions
 */
template< typename T_PC >
static void
deposit_charge (typename T_PC::ParIterType& pti,
                typename T_PC::RealVector& wp,
                const int * const ion_lev,
                amrex::MultiFab* rho, const int icomp, const int nc,
                const long offset, const long np_to_depose,
                amrex::FArrayBox& local_rho, const int lev, const int depos_lev,
                const amrex::Real charge, const int nox, const int noy, const int noz,
                const amrex::IntVect& ng_rho, const std::array<amrex::Real,3>& dx,
                const std::array<amrex::Real, 3>& xyzmin,
                const amrex::IntVect& ref_ratio,
                amrex::Real* cost, const int n_rz_azimuthal_modes,
                const long load_balance_costs_update_algo,
                const bool do_device_synchronize)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((depos_lev==(lev-1)) ||
                                     (depos_lev==(lev  )),
                                     "Deposition buffers only work for lev-1");

    // If no particles, do not do anything
    if (np_to_depose == 0) return;

    // Extract deposition order and check that particles shape fits within the guard cells.
    // NOTE: In specific situations where the staggering of rho and the charge deposition algorithm
    // are not trivial, this check might be too strict and we might need to relax it, as currently
    // done for the current deposition.

#if   defined(WARPX_DIM_1D_Z)
    amrex::ignore_unused(nox);
    amrex::ignore_unused(noy);
    const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(noz/2+1));
#elif   defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    amrex::ignore_unused(noy);
    const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(nox/2+1),
                                                       static_cast<int>(noz/2+1));
#elif defined(WARPX_DIM_3D)
    const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(nox/2+1),
                                                       static_cast<int>(noy/2+1),
                                                       static_cast<int>(noz/2+1));
#endif

    // On CPU: particles deposit on tile arrays, which have a small number of guard cells ng_rho
    // On GPU: particles deposit directly on the rho array, which usually have a larger number of guard cells
#ifndef AMREX_USE_GPU
    const amrex::IntVect range = ng_rho - shape_extent;
#else
    const amrex::IntVect range = rho->nGrowVect() - shape_extent;
#endif
    amrex::ignore_unused(range); // for release builds
    AMREX_ASSERT_WITH_MESSAGE(
        amrex::numParticlesOutOfRange(pti, range) == 0,
        "Particles shape does not fit within tile (CPU) or guard cells (GPU) used for charge deposition");

    ABLASTR_PROFILE_VAR_NS("ablastr::particles::deposit_charge::ChargeDeposition", blp_ppc_chd, do_device_synchronize);
    ABLASTR_PROFILE_VAR_NS("ablastr::particles::deposit_charge::Accumulate", blp_accumulate, do_device_synchronize);

    // Get tile box where charge is deposited.
    // The tile box is different when depositing in the buffers (depos_lev<lev)
    // or when depositing inside the level (depos_lev=lev)
    amrex::Box tilebox;
    if (lev == depos_lev) {
        tilebox = pti.tilebox();
    } else {
        tilebox = amrex::coarsen(pti.tilebox(),ref_ratio);
    }

#ifndef AMREX_USE_GPU
    // Staggered tile box
    amrex::Box tb = amrex::convert( tilebox, rho->ixType().toIntVect() );
#endif

    tilebox.grow(ng_rho);

#ifdef AMREX_USE_GPU
    amrex::ignore_unused(local_rho);
    // GPU, no tiling: rho_fab points to the full rho array
    amrex::MultiFab rhoi(*rho, amrex::make_alias, icomp*nc, nc);
    auto & rho_fab = rhoi.get(pti);
#else
    tb.grow(ng_rho);

    // CPU, tiling: rho_fab points to local_rho
    local_rho.resize(tb, nc);

    // local_rho is set to zero
    local_rho.setVal(0.0);

    auto & rho_fab = local_rho;
#endif

    const auto GetPosition = GetParticlePosition(pti, offset);

    // Indices of the lower bound
    const amrex::Dim3 lo = lbound(tilebox);

    ABLASTR_PROFILE_VAR_START(blp_ppc_chd, do_device_synchronize);

    if        (nox == 1){
        doChargeDepositionShapeN<1>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                    rho_fab, np_to_depose, dx, xyzmin, lo, charge,
                                    n_rz_azimuthal_modes, cost,
                                    load_balance_costs_update_algo);
    } else if (nox == 2){
        doChargeDepositionShapeN<2>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                    rho_fab, np_to_depose, dx, xyzmin, lo, charge,
                                    n_rz_azimuthal_modes, cost,
                                    load_balance_costs_update_algo);
    } else if (nox == 3){
        doChargeDepositionShapeN<3>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                    rho_fab, np_to_depose, dx, xyzmin, lo, charge,
                                    n_rz_azimuthal_modes, cost,
                                    load_balance_costs_update_algo);
    }
    ABLASTR_PROFILE_VAR_STOP(blp_ppc_chd, do_device_synchronize);

#ifndef AMREX_USE_GPU
    // CPU, tiling: atomicAdd local_rho into rho
    ABLASTR_PROFILE_VAR_START(blp_accumulate, do_device_synchronize);
    (*rho)[pti].atomicAdd(local_rho, tb, tb, 0, icomp*nc, nc);
    ABLASTR_PROFILE_VAR_STOP(blp_accumulate, do_device_synchronize);
#endif
}

} // namespace particles
} // namespace ablastr

#endif // ABLASTR_DEPOSIT_CHARGE_H_
