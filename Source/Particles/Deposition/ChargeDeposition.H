/* Copyright 2019 Axel Huebl, Andrew Myers, David Grote, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CHARGEDEPOSITION_H_
#define CHARGEDEPOSITION_H_

#include "Parallelization/KernelTimer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXProfilerWrapper.H"
#ifdef WARPX_DIM_RZ
#   include "Utils/WarpX_Complex.H"
#endif

#include <AMReX.H>

/* \brief Charge Deposition for thread thread_num
 * \param GetPosition : A functor for returning the particle position.
 * \param wp           : Pointer to array of particle weights.
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param rho_fab      : FArrayBox of charge density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param q            : species charge.
 * \param n_rz_azimuthal_modes: Number of azimuthal modes when using RZ geometry.
 * \param cost: Pointer to (load balancing) cost corresponding to box where present particles deposit current.
 * \param load_balance_costs_update_algo: Selected method for updating load balance costs.
 */
template <int depos_order>
void doChargeDepositionShapeN (const GetParticlePosition& GetPosition,
                               const amrex::ParticleReal * const wp,
                               const int * const ion_lev,
                               amrex::FArrayBox& rho_fab,
                               const long np_to_depose,
                               const std::array<amrex::Real,3>& dx,
                               const std::array<amrex::Real, 3> xyzmin,
                               const amrex::Dim3 lo,
                               const amrex::Real q,
                               const int n_rz_azimuthal_modes,
                               amrex::Real* cost,
                               const long load_balance_costs_update_algo)
{
    using namespace amrex;

#if !defined(AMREX_USE_GPU)
    amrex::ignore_unused(cost, load_balance_costs_update_algo);
#endif

    // Whether ion_lev is a null pointer (do_ionization=0) or a real pointer
    // (do_ionization=1)
    const bool do_ionization = ion_lev;
    const amrex::Real dzi = 1.0_rt/dx[2];
#if defined(WARPX_DIM_1D_Z)
    const amrex::Real invvol = dzi;
#endif
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const amrex::Real dxi = 1.0_rt/dx[0];
    const amrex::Real invvol = dxi*dzi;
#elif defined(WARPX_DIM_3D)
    const amrex::Real dxi = 1.0_rt/dx[0];
    const amrex::Real dyi = 1.0_rt/dx[1];
    const amrex::Real invvol = dxi*dyi*dzi;
#endif

#if (AMREX_SPACEDIM >= 2)
    const amrex::Real xmin = xyzmin[0];
#endif
#if defined(WARPX_DIM_3D)
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real zmin = xyzmin[2];

    amrex::Array4<amrex::Real> const& rho_arr = rho_fab.array();
    amrex::IntVect const rho_type = rho_fab.box().type();

    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // Loop over particles and deposit into rho_fab
#if defined(WARPX_USE_GPUCLOCK)
    amrex::Real* cost_real = nullptr;
    if( load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::GpuClock) {
        cost_real = (amrex::Real *) amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real));
        *cost_real = 0.;
    }
#endif
    amrex::ParallelFor(
        np_to_depose,
        [=] AMREX_GPU_DEVICE (long ip) {
#if defined(WARPX_USE_GPUCLOCK)
            KernelTimer kernelTimer(cost && load_balance_costs_update_algo
                                 == LoadBalanceCostsUpdateAlgo::GpuClock, cost_real);
#endif
            // --- Get particle quantities
            amrex::Real wq = q*wp[ip]*invvol;
            if (do_ionization){
                wq *= ion_lev[ip];
            }

            amrex::ParticleReal xp, yp, zp;
            GetPosition(ip, xp, yp, zp);

            // --- Compute shape factors
            Compute_shape_factor< depos_order > const compute_shape_factor;
#if (AMREX_SPACEDIM >= 2)
            // x direction
            // Get particle position in grid coordinates
#if defined(WARPX_DIM_RZ)
            const amrex::Real rp = std::sqrt(xp*xp + yp*yp);
            amrex::Real costheta;
            amrex::Real sintheta;
            if (rp > 0.) {
                costheta = xp/rp;
                sintheta = yp/rp;
            } else {
                costheta = 1._rt;
                sintheta = 0._rt;
            }
            const Complex xy0 = Complex{costheta, sintheta};
            const amrex::Real x = (rp - xmin)*dxi;
#else
            const amrex::Real x = (xp - xmin)*dxi;
#endif

            // Compute shape factor along x
            // i: leftmost grid point that the particle touches
            amrex::Real sx[depos_order + 1] = {0._rt};
            int i = 0;
            if (rho_type[0] == NODE) {
                i = compute_shape_factor(sx, x);
            } else if (rho_type[0] == CELL) {
                i = compute_shape_factor(sx, x - 0.5_rt);
            }
#endif //AMREX_SPACEDIM >= 2
#if defined(WARPX_DIM_3D)
            // y direction
            const amrex::Real y = (yp - ymin)*dyi;
            amrex::Real sy[depos_order + 1] = {0._rt};
            int j = 0;
            if (rho_type[1] == NODE) {
                j = compute_shape_factor(sy, y);
            } else if (rho_type[1] == CELL) {
                j = compute_shape_factor(sy, y - 0.5_rt);
            }
#endif
            // z direction
            const amrex::Real z = (zp - zmin)*dzi;
            amrex::Real sz[depos_order + 1] = {0._rt};
            int k = 0;
            if (rho_type[WARPX_ZINDEX] == NODE) {
                k = compute_shape_factor(sz, z);
            } else if (rho_type[WARPX_ZINDEX] == CELL) {
                k = compute_shape_factor(sz, z - 0.5_rt);
            }

            // Deposit charge into rho_arr
#if defined(WARPX_DIM_1D_Z)
            for (int iz=0; iz<=depos_order; iz++){
                amrex::Gpu::Atomic::AddNoRet(
                    &rho_arr(lo.x+k+iz, 0, 0, 0),
                    sz[iz]*wq);
            }
#endif
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    amrex::Gpu::Atomic::AddNoRet(
                        &rho_arr(lo.x+i+ix, lo.y+k+iz, 0, 0),
                        sx[ix]*sz[iz]*wq);
#if defined(WARPX_DIM_RZ)
                    Complex xy = xy0; // Throughout the following loop, xy takes the value e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 on the weighting comes from the normalization of the modes
                        amrex::Gpu::Atomic::AddNoRet( &rho_arr(lo.x+i+ix, lo.y+k+iz, 0, 2*imode-1), 2._rt*sx[ix]*sz[iz]*wq*xy.real());
                        amrex::Gpu::Atomic::AddNoRet( &rho_arr(lo.x+i+ix, lo.y+k+iz, 0, 2*imode  ), 2._rt*sx[ix]*sz[iz]*wq*xy.imag());
                        xy = xy*xy0;
                    }
#endif
                }
            }
#elif defined(WARPX_DIM_3D)
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        amrex::Gpu::Atomic::AddNoRet(
                            &rho_arr(lo.x+i+ix, lo.y+j+iy, lo.z+k+iz),
                            sx[ix]*sy[iy]*sz[iz]*wq);
                    }
                }
            }
#endif
        }
        );
#if defined(WARPX_USE_GPUCLOCK)
        if( load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::GpuClock) {
            amrex::Gpu::streamSynchronize();
            *cost += *cost_real;
            amrex::The_Managed_Arena()->free(cost_real);
        }
#endif

#ifndef WARPX_DIM_RZ
        amrex::ignore_unused(n_rz_azimuthal_modes);
#endif
}

template <typename PC>
void DepositChargeDoIt (typename PC::ParIterType& pti,
                        typename PC::RealVector& wp,
                        const int * const ion_lev,
                        amrex::MultiFab* rho, const int icomp,
                        const long offset, const long np_to_depose,
                        amrex::FArrayBox& local_rho, const int lev, const int depos_lev,
                        const amrex::Real charge, const int nox, const int noy, const int noz,
                        const amrex::IntVect& ng_rho, const std::array<amrex::Real,3>& dx,
                        const amrex::IntVect& ref_ratio, const int nc,
                        const std::array<amrex::Real, 3>& xyzmin,
                        amrex::Real* cost, const int n_rz_azimuthal_modes,
                        const long load_balance_costs_update_algo)
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
    const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(noz/2+1));
#elif   defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::numParticlesOutOfRange(pti, range) == 0,
        "Particles shape does not fit within tile (CPU) or guard cells (GPU) used for charge deposition");

    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCharge::ChargeDeposition", blp_ppc_chd);
    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCharge::Accumulate", blp_accumulate);

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

    WARPX_PROFILE_VAR_START(blp_ppc_chd);

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
    WARPX_PROFILE_VAR_STOP(blp_ppc_chd);

#ifndef AMREX_USE_GPU
    // CPU, tiling: atomicAdd local_rho into rho
    WARPX_PROFILE_VAR_START(blp_accumulate);
    (*rho)[pti].atomicAdd(local_rho, tb, tb, 0, icomp*nc, nc);
    WARPX_PROFILE_VAR_STOP(blp_accumulate);
#endif
}

#endif // CHARGEDEPOSITION_H_
