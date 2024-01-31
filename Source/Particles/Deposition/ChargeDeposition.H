/* Copyright 2019 Axel Huebl, Andrew Myers, David Grote, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CHARGEDEPOSITION_H_
#define CHARGEDEPOSITION_H_

#include "Particles/Deposition/SharedDepositionUtils.H"
#include "ablastr/parallelization/KernelTimer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"
#include "Utils/WarpXAlgorithmSelection.H"
#ifdef WARPX_DIM_RZ
#   include "Utils/WarpX_Complex.H"
#endif

#include <AMReX.H>

/* \brief Perform charge deposition on a tile
 * \param GetPosition A functor for returning the particle position.
 * \param wp           Pointer to array of particle weights.
 * \param ion_lev      Pointer to array of particle ionization level. This is
                       required to have the charge of each macroparticle
                       since q is a scalar. For non-ionizable species,
                       ion_lev is a null pointer.
 * \param rho_fab      FArrayBox of charge density, either full array or tile.
 * \param np_to_deposit Number of particles for which current is deposited.
 * \param dx           3D cell size
 * \param xyzmin       Physical lower bounds of domain.
 * \param lo           Index lower bounds of domain.
 * \param q            species charge.
 * \param n_rz_azimuthal_modes Number of azimuthal modes when using RZ geometry.
 * \param cost: Pointer to (load balancing) cost corresponding to box where present particles deposit current.
 * \param load_balance_costs_update_algo Selected method for updating load balance costs.
 */
template <int depos_order>
void doChargeDepositionShapeN (const GetParticlePosition<PIdx>& GetPosition,
                               const amrex::ParticleReal * const wp,
                               const int* ion_lev,
                               amrex::FArrayBox& rho_fab,
                               long np_to_deposit,
                               const std::array<amrex::Real,3>& dx,
                               const std::array<amrex::Real, 3> xyzmin,
                               amrex::Dim3 lo,
                               amrex::Real q,
                               int n_rz_azimuthal_modes,
                               amrex::Real* cost,
                               long load_balance_costs_update_algo)
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

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_3D)
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
            np_to_deposit,
            [=] AMREX_GPU_DEVICE (long ip) {
#if defined(WARPX_USE_GPUCLOCK)
            const auto KernelTimer = ablastr::parallelization::KernelTimer(
                cost && (load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::GpuClock),
                cost_real);
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
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_3D)
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
#endif //defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_3D)
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
    if (cost && load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::GpuClock) {
        amrex::Gpu::streamSynchronize();
        *cost += *cost_real;
        amrex::The_Managed_Arena()->free(cost_real);
    }
#endif

#ifndef WARPX_DIM_RZ
    amrex::ignore_unused(n_rz_azimuthal_modes);
#endif
}

/* \brief Perform charge deposition on a tile using shared memory
 * \param GetPosition   A functor for returning the particle position.
 * \param wp            Pointer to array of particle weights.
 * \param ion_lev       Pointer to array of particle ionization level. This is
                        required to have the charge of each macroparticle
                        since q is a scalar. For non-ionizable species,
                        ion_lev is a null pointer.
 * \param rho_fab       FArrayBox of charge density, either full array or tile.
 * \param ix_type
 * \param np_to_deposit Number of particles for which charge is deposited.
 * \param dx            3D cell size
 * \param xyzmin        Physical lower bounds of domain.
 * \param lo            Index lower bounds of domain.
 * \param q             species charge.
 * \param n_rz_azimuthal_modes Number of azimuthal modes when using RZ geometry.
 * \param cost Pointer to (load balancing) cost corresponding to box where present particles deposit current.
 * \param load_balance_costs_update_algo Selected method for updating load balance costs.
 * \param a_bins
 * \param box
 * \param geom
 * \param a_tbox_max_size
 * \param bin_size tile size to use for shared current deposition operations
 */
template <int depos_order>
void doChargeDepositionSharedShapeN (const GetParticlePosition<PIdx>& GetPosition,
                                     const amrex::ParticleReal * const wp,
                                     const int* ion_lev,
                                     amrex::FArrayBox& rho_fab,
                                     const amrex::IntVect& ix_type,
                                     const long np_to_deposit,
                                     const std::array<amrex::Real,3>& dx,
                                     const std::array<amrex::Real, 3> xyzmin,
                                     const amrex::Dim3 lo,
                                     const amrex::Real q,
                                     const int n_rz_azimuthal_modes,
                                     amrex::Real* cost,
                                     const long load_balance_costs_update_algo,
                                     const amrex::DenseBins<WarpXParticleContainer::ParticleType>& a_bins,
                                     const amrex::Box& box,
                                     const amrex::Geometry& geom,
                                     const amrex::IntVect& a_tbox_max_size,
                                     const amrex::IntVect bin_size
                                     )
{
    using namespace amrex;

    const auto *permutation = a_bins.permutationPtr();

#if !defined(AMREX_USE_GPU)
    amrex::ignore_unused(ix_type, cost, load_balance_costs_update_algo, a_bins, box, geom, a_tbox_max_size, bin_size);
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

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_3D)
    const amrex::Real xmin = xyzmin[0];
#endif
#if defined(WARPX_DIM_3D)
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real zmin = xyzmin[2];

    amrex::Array4<amrex::Real> const& rho_arr = rho_fab.array();
    auto rho_box = rho_fab.box();
    amrex::IntVect const rho_type = rho_box.type();

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

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    const auto dxiarr = geom.InvCellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    const amrex::Box sample_tbox(IntVect(AMREX_D_DECL(0, 0, 0)), a_tbox_max_size - 1);
    amrex::Box sample_tbox_x = convert(sample_tbox, ix_type);

    sample_tbox_x.grow(depos_order);

    const auto npts = sample_tbox_x.numPts();

    const int nblocks = a_bins.numBins();
    const auto offsets_ptr = a_bins.offsetsPtr();
    const int threads_per_block = 256;

    std::size_t shared_mem_bytes = npts*sizeof(amrex::Real);

    const std::size_t max_shared_mem_bytes = amrex::Gpu::Device::sharedMemPerBlock();

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(shared_mem_bytes <= max_shared_mem_bytes,
                                     "Tile size too big for GPU shared memory charge deposition");
#endif

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    amrex::ignore_unused(np_to_deposit);
    // Loop with one block per tile (the shared memory is allocated on a per-block basis)
    // The threads within each block loop over the particles of its tile
    // (Each threads processes a different set of particles.)
    amrex::launch(
                  nblocks, threads_per_block, shared_mem_bytes, amrex::Gpu::gpuStream(),
                  [=] AMREX_GPU_DEVICE () noexcept
#else // defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    amrex::ParallelFor(np_to_deposit, [=] AMREX_GPU_DEVICE (long ip_orig) noexcept
#endif
      {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        const int bin_id = blockIdx.x;
        const unsigned int bin_start = offsets_ptr[bin_id];
        const unsigned int bin_stop = offsets_ptr[bin_id+1];

        if (bin_start == bin_stop) { return; }

        amrex::Box buffer_box;
        {
          ParticleReal xp, yp, zp;
          GetPosition(permutation[bin_start], xp, yp, zp);
#if defined(WARPX_DIM_3D)
          IntVect iv = IntVect(int(amrex::Math::floor((xp-plo[0])*dxiarr[0])),
                               int(amrex::Math::floor((yp-plo[1])*dxiarr[1])),
                               int(amrex::Math::floor((zp-plo[2])*dxiarr[2])));
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
          IntVect iv = IntVect(
                               int(amrex::Math::floor((xp-plo[0])*dxiarr[0])),
                               int(amrex::Math::floor((zp-plo[1])*dxiarr[1])));
#elif defined(WARPX_DIM_1D_Z)
          IntVect iv = IntVect(int(amrex::Math::floor((zp-plo[0])*dxiarr[0])));
#endif
          iv += domain.smallEnd();
          getTileIndex(iv, box, true, bin_size, buffer_box);
        }

        Box tbx = convert( buffer_box, ix_type);
        tbx.grow(depos_order);

        Gpu::SharedMemory<amrex::Real> gsm;
        amrex::Real* const shared = gsm.dataPtr();

        amrex::Array4<amrex::Real> buf(shared, amrex::begin(tbx), amrex::end(tbx), 1);

        // Zero-initialize the temporary array in shared memory
        volatile amrex::Real* vs = shared;
        for (int i = threadIdx.x; i < tbx.numPts(); i += blockDim.x) {
          vs[i] = 0.0;
        }
        __syncthreads();
#else
        amrex::Array4<amrex::Real> const &buf = rho_arr;
#endif // defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        // Loop over macroparticles: each threads loops over particles with a stride of `blockDim.x`
        for (unsigned int ip_orig = bin_start + threadIdx.x; ip_orig < bin_stop; ip_orig += blockDim.x)
#endif
        {
            const unsigned int ip = permutation[ip_orig];

#if defined(WARPX_USE_GPUCLOCK)
            const auto KernelTimer = ablastr::parallelization::KernelTimer(
                cost && (load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::GpuClock),
                cost_real);
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
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_3D)
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
#endif //defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_3D)
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

            // Deposit charge into buf
#if defined(WARPX_DIM_1D_Z)
            for (int iz=0; iz<=depos_order; iz++){
                amrex::Gpu::Atomic::AddNoRet(
                    &buf(lo.x+k+iz, 0, 0, 0),
                    sz[iz]*wq);
            }
#endif
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    amrex::Gpu::Atomic::AddNoRet(
                        &buf(lo.x+i+ix, lo.y+k+iz, 0, 0),
                        sx[ix]*sz[iz]*wq);
#if defined(WARPX_DIM_RZ)
                    Complex xy = xy0; // Throughout the following loop, xy takes the value e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 on the weighting comes from the normalization of the modes
                        amrex::Gpu::Atomic::AddNoRet( &buf(lo.x+i+ix, lo.y+k+iz, 0, 2*imode-1), 2._rt*sx[ix]*sz[iz]*wq*xy.real());
                        amrex::Gpu::Atomic::AddNoRet( &buf(lo.x+i+ix, lo.y+k+iz, 0, 2*imode  ), 2._rt*sx[ix]*sz[iz]*wq*xy.imag());
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
                            &buf(lo.x+i+ix, lo.y+j+iy, lo.z+k+iz),
                            sx[ix]*sy[iy]*sz[iz]*wq);
                    }
                }
            }
#endif
        }

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
          __syncthreads();

          addLocalToGlobal(tbx, rho_arr, buf);
#endif // defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
      }
      );
#if defined(WARPX_USE_GPUCLOCK)
        if(cost && load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::GpuClock) {
            amrex::Gpu::streamSynchronize();
            *cost += *cost_real;
            amrex::The_Managed_Arena()->free(cost_real);
        }
#endif

#ifndef WARPX_DIM_RZ
        amrex::ignore_unused(n_rz_azimuthal_modes);
#endif
}

#endif // CHARGEDEPOSITION_H_
