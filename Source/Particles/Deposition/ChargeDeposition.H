/* Copyright 2019 Axel Huebl, David Grote, Maxence Thevenet
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
#include "Utils/FiniteCentering.H"
#include "Utils/WarpXAlgorithmSelection.H"
#ifdef WARPX_DIM_RZ
#   include "Utils/WarpX_Complex.H"
#endif

#include <AMReX.H>

using namespace amrex::literals;

// TODO Add Doxygen
AMREX_FORCE_INLINE
void doChargeDepositionFiniteCentering (const GetParticlePosition& GetPosition,
                                        const amrex::ParticleReal * const wp,
                                        const int * const ion_lev,
                                        amrex::FArrayBox& rho_fab,
                                        const long np_to_depose,
                                        const std::array<amrex::Real,3>& dxyz,
                                        const std::array<amrex::Real,3>& xyz_min,
                                        const amrex::Dim3 lo,
                                        const amrex::Real q,
                                        const int nox,
                                        const int noy,
                                        const int noz)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(noy);
#endif

    // Whether ion_lev is a null pointer (do_ionization = 0) or a real pointer (do_ionization = 1)
    const bool do_ionization = ion_lev;

    const amrex::Real dxi = 1.0_rt / dxyz[0];
    const amrex::Real dzi = 1.0_rt / dxyz[2];
#if defined(WARPX_DIM_XZ)
    const amrex::Real invvol = dxi * dzi;
#elif defined(WARPX_DIM_3D)
    const amrex::Real dyi = 1.0_rt / dxyz[1];
    const amrex::Real invvol = dxi * dyi * dzi;
#endif

    amrex::Array4<amrex::Real> const& rho_arr = rho_fab.array();

    amrex::ParallelFor(np_to_depose, [=] AMREX_GPU_DEVICE (long idxp)
    {
        amrex::Real wq = (do_ionization) ? q * wp[idxp] * ion_lev[idxp] * invvol : q * wp[idxp] * invvol;

        amrex::ParticleReal xp, yp, zp;
        GetPosition(idxp, xp, yp, zp);

        const amrex::Real x = ((xp - xyz_min[0])) * dxi;
#if (AMREX_SPACEDIM == 3)
        const amrex::Real y = ((yp - xyz_min[1])) * dyi;
#endif
        const amrex::Real z = ((zp - xyz_min[2])) * dzi;

        // Compute largest nearest integer index
        const int ip = static_cast<int>(std::ceil(x));
#if   (AMREX_SPACEDIM == 2)
        const int jp = static_cast<int>(std::ceil(z));
        const int kp = 0;
#elif (AMREX_SPACEDIM == 3)
        const int jp = static_cast<int>(std::ceil(y));
        const int kp = static_cast<int>(std::ceil(z));
#endif

        const int noi = nox;
#if   (AMREX_SPACEDIM == 2)
        const int noj = noz;
#elif (AMREX_SPACEDIM == 3)
        const int noj = noy;
        const int nok = noz;
#endif

        // Min and max for interpolation loop along i
        const int imin = ip - noi/2;
        const int imax = ip + noi/2 - 1;

        // Min and max for interpolation loop along j
        const int jmin = jp - noj/2;
        const int jmax = jp + noj/2 - 1;

        // Min and max for interpolation loop along k
#if   (AMREX_SPACEDIM == 2)
        // kp = 0 always
        const int kmin = kp;
        const int kmax = kp;
#elif (AMREX_SPACEDIM == 3)
        const int kmin = kp - nok/2;
        const int kmax = kp + nok/2 - 1;
#endif

        // Number of interpolation points
        const int ni = imax - imin;
        const int nj = jmax - jmin;
        const int nk = kmax - kmin;

        // Compute difference r in each direction
        // (r always negative, by construction)
        const amrex::Real ri = x - ip;
#if   (AMREX_SPACEDIM == 2)
        const amrex::Real rj = z - jp;
#elif (AMREX_SPACEDIM == 3)
        const amrex::Real rj = y - jp;
        const amrex::Real rk = z - kp;
#endif

        // Compute centering coefficients
        // TODO Currently using maximum size to fix GPU builds, needs to be optimized
        amrex::GpuArray<amrex::Real, 8> sci;
        amrex::GpuArray<amrex::Real, 8> scj;
#if (AMREX_SPACEDIM == 3)
        amrex::GpuArray<amrex::Real, 8> sck;
#endif

        FiniteCentering::compute_centering_coefficients(sci, ri, noi);
        FiniteCentering::compute_centering_coefficients(scj, rj, noj);
#if (AMREX_SPACEDIM == 3)
        FiniteCentering::compute_centering_coefficients(sck, rk, nok);
#endif

        amrex::Real ci = 1.0_rt;
        amrex::Real cj = 1.0_rt;
        amrex::Real ck = 1.0_rt;

        // TODO Find how to do this better
        bool not_out_of_bound = false;
        const amrex::Dim3 lo_rho = lbound(rho_arr);
        const amrex::Dim3 hi_rho = ubound(rho_arr);

        int i, j, k;
        for (int kk = 0; kk <= nk; kk++)
        {
            k = lo.z + kmin + kk;
#if (AMREX_SPACEDIM == 3)
            ck = sck[kk];
#endif
            // TODO Find how to do this better
            not_out_of_bound = (k >= lo_rho.z && k <= hi_rho.z);

            for (int jj = 0; jj <= nj; jj++)
            {
                j = lo.y + jmin + jj;
                cj = scj[jj];

                // TODO Find how to do this better
                not_out_of_bound = not_out_of_bound && (j >= lo_rho.y && j <= hi_rho.y);

                for (int ii = 0; ii <= ni; ii++)
                {
                    i = lo.x + imin + ii;
                    ci = sci[ii];

                    // TODO Find how to do this better
                    not_out_of_bound =  not_out_of_bound && (i >= lo_rho.x && i <= hi_rho.x);

                    if (not_out_of_bound)
                    {
                        amrex::Gpu::Atomic::AddNoRet(&rho_arr(lo.x+imin+ii,lo.y+jmin+jj,lo.z+kmin+kk),
                                                     ci * cj * ck * wq);
                    }
                }
            }
        }
    });
}

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

#endif // CHARGEDEPOSITION_H_
