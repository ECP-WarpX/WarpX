/* Copyright 2019 Axel Huebl, David Grote, Maxence Thevenet
 * Remi Lehe, Weiqun Zhang, Michael Rowan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"
#include "Utils/WarpX_Complex.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

/**
 * \brief Current Deposition for thread thread_num
 * /param GetPosition : A functor for returning the particle position.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param jx_fab       : FArrayBox of current density, either full array or tile.
 * \param jy_fab       : FArrayBox of current density, either full array or tile.
 * \param jz_fab       : FArrayBox of current density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dt           : Time step for particle level
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * /param q            : species charge.
 */
template <int depos_order>
void doDepositionShapeN(const GetParticlePosition& GetPosition,
                        const amrex::ParticleReal * const wp,
                        const amrex::ParticleReal * const uxp,
                        const amrex::ParticleReal * const uyp,
                        const amrex::ParticleReal * const uzp,
                        const int * const ion_lev,
                        amrex::FArrayBox& jx_fab,
                        amrex::FArrayBox& jy_fab,
                        amrex::FArrayBox& jz_fab,
                        const long np_to_depose, const amrex::Real dt,
                        const std::array<amrex::Real,3>& dx,
                        const std::array<amrex::Real,3>& xyzmin,
                        const amrex::Dim3 lo,
                        const amrex::Real q,
                        const int n_rz_azimuthal_modes)
{
#if !defined(WARPX_DIM_RZ)
    amrex::ignore_unused(n_rz_azimuthal_modes);
#endif

    // Whether ion_lev is a null pointer (do_ionization=0) or a real pointer
    // (do_ionization=1)
    const bool do_ionization = ion_lev;
    const amrex::Real dxi = 1.0_rt/dx[0];
    const amrex::Real dzi = 1.0_rt/dx[2];
#if !(defined WARPX_DIM_RZ)
    const amrex::Real dts2dx = 0.5_rt*dt*dxi;
#endif
    const amrex::Real dts2dz = 0.5_rt*dt*dzi;
#if (AMREX_SPACEDIM == 2)
    const amrex::Real invvol = dxi*dzi;
#elif (defined WARPX_DIM_3D)
    const amrex::Real dyi = 1.0_rt/dx[1];
    const amrex::Real dts2dy = 0.5_rt*dt*dyi;
    const amrex::Real invvol = dxi*dyi*dzi;
#endif

    const amrex::Real xmin = xyzmin[0];
#if (defined WARPX_DIM_3D)
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real zmin = xyzmin[2];

    const amrex::Real clightsq = 1.0_rt/PhysConst::c/PhysConst::c;

    amrex::Array4<amrex::Real> const& jx_arr = jx_fab.array();
    amrex::Array4<amrex::Real> const& jy_arr = jy_fab.array();
    amrex::Array4<amrex::Real> const& jz_arr = jz_fab.array();
    amrex::IntVect const jx_type = jx_fab.box().type();
    amrex::IntVect const jy_type = jy_fab.box().type();
    amrex::IntVect const jz_type = jz_fab.box().type();

    constexpr int zdir = (AMREX_SPACEDIM - 1);
    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // Loop over particles and deposit into jx_fab, jy_fab and jz_fab
    amrex::ParallelFor(
        np_to_depose,
        [=] AMREX_GPU_DEVICE (long ip) {
            // --- Get particle quantities
            const amrex::Real gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
                                                         + uyp[ip]*uyp[ip]*clightsq
                                                         + uzp[ip]*uzp[ip]*clightsq);
            amrex::Real wq  = q*wp[ip];
            if (do_ionization){
                wq *= ion_lev[ip];
            }

            amrex::ParticleReal xp, yp, zp;
            GetPosition(ip, xp, yp, zp);

            const amrex::Real vx  = uxp[ip]*gaminv;
            const amrex::Real vy  = uyp[ip]*gaminv;
            const amrex::Real vz  = uzp[ip]*gaminv;
            // wqx, wqy wqz are particle current in each direction
#if (defined WARPX_DIM_RZ)
            // In RZ, wqx is actually wqr, and wqy is wqtheta
            // Convert to cylinderical at the mid point
            const amrex::Real xpmid = xp - 0.5_rt*dt*vx;
            const amrex::Real ypmid = yp - 0.5_rt*dt*vy;
            const amrex::Real rpmid = std::sqrt(xpmid*xpmid + ypmid*ypmid);
            amrex::Real costheta;
            amrex::Real sintheta;
            if (rpmid > 0.) {
                costheta = xpmid/rpmid;
                sintheta = ypmid/rpmid;
            } else {
                costheta = 1._rt;
                sintheta = 0._rt;
            }
            const Complex xy0 = Complex{costheta, sintheta};
            const amrex::Real wqx = wq*invvol*(+vx*costheta + vy*sintheta);
            const amrex::Real wqy = wq*invvol*(-vx*sintheta + vy*costheta);
#else
            const amrex::Real wqx = wq*invvol*vx;
            const amrex::Real wqy = wq*invvol*vy;
#endif
            const amrex::Real wqz = wq*invvol*vz;

            // --- Compute shape factors
            // x direction
            // Get particle position after 1/2 push back in position
#if (defined WARPX_DIM_RZ)
            // Keep these double to avoid bug in single precision
            const double xmid = (rpmid - xmin)*dxi;
#else
            const double xmid = (xp - xmin)*dxi - dts2dx*vx;
#endif
            // j_j[xyz] leftmost grid point in x that the particle touches for the centering of each current
            // sx_j[xyz] shape factor along x for the centering of each current
            // There are only two possible centerings, node or cell centered, so at most only two shape factor
            // arrays will be needed.
            // Keep these double to avoid bug in single precision
            double sx_node[depos_order + 1];
            double sx_cell[depos_order + 1];
            int j_node = 0;
            int j_cell = 0;
            Compute_shape_factor< depos_order > const compute_shape_factor;
            if (jx_type[0] == NODE || jy_type[0] == NODE || jz_type[0] == NODE) {
                j_node = compute_shape_factor(sx_node, xmid);
            }
            if (jx_type[0] == CELL || jy_type[0] == CELL || jz_type[0] == CELL) {
                j_cell = compute_shape_factor(sx_cell, xmid - 0.5);
            }

            amrex::Real sx_jx[depos_order + 1] = {0._rt};
            amrex::Real sx_jy[depos_order + 1] = {0._rt};
            amrex::Real sx_jz[depos_order + 1] = {0._rt};
            for (int ix=0; ix<=depos_order; ix++)
            {
                sx_jx[ix] = ((jx_type[0] == NODE) ? amrex::Real(sx_node[ix]) : amrex::Real(sx_cell[ix]));
                sx_jy[ix] = ((jy_type[0] == NODE) ? amrex::Real(sx_node[ix]) : amrex::Real(sx_cell[ix]));
                sx_jz[ix] = ((jz_type[0] == NODE) ? amrex::Real(sx_node[ix]) : amrex::Real(sx_cell[ix]));
            }

            int const j_jx = ((jx_type[0] == NODE) ? j_node : j_cell);
            int const j_jy = ((jy_type[0] == NODE) ? j_node : j_cell);
            int const j_jz = ((jz_type[0] == NODE) ? j_node : j_cell);

#if (defined WARPX_DIM_3D)
            // y direction
            // Keep these double to avoid bug in single precision
            const double ymid = (yp - ymin)*dyi - dts2dy*vy;
            double sy_node[depos_order + 1];
            double sy_cell[depos_order + 1];
            int k_node = 0;
            int k_cell = 0;
            if (jx_type[1] == NODE || jy_type[1] == NODE || jz_type[1] == NODE) {
                k_node = compute_shape_factor(sy_node, ymid);
            }
            if (jx_type[1] == CELL || jy_type[1] == CELL || jz_type[1] == CELL) {
                k_cell = compute_shape_factor(sy_cell, ymid - 0.5);
            }
            amrex::Real sy_jx[depos_order + 1] = {0.};
            amrex::Real sy_jy[depos_order + 1] = {0.};
            amrex::Real sy_jz[depos_order + 1] = {0.};
            for (int iy=0; iy<=depos_order; iy++)
            {
                sy_jx[iy] = ((jx_type[1] == NODE) ? amrex::Real(sy_node[iy]) : amrex::Real(sy_cell[iy]));
                sy_jy[iy] = ((jy_type[1] == NODE) ? amrex::Real(sy_node[iy]) : amrex::Real(sy_cell[iy]));
                sy_jz[iy] = ((jz_type[1] == NODE) ? amrex::Real(sy_node[iy]) : amrex::Real(sy_cell[iy]));
            }
            int const k_jx = ((jx_type[1] == NODE) ? k_node : k_cell);
            int const k_jy = ((jy_type[1] == NODE) ? k_node : k_cell);
            int const k_jz = ((jz_type[1] == NODE) ? k_node : k_cell);
#endif

            // z direction
            // Keep these double to avoid bug in single precision
            const double zmid = (zp - zmin)*dzi - dts2dz*vz;
            double sz_node[depos_order + 1];
            double sz_cell[depos_order + 1];
            int l_node = 0;
            int l_cell = 0;
            if (jx_type[zdir] == NODE || jy_type[zdir] == NODE || jz_type[zdir] == NODE) {
                l_node = compute_shape_factor(sz_node, zmid);
            }
            if (jx_type[zdir] == CELL || jy_type[zdir] == CELL || jz_type[zdir] == CELL) {
                l_cell = compute_shape_factor(sz_cell, zmid - 0.5);
            }
            amrex::Real sz_jx[depos_order + 1] = {0.};
            amrex::Real sz_jy[depos_order + 1] = {0.};
            amrex::Real sz_jz[depos_order + 1] = {0.};
            for (int iz=0; iz<=depos_order; iz++)
            {
                sz_jx[iz] = ((jx_type[zdir] == NODE) ? amrex::Real(sz_node[iz]) : amrex::Real(sz_cell[iz]));
                sz_jy[iz] = ((jy_type[zdir] == NODE) ? amrex::Real(sz_node[iz]) : amrex::Real(sz_cell[iz]));
                sz_jz[iz] = ((jz_type[zdir] == NODE) ? amrex::Real(sz_node[iz]) : amrex::Real(sz_cell[iz]));
            }
            int const l_jx = ((jx_type[zdir] == NODE) ? l_node : l_cell);
            int const l_jy = ((jy_type[zdir] == NODE) ? l_node : l_cell);
            int const l_jz = ((jz_type[zdir] == NODE) ? l_node : l_cell);

            // Deposit current into jx_arr, jy_arr and jz_arr
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    amrex::Gpu::Atomic::AddNoRet(
                        &jx_arr(lo.x+j_jx+ix, lo.y+l_jx+iz, 0, 0),
                        sx_jx[ix]*sz_jx[iz]*wqx);
                    amrex::Gpu::Atomic::AddNoRet(
                        &jy_arr(lo.x+j_jy+ix, lo.y+l_jy+iz, 0, 0),
                        sx_jy[ix]*sz_jy[iz]*wqy);
                    amrex::Gpu::Atomic::AddNoRet(
                        &jz_arr(lo.x+j_jz+ix, lo.y+l_jz+iz, 0, 0),
                        sx_jz[ix]*sz_jz[iz]*wqz);
#if (defined WARPX_DIM_RZ)
                    Complex xy = xy0; // Note that xy is equal to e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 on the weighting comes from the normalization of the modes
                        amrex::Gpu::Atomic::AddNoRet( &jx_arr(lo.x+j_jx+ix, lo.y+l_jx+iz, 0, 2*imode-1), 2._rt*sx_jx[ix]*sz_jx[iz]*wqx*xy.real());
                        amrex::Gpu::Atomic::AddNoRet( &jx_arr(lo.x+j_jx+ix, lo.y+l_jx+iz, 0, 2*imode  ), 2._rt*sx_jx[ix]*sz_jx[iz]*wqx*xy.imag());
                        amrex::Gpu::Atomic::AddNoRet( &jy_arr(lo.x+j_jy+ix, lo.y+l_jy+iz, 0, 2*imode-1), 2._rt*sx_jy[ix]*sz_jy[iz]*wqy*xy.real());
                        amrex::Gpu::Atomic::AddNoRet( &jy_arr(lo.x+j_jy+ix, lo.y+l_jy+iz, 0, 2*imode  ), 2._rt*sx_jy[ix]*sz_jy[iz]*wqy*xy.imag());
                        amrex::Gpu::Atomic::AddNoRet( &jz_arr(lo.x+j_jz+ix, lo.y+l_jz+iz, 0, 2*imode-1), 2._rt*sx_jz[ix]*sz_jz[iz]*wqz*xy.real());
                        amrex::Gpu::Atomic::AddNoRet( &jz_arr(lo.x+j_jz+ix, lo.y+l_jz+iz, 0, 2*imode  ), 2._rt*sx_jz[ix]*sz_jz[iz]*wqz*xy.imag());
                        xy = xy*xy0;
                    }
#endif
                }
            }
#elif (defined WARPX_DIM_3D)
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        amrex::Gpu::Atomic::AddNoRet(
                            &jx_arr(lo.x+j_jx+ix, lo.y+k_jx+iy, lo.z+l_jx+iz),
                            sx_jx[ix]*sy_jx[iy]*sz_jx[iz]*wqx);
                        amrex::Gpu::Atomic::AddNoRet(
                            &jy_arr(lo.x+j_jy+ix, lo.y+k_jy+iy, lo.z+l_jy+iz),
                            sx_jy[ix]*sy_jy[iy]*sz_jy[iz]*wqy);
                        amrex::Gpu::Atomic::AddNoRet(
                            &jz_arr(lo.x+j_jz+ix, lo.y+k_jz+iy, lo.z+l_jz+iz),
                            sx_jz[ix]*sy_jz[iy]*sz_jz[iz]*wqz);
                    }
                }
            }
#endif
        }
        );
}

/**
 * \brief Esirkepov Current Deposition for thread thread_num
 *
 * /param GetPosition : A functor for returning the particle position.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param Jx_arr       : Array4 of current density, either full array or tile.
 * \param Jy_arr       : Array4 of current density, either full array or tile.
 * \param Jz_arr       : Array4 of current density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dt           : Time step for particle level
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param q            : species charge.
 * \param n_rz_azimuthal_modes: Number of azimuthal modes when using RZ geometry
 */
template <int depos_order>
void doEsirkepovDepositionShapeN (const GetParticlePosition& GetPosition,
                                  const amrex::ParticleReal * const wp,
                                  const amrex::ParticleReal * const uxp,
                                  const amrex::ParticleReal * const uyp,
                                  const amrex::ParticleReal * const uzp,
                                  const int * ion_lev,
                                  const amrex::Array4<amrex::Real>& Jx_arr,
                                  const amrex::Array4<amrex::Real>& Jy_arr,
                                  const amrex::Array4<amrex::Real>& Jz_arr,
                                  const long np_to_depose,
                                  const amrex::Real dt,
                                  const std::array<amrex::Real,3>& dx,
                                  const std::array<amrex::Real, 3> xyzmin,
                                  const amrex::Dim3 lo,
                                  const amrex::Real q,
                                  const int n_rz_azimuthal_modes)
{
    using namespace amrex;
#if !defined(WARPX_DIM_RZ)
    ignore_unused(n_rz_azimuthal_modes);
#endif

    // Whether ion_lev is a null pointer (do_ionization=0) or a real pointer
    // (do_ionization=1)
    bool const do_ionization = ion_lev;
    Real const dxi = 1.0_rt / dx[0];
#if !(defined WARPX_DIM_RZ)
    Real const dtsdx0 = dt*dxi;
#endif
    Real const xmin = xyzmin[0];
#if (defined WARPX_DIM_3D)
    Real const dyi = 1.0_rt / dx[1];
    Real const dtsdy0 = dt*dyi;
    Real const ymin = xyzmin[1];
#endif
    Real const dzi = 1.0_rt / dx[2];
    Real const dtsdz0 = dt*dzi;
    Real const zmin = xyzmin[2];

#if (defined WARPX_DIM_3D)
    Real const invdtdx = 1.0_rt / (dt*dx[1]*dx[2]);
    Real const invdtdy = 1.0_rt / (dt*dx[0]*dx[2]);
    Real const invdtdz = 1.0_rt / (dt*dx[0]*dx[1]);
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    Real const invdtdx = 1.0_rt / (dt*dx[2]);
    Real const invdtdz = 1.0_rt / (dt*dx[0]);
    Real const invvol = 1.0_rt / (dx[0]*dx[2]);
#endif

#if (defined WARPX_DIM_RZ)
    Complex const I = Complex{0._rt, 1._rt};
#endif

    Real const clightsq = 1.0_rt / ( PhysConst::c * PhysConst::c );

    // Loop over particles and deposit into Jx_arr, Jy_arr and Jz_arr
    amrex::ParallelFor(
        np_to_depose,
        [=] AMREX_GPU_DEVICE (long const ip) {

            // --- Get particle quantities
      Real const gaminv = 1.0_rt/std::sqrt(1.0_rt + uxp[ip]*uxp[ip]*clightsq
                           + uyp[ip]*uyp[ip]*clightsq
                           + uzp[ip]*uzp[ip]*clightsq);

            // wqx, wqy wqz are particle current in each direction
            Real wq = q*wp[ip];
            if (do_ionization){
                wq *= ion_lev[ip];
            }

            ParticleReal xp, yp, zp;
            GetPosition(ip, xp, yp, zp);

            Real const wqx = wq*invdtdx;
#if (defined WARPX_DIM_3D)
            Real const wqy = wq*invdtdy;
#endif
            Real const wqz = wq*invdtdz;

            // computes current and old position in grid units
#if (defined WARPX_DIM_RZ)
            Real const xp_mid = xp - 0.5_rt * dt*uxp[ip]*gaminv;
            Real const yp_mid = yp - 0.5_rt * dt*uyp[ip]*gaminv;
            Real const xp_old = xp - dt*uxp[ip]*gaminv;
            Real const yp_old = yp - dt*uyp[ip]*gaminv;
            Real const rp_new = std::sqrt(xp*xp
                                        + yp*yp);
            Real const rp_mid = std::sqrt(xp_mid*xp_mid + yp_mid*yp_mid);
            Real const rp_old = std::sqrt(xp_old*xp_old + yp_old*yp_old);
            Real costheta_new, sintheta_new;
            if (rp_new > 0._rt) {
                costheta_new = xp/rp_new;
                sintheta_new = yp/rp_new;
            } else {
                costheta_new = 1._rt;
                sintheta_new = 0._rt;
            }
            amrex::Real costheta_mid, sintheta_mid;
            if (rp_mid > 0._rt) {
                costheta_mid = xp_mid/rp_mid;
                sintheta_mid = yp_mid/rp_mid;
            } else {
                costheta_mid = 1._rt;
                sintheta_mid = 0._rt;
            }
            amrex::Real costheta_old, sintheta_old;
            if (rp_old > 0._rt) {
                costheta_old = xp_old/rp_old;
                sintheta_old = yp_old/rp_old;
            } else {
                costheta_old = 1._rt;
                sintheta_old = 0._rt;
            }
            const Complex xy_new0 = Complex{costheta_new, sintheta_new};
            const Complex xy_mid0 = Complex{costheta_mid, sintheta_mid};
            const Complex xy_old0 = Complex{costheta_old, sintheta_old};
            // Keep these double to avoid bug in single precision
            double const x_new = (rp_new - xmin)*dxi;
            double const x_old = (rp_old - xmin)*dxi;
#else
            // Keep these double to avoid bug in single precision
            double const x_new = (xp - xmin)*dxi;
            double const x_old = x_new - dtsdx0*uxp[ip]*gaminv;
#endif
#if (defined WARPX_DIM_3D)
            // Keep these double to avoid bug in single precision
            double const y_new = (yp - ymin)*dyi;
            double const y_old = y_new - dtsdy0*uyp[ip]*gaminv;
#endif
            // Keep these double to avoid bug in single precision
            double const z_new = (zp - zmin)*dzi;
            double const z_old = z_new - dtsdz0*uzp[ip]*gaminv;

#if (defined WARPX_DIM_RZ)
            Real const vy = (-uxp[ip]*sintheta_mid + uyp[ip]*costheta_mid)*gaminv;
#elif (defined WARPX_DIM_XZ)
            Real const vy = uyp[ip]*gaminv;
#endif

            // Shape factor arrays
            // Note that there are extra values above and below
            // to possibly hold the factor for the old particle
            // which can be at a different grid location.
            // Keep these double to avoid bug in single precision
            double sx_new[depos_order + 3] = {0.};
            double sx_old[depos_order + 3] = {0.};
#if (defined WARPX_DIM_3D)
            // Keep these double to avoid bug in single precision
            double sy_new[depos_order + 3] = {0.};
            double sy_old[depos_order + 3] = {0.};
#endif
            // Keep these double to avoid bug in single precision
            double sz_new[depos_order + 3] = {0.};
            double sz_old[depos_order + 3] = {0.};

            // --- Compute shape factors
            // Compute shape factors for position as they are now and at old positions
            // [ijk]_new: leftmost grid point that the particle touches
            Compute_shape_factor< depos_order > compute_shape_factor;
            Compute_shifted_shape_factor< depos_order > compute_shifted_shape_factor;

            const int i_new = compute_shape_factor(sx_new+1, x_new);
            const int i_old = compute_shifted_shape_factor(sx_old, x_old, i_new);
#if (defined WARPX_DIM_3D)
            const int j_new = compute_shape_factor(sy_new+1, y_new);
            const int j_old = compute_shifted_shape_factor(sy_old, y_old, j_new);
#endif
            const int k_new = compute_shape_factor(sz_new+1, z_new);
            const int k_old = compute_shifted_shape_factor(sz_old, z_old, k_new);

            // computes min/max positions of current contributions
            int dil = 1, diu = 1;
            if (i_old < i_new) dil = 0;
            if (i_old > i_new) diu = 0;
#if (defined WARPX_DIM_3D)
            int djl = 1, dju = 1;
            if (j_old < j_new) djl = 0;
            if (j_old > j_new) dju = 0;
#endif
            int dkl = 1, dku = 1;
            if (k_old < k_new) dkl = 0;
            if (k_old > k_new) dku = 0;

#if (defined WARPX_DIM_3D)

            for (int k=dkl; k<=depos_order+2-dku; k++) {
                for (int j=djl; j<=depos_order+2-dju; j++) {
                    amrex::Real sdxi = 0._rt;
                    for (int i=dil; i<=depos_order+1-diu; i++) {
                        sdxi += wqx*(sx_old[i] - sx_new[i])*((sy_new[j] + 0.5_rt*(sy_old[j] - sy_new[j]))*sz_new[k] +
                                                             (0.5_rt*sy_new[j] + 1._rt/3._rt*(sy_old[j] - sy_new[j]))*(sz_old[k] - sz_new[k]));
                        amrex::Gpu::Atomic::AddNoRet( &Jx_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdxi);
                    }
                }
            }
            for (int k=dkl; k<=depos_order+2-dku; k++) {
                for (int i=dil; i<=depos_order+2-diu; i++) {
                    amrex::Real sdyj = 0._rt;
                    for (int j=djl; j<=depos_order+1-dju; j++) {
                        sdyj += wqy*(sy_old[j] - sy_new[j])*((sz_new[k] + 0.5_rt*(sz_old[k] - sz_new[k]))*sx_new[i] +
                                                             (0.5_rt*sz_new[k] + 1._rt/3._rt*(sz_old[k] - sz_new[k]))*(sx_old[i] - sx_new[i]));
                        amrex::Gpu::Atomic::AddNoRet( &Jy_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdyj);
                    }
                }
            }
            for (int j=djl; j<=depos_order+2-dju; j++) {
                for (int i=dil; i<=depos_order+2-diu; i++) {
                    amrex::Real sdzk = 0._rt;
                    for (int k=dkl; k<=depos_order+1-dku; k++) {
                        sdzk += wqz*(sz_old[k] - sz_new[k])*((sx_new[i] + 0.5_rt*(sx_old[i] - sx_new[i]))*sy_new[j] +
                                                             (0.5_rt*sx_new[i] + 1._rt/3._rt*(sx_old[i] - sx_new[i]))*(sy_old[j] - sy_new[j]));
                        amrex::Gpu::Atomic::AddNoRet( &Jz_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdzk);
                    }
                }
            }

#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)

            for (int k=dkl; k<=depos_order+2-dku; k++) {
                amrex::Real sdxi = 0._rt;
                for (int i=dil; i<=depos_order+1-diu; i++) {
                    sdxi += wqx*(sx_old[i] - sx_new[i])*(sz_new[k] + 0.5_rt*(sz_old[k] - sz_new[k]));
                    amrex::Gpu::Atomic::AddNoRet( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 0), sdxi);
#if (defined WARPX_DIM_RZ)
                    Complex xy_mid = xy_mid0; // Throughout the following loop, xy_mid takes the value e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 comes from the normalization of the modes
                        const Complex djr_cmplx = 2._rt *sdxi*xy_mid;
                        amrex::Gpu::Atomic::AddNoRet( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode-1), djr_cmplx.real());
                        amrex::Gpu::Atomic::AddNoRet( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode), djr_cmplx.imag());
                        xy_mid = xy_mid*xy_mid0;
                    }
#endif
                }
            }
            for (int k=dkl; k<=depos_order+2-dku; k++) {
                for (int i=dil; i<=depos_order+2-diu; i++) {
                    Real const sdyj = wq*vy*invvol*((sz_new[k] + 0.5_rt * (sz_old[k] - sz_new[k]))*sx_new[i] +
                                                           (0.5_rt * sz_new[k] + 1._rt / 3._rt *(sz_old[k] - sz_new[k]))*(sx_old[i] - sx_new[i]));
                    amrex::Gpu::Atomic::AddNoRet( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 0), sdyj);
#if (defined WARPX_DIM_RZ)
                    Complex xy_new = xy_new0;
                    Complex xy_mid = xy_mid0;
                    Complex xy_old = xy_old0;
                    // Throughout the following loop, xy_ takes the value e^{i m theta_}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 comes from the normalization of the modes
                        // The minus sign comes from the different convention with respect to Davidson et al.
                        const Complex djt_cmplx = -2._rt * I*(i_new-1 + i + xmin*dxi)*wq*invdtdx/(amrex::Real)imode
                                                  *(Complex(sx_new[i]*sz_new[k], 0.)*(xy_new - xy_mid)
                                                  + Complex(sx_old[i]*sz_old[k], 0.)*(xy_mid - xy_old));
                        amrex::Gpu::Atomic::AddNoRet( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode-1), djt_cmplx.real());
                        amrex::Gpu::Atomic::AddNoRet( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode), djt_cmplx.imag());
                        xy_new = xy_new*xy_new0;
                        xy_mid = xy_mid*xy_mid0;
                        xy_old = xy_old*xy_old0;
                    }
#endif
                }
            }
            for (int i=dil; i<=depos_order+2-diu; i++) {
                Real sdzk = 0._rt;
                for (int k=dkl; k<=depos_order+1-dku; k++) {
                    sdzk += wqz*(sz_old[k] - sz_new[k])*(sx_new[i] + 0.5_rt * (sx_old[i] - sx_new[i]));
                    amrex::Gpu::Atomic::AddNoRet( &Jz_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 0), sdzk);
#if (defined WARPX_DIM_RZ)
                    Complex xy_mid = xy_mid0; // Throughout the following loop, xy_mid takes the value e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 comes from the normalization of the modes
                        const Complex djz_cmplx = 2._rt * sdzk * xy_mid;
                        amrex::Gpu::Atomic::AddNoRet( &Jz_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode-1), djz_cmplx.real());
                        amrex::Gpu::Atomic::AddNoRet( &Jz_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode), djz_cmplx.imag());
                        xy_mid = xy_mid*xy_mid0;
                    }
#endif
                }
            }

#endif
        }
        );
}

/**
 * \brief Vay current deposition
 * (<a href="https://doi.org/10.1016/j.jcp.2013.03.010"> Vay et al, 2013</a>)
 * for thread \c thread_num: deposit \c D in real space and store the result in
 * \c jx_fab, \c jy_fab, \c jz_fab
 *
 * \param[in] GetPosition  Functor that returns the particle position
 * \param[in] wp           Pointer to array of particle weights
 * \param[in] uxp          Pointer to arrays of particle momentum along \c x
 * \param[in] uyp          Pointer to arrays of particle momentum along \c y
 * \param[in] uzp          Pointer to arrays of particle momentum along \c z
 * \param[in] ion_lev      Pointer to array of particle ionization level. This is
                           required to have the charge of each macroparticle since \c q
                           is a scalar. For non-ionizable species, \c ion_lev is \c null
 * \param[in,out] jx_fab   FArrayBox of current density, either full array or tile
 * \param[in,out] jy_fab   FArrayBox of current density, either full array or tile
 * \param[in,out] jz_fab   FArrayBox of current density, either full array or tile
 * \param[in] np_to_depose Number of particles for which current is deposited
 * \param[in] dt           Time step for particle level
 * \param[in] dx           3D cell size
 * \param[in] xyzmin       3D lower bounds of physical domain
 * \param[in] lo           Dimension-agnostic lower bounds of index domain
 * \param[in] q            Species charge
 * \param[in] n_rz_azimuthal_modes Number of azimuthal modes in RZ geometry
 */
template <int depos_order>
void doVayDepositionShapeN (const GetParticlePosition& GetPosition,
                            const amrex::ParticleReal* const wp,
                            const amrex::ParticleReal* const uxp,
                            const amrex::ParticleReal* const uyp,
                            const amrex::ParticleReal* const uzp,
                            const int* const ion_lev,
                            amrex::FArrayBox& jx_fab,
                            amrex::FArrayBox& jy_fab,
                            amrex::FArrayBox& jz_fab,
                            const long np_to_depose,
                            const amrex::Real dt,
                            const std::array<amrex::Real,3>& dx,
                            const std::array<amrex::Real,3>& xyzmin,
                            const amrex::Dim3 lo,
                            const amrex::Real q,
                            const int n_rz_azimuthal_modes)
{
#if (defined WARPX_DIM_RZ)
    amrex::ignore_unused(GetPosition,
        wp, uxp, uyp, uzp, ion_lev, jx_fab, jy_fab, jz_fab,
        np_to_depose, dt, dx, xyzmin, lo, q, n_rz_azimuthal_modes);
    amrex::Abort("Vay deposition not implemented in RZ geometry");
#endif

#if !(defined WARPX_DIM_RZ)
    amrex::ignore_unused(n_rz_azimuthal_modes);

    // If ion_lev is a null pointer, then do_ionization=0, else do_ionization=1
    const bool do_ionization = ion_lev;

    // Inverse cell volume in each direction
    const amrex::Real dxi = 1._rt / dx[0];
    const amrex::Real dzi = 1._rt / dx[2];
#if (defined WARPX_DIM_3D)
    const amrex::Real dyi = 1._rt / dx[1];
#endif

    // Inverse of time step
    const amrex::Real invdt = 1._rt / dt;

    // Total inverse cell volume
#if   (defined WARPX_DIM_XZ)
    const amrex::Real invvol = dxi * dzi;
#elif (defined WARPX_DIM_3D)
    const amrex::Real invvol = dxi * dyi * dzi;
#endif

    // Lower bound of physical domain in each direction
    const amrex::Real xmin = xyzmin[0];
    const amrex::Real zmin = xyzmin[2];
#if (defined WARPX_DIM_3D)
    const amrex::Real ymin = xyzmin[1];
#endif

    // Auxiliary constants
#if (defined WARPX_DIM_3D)
    const amrex::Real onethird = 1._rt / 3._rt;
    const amrex::Real onesixth = 1._rt / 6._rt;
#endif

    // Inverse of light speed squared
    const amrex::Real invcsq = 1._rt / (PhysConst::c * PhysConst::c);

    // Arrays where D will be stored
    amrex::Array4<amrex::Real> const& jx_arr = jx_fab.array();
    amrex::Array4<amrex::Real> const& jy_arr = jy_fab.array();
    amrex::Array4<amrex::Real> const& jz_arr = jz_fab.array();

    // Loop over particles and deposit (Dx,Dy,Dz) into jx_fab, jy_fab and jz_fab
    amrex::ParallelFor(np_to_depose, [=] AMREX_GPU_DEVICE (long ip)
    {
        // Inverse of Lorentz factor gamma
        const amrex::Real invgam = 1._rt / std::sqrt(1._rt + uxp[ip] * uxp[ip] * invcsq
                                                           + uyp[ip] * uyp[ip] * invcsq
                                                           + uzp[ip] * uzp[ip] * invcsq);
        // Product of particle charges and weights
        amrex::Real wq = q * wp[ip];
        if (do_ionization) wq *= ion_lev[ip];

        // Current particle positions (in physical units)
        amrex::ParticleReal xp, yp, zp;
        GetPosition(ip, xp, yp, zp);

        // Particle velocities
        const amrex::Real vx = uxp[ip] * invgam;
        const amrex::Real vy = uyp[ip] * invgam;
        const amrex::Real vz = uzp[ip] * invgam;

        // Particle current densities
#if (defined WARPX_DIM_XZ)
        const amrex::Real wqy = wq * vy * invvol;
#endif

        // Current and old particle positions in grid units
        // Keep these double to avoid bug in single precision.
        double const x_new = (xp - xmin) * dxi;
        double const x_old = x_new - vx * dt * dxi;
#if (defined WARPX_DIM_3D)
        // Keep these double to avoid bug in single precision.
        double const y_new = (yp - ymin) * dyi;
        double const y_old = y_new - vy * dt * dyi;
#endif
        // Keep these double to avoid bug in single precision.
        double const z_new = (zp - zmin) * dzi;
        double const z_old = z_new - vz * dt * dzi;

        // Shape factor arrays for current and old positions (nodal)
        // Keep these double to avoid bug in single precision.
        double sx_new[depos_order+1] = {0.};
        double sx_old[depos_order+1] = {0.};
#if (defined WARPX_DIM_3D)
        // Keep these double to avoid bug in single precision.
        double sy_new[depos_order+1] = {0.};
        double sy_old[depos_order+1] = {0.};
#endif
        // Keep these double to avoid bug in single precision.
        double sz_new[depos_order+1] = {0.};
        double sz_old[depos_order+1] = {0.};

        // Compute shape factors for current positions

        // i_new leftmost grid point in x that the particle touches
        // sx_new shape factor along x for the centering of each current
        Compute_shape_factor< depos_order > const compute_shape_factor;
        const int i_new = compute_shape_factor(sx_new, x_new);
#if (defined WARPX_DIM_3D)
        // j_new leftmost grid point in y that the particle touches
        // sy_new shape factor along y for the centering of each current
        const int j_new = compute_shape_factor(sy_new, y_new);
#endif
        // k_new leftmost grid point in z that the particle touches
        // sz_new shape factor along z for the centering of each current
        const int k_new = compute_shape_factor(sz_new, z_new);

        // Compute shape factors for old positions

        // i_old leftmost grid point in x that the particle touches
        // sx_old shape factor along x for the centering of each current
        const int i_old = compute_shape_factor(sx_old, x_old);
#if (defined WARPX_DIM_3D)
        // j_old leftmost grid point in y that the particle touches
        // sy_old shape factor along y for the centering of each current
        const int j_old = compute_shape_factor(sy_old, y_old);
#endif
        // k_old leftmost grid point in z that the particle touches
        // sz_old shape factor along z for the centering of each current
        const int k_old = compute_shape_factor(sz_old, z_old);

        // Deposit current into jx_arr, jy_arr and jz_arr
#if (defined WARPX_DIM_XZ)

        for (int k=0; k<=depos_order; k++) {
            for (int i=0; i<=depos_order; i++) {

                // Re-casting sx_new and sz_new from double to amrex::Real so that
                // Atomic::Add has consistent types in its argument
                amrex::Real const sxn_szn = sx_new[i] * sz_new[k];
                amrex::Real const sxo_szn = sx_old[i] * sz_new[k];
                amrex::Real const sxn_szo = sx_new[i] * sz_old[k];
                amrex::Real const sxo_szo = sx_old[i] * sz_old[k];

                // Jx
                amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_new + i, lo.y + k_new + k, 0, 0),
                    wq * invvol * invdt * 0.5_rt * sxn_szn);

                amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_old + i, lo.y + k_new + k, 0, 0),
                    - wq * invvol * invdt * 0.5_rt * sxo_szn);

                amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_new + i, lo.y + k_old + k, 0, 0),
                    wq * invvol * invdt * 0.5_rt * sxn_szo);

                amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_old + i, lo.y + k_old + k, 0, 0),
                    - wq * invvol * invdt * 0.5_rt * sxo_szo);

                // Jy
                amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_new + i, lo.y + k_new + k, 0, 0),
                    wqy * 0.25_rt * sxn_szn);

                amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_new + i, lo.y + k_old + k, 0, 0),
                    wqy * 0.25_rt * sxn_szo);

                amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_old + i, lo.y + k_new + k, 0, 0),
                    wqy * 0.25_rt * sxo_szn);

                amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_old + i, lo.y + k_old + k, 0, 0),
                    wqy * 0.25_rt * sxo_szo);

                // Jz
                amrex::Gpu::Atomic::AddNoRet(&jz_arr(lo.x + i_new + i, lo.y + k_new + k, 0, 0),
                    wq * invvol * invdt * 0.5_rt * sxn_szn);

                amrex::Gpu::Atomic::AddNoRet(&jz_arr(lo.x+i_new+i,lo.y+k_old+k,0,0),
                    - wq * invvol * invdt * 0.5_rt * sxn_szo);

                amrex::Gpu::Atomic::AddNoRet(&jz_arr(lo.x+i_old+i,lo.y+k_new+k,0,0),
                    wq * invvol * invdt * 0.5_rt * sxo_szn);

                amrex::Gpu::Atomic::AddNoRet(&jz_arr(lo.x + i_old + i, lo.y + k_old + k, 0, 0),
                    - wq * invvol * invdt * 0.5_rt * sxo_szo);
            }
        }

#elif (defined WARPX_DIM_3D)

        for (int k=0; k<=depos_order; k++) {
            for (int j=0; j<=depos_order; j++) {

                amrex::Real const syn_szn = sy_new[j] * sz_new[k];
                amrex::Real const syo_szn = sy_old[j] * sz_new[k];
                amrex::Real const syn_szo = sy_new[j] * sz_old[k];
                amrex::Real const syo_szo = sy_old[j] * sz_old[k];

                for (int i=0; i<=depos_order; i++) {

                    amrex::Real const sxn_syn_szn = sx_new[i] * syn_szn;
                    amrex::Real const sxo_syn_szn = sx_old[i] * syn_szn;
                    amrex::Real const sxn_syo_szn = sx_new[i] * syo_szn;
                    amrex::Real const sxo_syo_szn = sx_old[i] * syo_szn;
                    amrex::Real const sxn_syn_szo = sx_new[i] * syn_szo;
                    amrex::Real const sxo_syn_szo = sx_old[i] * syn_szo;
                    amrex::Real const sxn_syo_szo = sx_new[i] * syo_szo;
                    amrex::Real const sxo_syo_szo = sx_old[i] * syo_szo;

                    // Jx
                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_new + i, lo.y + j_new + j, lo.z + k_new + k),
                        wq * invvol * invdt * onethird * sxn_syn_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_old + i, lo.y + j_new + j, lo.z + k_new + k),
                        - wq * invvol * invdt * onethird * sxo_syn_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_new + i, lo.y + j_old + j, lo.z + k_new + k),
                        wq * invvol * invdt * onesixth * sxn_syo_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_old + i, lo.y + j_old + j,lo.z + k_new + k),
                        - wq * invvol * invdt * onesixth * sxo_syo_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_new + i, lo.y + j_new + j, lo.z + k_old + k),
                        wq * invvol * invdt * onesixth * sxn_syn_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_old + i, lo.y + j_new + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onesixth * sxo_syn_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_new + i, lo.y + j_old + j, lo.z + k_old + k),
                        wq * invvol * invdt * onethird * sxn_syo_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jx_arr(lo.x + i_old + i, lo.y + j_old + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onethird * sxo_syo_szo);

                    // Jy
                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_new + i, lo.y + j_new + j, lo.z + k_new + k),
                        wq * invvol * invdt * onethird * sxn_syn_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_new + i, lo.y + j_old + j, lo.z + k_new + k),
                        - wq * invvol * invdt * onethird * sxn_syo_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_old + i, lo.y + j_new + j, lo.z + k_new + k),
                        wq * invvol * invdt * onesixth * sxo_syn_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_old + i, lo.y + j_old + j, lo.z + k_new + k),
                        - wq * invvol * invdt * onesixth * sxo_syo_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_new + i, lo.y + j_new + j, lo.z + k_old + k),
                        wq * invvol * invdt * onesixth * sxn_syn_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_new + i, lo.y + j_old + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onesixth * sxn_syo_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_old + i, lo.y + j_new + j, lo.z + k_old + k),
                        wq * invvol * invdt * onethird * sxo_syn_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jy_arr(lo.x + i_old + i, lo.y + j_old + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onethird * sxo_syo_szo);

                    // Jz
                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_new + i, lo.y + j_new + j, lo.z + k_new + k),
                        wq * invvol * invdt * onethird * sxn_syn_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_new + i, lo.y + j_new + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onethird * sxn_syn_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_old + i, lo.y + j_new + j, lo.z + k_new + k),
                        wq * invvol * invdt * onesixth * sxo_syn_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_old + i, lo.y + j_new + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onesixth * sxo_syn_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_new + i, lo.y + j_old + j, lo.z + k_new + k),
                        wq * invvol * invdt * onesixth * sxn_syo_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_new + i, lo.y + j_old + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onesixth * sxn_syo_szo);

                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_old + i, lo.y + j_old + j, lo.z + k_new + k),
                        wq * invvol * invdt * onethird * sxo_syo_szn);

                    amrex::Gpu::Atomic::AddNoRet(&jz_arr( lo.x + i_old + i, lo.y + j_old + j, lo.z + k_old + k),
                        - wq * invvol * invdt * onethird * sxo_syo_szo);
                }
            }
        }
#endif
    } );
#endif // #if !(defined WARPX_DIM_RZ)
}
#endif // CURRENTDEPOSITION_H_
