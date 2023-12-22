/* Copyright 2019 Axel Huebl, David Grote, Maxence Thevenet
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef FIELDGATHER_H_
#define FIELDGATHER_H_

#include "Particles/Gather/GetExternalFields.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"
#include "Utils/WarpX_Complex.H"

#include <AMReX.H>

/**
 * \brief Field gather for a single particle
 *
 * \tparam depos_order              Particle shape order
 * \tparam galerkin_interpolation   Lower the order of the particle shape by
 *                                  this value (0/1) for the parallel field component
 * \param xp,yp,zp                        Particle position coordinates
 * \param Exp,Eyp,Ezp                     Electric field on particles.
 * \param Bxp,Byp,Bzp                     Magnetic field on particles.
 * \param ex_arr,ey_arr,ez_arr            Array4 of the electric field, either full array or tile.
 * \param bx_arr,by_arr,bz_arr            Array4 of the magnetic field, either full array or tile.
 * \param ex_type,ey_type,ez_type         IndexType of the electric field
 * \param bx_type,by_type,bz_type         IndexType of the magnetic field
 * \param dx                        3D cell spacing
 * \param xyzmin                    Physical lower bounds of domain in x, y, z.
 * \param lo                        Index lower bounds of domain.
 * \param n_rz_azimuthal_modes       Number of azimuthal modes when using RZ geometry
 */
template <int depos_order, int galerkin_interpolation>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void doGatherShapeN (const amrex::ParticleReal xp,
                     const amrex::ParticleReal yp,
                     const amrex::ParticleReal zp,
                     amrex::ParticleReal& Exp,
                     amrex::ParticleReal& Eyp,
                     amrex::ParticleReal& Ezp,
                     amrex::ParticleReal& Bxp,
                     amrex::ParticleReal& Byp,
                     amrex::ParticleReal& Bzp,
                     amrex::Array4<amrex::Real const> const& ex_arr,
                     amrex::Array4<amrex::Real const> const& ey_arr,
                     amrex::Array4<amrex::Real const> const& ez_arr,
                     amrex::Array4<amrex::Real const> const& bx_arr,
                     amrex::Array4<amrex::Real const> const& by_arr,
                     amrex::Array4<amrex::Real const> const& bz_arr,
                     const amrex::IndexType ex_type,
                     const amrex::IndexType ey_type,
                     const amrex::IndexType ez_type,
                     const amrex::IndexType bx_type,
                     const amrex::IndexType by_type,
                     const amrex::IndexType bz_type,
                     const amrex::GpuArray<amrex::Real, 3>& dx,
                     const amrex::GpuArray<amrex::Real, 3>& xyzmin,
                     const amrex::Dim3& lo,
                     const int n_rz_azimuthal_modes)
{
    using namespace amrex;

#if defined(WARPX_DIM_XZ)
    amrex::ignore_unused(yp);
#endif

#if defined(WARPX_DIM_1D_Z)
    amrex::ignore_unused(xp,yp);
#endif

#ifndef WARPX_DIM_RZ
    amrex::ignore_unused(n_rz_azimuthal_modes);
#endif

#if (AMREX_SPACEDIM >= 2)
    const amrex::Real dxi = 1.0_rt/dx[0];
#endif
    const amrex::Real dzi = 1.0_rt/dx[2];
#if defined(WARPX_DIM_3D)
    const amrex::Real dyi = 1.0_rt/dx[1];
#endif

#if (AMREX_SPACEDIM >= 2)
    const amrex::Real xmin = xyzmin[0];
#endif
#if defined(WARPX_DIM_3D)
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real zmin = xyzmin[2];

    constexpr int zdir = WARPX_ZINDEX;
    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // --- Compute shape factors

    Compute_shape_factor< depos_order > const compute_shape_factor;
    Compute_shape_factor<depos_order - galerkin_interpolation > const compute_shape_factor_galerkin;

#if (AMREX_SPACEDIM >= 2)
    // x direction
    // Get particle position
#ifdef WARPX_DIM_RZ
    const amrex::Real rp = std::sqrt(xp*xp + yp*yp);
    const amrex::Real x = (rp - xmin)*dxi;
#else
    const amrex::Real x = (xp-xmin)*dxi;
#endif

    // j_[eb][xyz] leftmost grid point in x that the particle touches for the centering of each current
    // sx_[eb][xyz] shape factor along x for the centering of each current
    // There are only two possible centerings, node or cell centered, so at most only two shape factor
    // arrays will be needed.
    amrex::Real sx_node[depos_order + 1];
    amrex::Real sx_cell[depos_order + 1];
    amrex::Real sx_node_galerkin[depos_order + 1 - galerkin_interpolation] = {0._rt};
    amrex::Real sx_cell_galerkin[depos_order + 1 - galerkin_interpolation] = {0._rt};

    int j_node = 0;
    int j_cell = 0;
    int j_node_v = 0;
    int j_cell_v = 0;
    if ((ey_type[0] == NODE) || (ez_type[0] == NODE) || (bx_type[0] == NODE)) {
        j_node = compute_shape_factor(sx_node, x);
    }
    if ((ey_type[0] == CELL) || (ez_type[0] == CELL) || (bx_type[0] == CELL)) {
        j_cell = compute_shape_factor(sx_cell, x - 0.5_rt);
    }
    if ((ex_type[0] == NODE) || (by_type[0] == NODE) || (bz_type[0] == NODE)) {
        j_node_v = compute_shape_factor_galerkin(sx_node_galerkin, x);
    }
    if ((ex_type[0] == CELL) || (by_type[0] == CELL) || (bz_type[0] == CELL)) {
        j_cell_v = compute_shape_factor_galerkin(sx_cell_galerkin, x - 0.5_rt);
    }
    const amrex::Real (&sx_ex)[depos_order + 1 - galerkin_interpolation] = ((ex_type[0] == NODE) ? sx_node_galerkin : sx_cell_galerkin);
    const amrex::Real (&sx_ey)[depos_order + 1             ] = ((ey_type[0] == NODE) ? sx_node   : sx_cell  );
    const amrex::Real (&sx_ez)[depos_order + 1             ] = ((ez_type[0] == NODE) ? sx_node   : sx_cell  );
    const amrex::Real (&sx_bx)[depos_order + 1             ] = ((bx_type[0] == NODE) ? sx_node   : sx_cell  );
    const amrex::Real (&sx_by)[depos_order + 1 - galerkin_interpolation] = ((by_type[0] == NODE) ? sx_node_galerkin : sx_cell_galerkin);
    const amrex::Real (&sx_bz)[depos_order + 1 - galerkin_interpolation] = ((bz_type[0] == NODE) ? sx_node_galerkin : sx_cell_galerkin);
    int const j_ex = ((ex_type[0] == NODE) ? j_node_v : j_cell_v);
    int const j_ey = ((ey_type[0] == NODE) ? j_node   : j_cell  );
    int const j_ez = ((ez_type[0] == NODE) ? j_node   : j_cell  );
    int const j_bx = ((bx_type[0] == NODE) ? j_node   : j_cell  );
    int const j_by = ((by_type[0] == NODE) ? j_node_v : j_cell_v);
    int const j_bz = ((bz_type[0] == NODE) ? j_node_v : j_cell_v);
#endif

#if defined(WARPX_DIM_3D)
    // y direction
    const amrex::Real y = (yp-ymin)*dyi;
    amrex::Real sy_node[depos_order + 1];
    amrex::Real sy_cell[depos_order + 1];
    amrex::Real sy_node_v[depos_order + 1 - galerkin_interpolation];
    amrex::Real sy_cell_v[depos_order + 1 - galerkin_interpolation];
    int k_node = 0;
    int k_cell = 0;
    int k_node_v = 0;
    int k_cell_v = 0;
    if ((ex_type[1] == NODE) || (ez_type[1] == NODE) || (by_type[1] == NODE)) {
        k_node = compute_shape_factor(sy_node, y);
    }
    if ((ex_type[1] == CELL) || (ez_type[1] == CELL) || (by_type[1] == CELL)) {
        k_cell = compute_shape_factor(sy_cell, y - 0.5_rt);
    }
    if ((ey_type[1] == NODE) || (bx_type[1] == NODE) || (bz_type[1] == NODE)) {
        k_node_v = compute_shape_factor_galerkin(sy_node_v, y);
    }
    if ((ey_type[1] == CELL) || (bx_type[1] == CELL) || (bz_type[1] == CELL)) {
        k_cell_v = compute_shape_factor_galerkin(sy_cell_v, y - 0.5_rt);
    }
    const amrex::Real (&sy_ex)[depos_order + 1             ] = ((ex_type[1] == NODE) ? sy_node   : sy_cell  );
    const amrex::Real (&sy_ey)[depos_order + 1 - galerkin_interpolation] = ((ey_type[1] == NODE) ? sy_node_v : sy_cell_v);
    const amrex::Real (&sy_ez)[depos_order + 1             ] = ((ez_type[1] == NODE) ? sy_node   : sy_cell  );
    const amrex::Real (&sy_bx)[depos_order + 1 - galerkin_interpolation] = ((bx_type[1] == NODE) ? sy_node_v : sy_cell_v);
    const amrex::Real (&sy_by)[depos_order + 1             ] = ((by_type[1] == NODE) ? sy_node   : sy_cell  );
    const amrex::Real (&sy_bz)[depos_order + 1 - galerkin_interpolation] = ((bz_type[1] == NODE) ? sy_node_v : sy_cell_v);
    int const k_ex = ((ex_type[1] == NODE) ? k_node   : k_cell  );
    int const k_ey = ((ey_type[1] == NODE) ? k_node_v : k_cell_v);
    int const k_ez = ((ez_type[1] == NODE) ? k_node   : k_cell  );
    int const k_bx = ((bx_type[1] == NODE) ? k_node_v : k_cell_v);
    int const k_by = ((by_type[1] == NODE) ? k_node   : k_cell  );
    int const k_bz = ((bz_type[1] == NODE) ? k_node_v : k_cell_v);

#endif
    // z direction
    const amrex::Real z = (zp-zmin)*dzi;
    amrex::Real sz_node[depos_order + 1];
    amrex::Real sz_cell[depos_order + 1];
    amrex::Real sz_node_v[depos_order + 1 - galerkin_interpolation];
    amrex::Real sz_cell_v[depos_order + 1 - galerkin_interpolation];
    int l_node = 0;
    int l_cell = 0;
    int l_node_v = 0;
    int l_cell_v = 0;
    if ((ex_type[zdir] == NODE) || (ey_type[zdir] == NODE) || (bz_type[zdir] == NODE)) {
        l_node = compute_shape_factor(sz_node, z);
    }
    if ((ex_type[zdir] == CELL) || (ey_type[zdir] == CELL) || (bz_type[zdir] == CELL)) {
        l_cell = compute_shape_factor(sz_cell, z - 0.5_rt);
    }
    if ((ez_type[zdir] == NODE) || (bx_type[zdir] == NODE) || (by_type[zdir] == NODE)) {
        l_node_v = compute_shape_factor_galerkin(sz_node_v, z);
    }
    if ((ez_type[zdir] == CELL) || (bx_type[zdir] == CELL) || (by_type[zdir] == CELL)) {
        l_cell_v = compute_shape_factor_galerkin(sz_cell_v, z - 0.5_rt);
    }
    const amrex::Real (&sz_ex)[depos_order + 1             ] = ((ex_type[zdir] == NODE) ? sz_node   : sz_cell  );
    const amrex::Real (&sz_ey)[depos_order + 1             ] = ((ey_type[zdir] == NODE) ? sz_node   : sz_cell  );
    const amrex::Real (&sz_ez)[depos_order + 1 - galerkin_interpolation] = ((ez_type[zdir] == NODE) ? sz_node_v : sz_cell_v);
    const amrex::Real (&sz_bx)[depos_order + 1 - galerkin_interpolation] = ((bx_type[zdir] == NODE) ? sz_node_v : sz_cell_v);
    const amrex::Real (&sz_by)[depos_order + 1 - galerkin_interpolation] = ((by_type[zdir] == NODE) ? sz_node_v : sz_cell_v);
    const amrex::Real (&sz_bz)[depos_order + 1             ] = ((bz_type[zdir] == NODE) ? sz_node   : sz_cell  );
    int const l_ex = ((ex_type[zdir] == NODE) ? l_node   : l_cell  );
    int const l_ey = ((ey_type[zdir] == NODE) ? l_node   : l_cell  );
    int const l_ez = ((ez_type[zdir] == NODE) ? l_node_v : l_cell_v);
    int const l_bx = ((bx_type[zdir] == NODE) ? l_node_v : l_cell_v);
    int const l_by = ((by_type[zdir] == NODE) ? l_node_v : l_cell_v);
    int const l_bz = ((bz_type[zdir] == NODE) ? l_node   : l_cell  );


    // Each field is gathered in a separate block of
    // AMREX_SPACEDIM nested loops because the deposition
    // order can differ for each component of each field
    // when galerkin_interpolation is set to 1

#if defined(WARPX_DIM_1D_Z)
    // Gather field on particle Eyp from field on grid ey_arr
    // Gather field on particle Exp from field on grid ex_arr
    // Gather field on particle Bzp from field on grid bz_arr
    for (int iz=0; iz<=depos_order; iz++){
        Eyp += sz_ey[iz]*
            ey_arr(lo.x+l_ey+iz, 0, 0, 0);
        Exp += sz_ex[iz]*
            ex_arr(lo.x+l_ex+iz, 0, 0, 0);
        Bzp += sz_bz[iz]*
            bz_arr(lo.x+l_bz+iz, 0, 0, 0);
    }

    // Gather field on particle Byp from field on grid by_arr
    // Gather field on particle Ezp from field on grid ez_arr
    // Gather field on particle Bxp from field on grid bx_arr
    for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
        Ezp += sz_ez[iz]*
            ez_arr(lo.x+l_ez+iz, 0, 0, 0);
        Bxp += sz_bx[iz]*
            bx_arr(lo.x+l_bx+iz, 0, 0, 0);
        Byp += sz_by[iz]*
            by_arr(lo.x+l_by+iz, 0, 0, 0);
    }

#elif defined(WARPX_DIM_XZ)
    // Gather field on particle Eyp from field on grid ey_arr
    for (int iz=0; iz<=depos_order; iz++){
        for (int ix=0; ix<=depos_order; ix++){
            Eyp += sx_ey[ix]*sz_ey[iz]*
                ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 0);
        }
    }
    // Gather field on particle Exp from field on grid ex_arr
    // Gather field on particle Bzp from field on grid bz_arr
    for (int iz=0; iz<=depos_order; iz++){
        for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
            Exp += sx_ex[ix]*sz_ex[iz]*
                ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 0);
            Bzp += sx_bz[ix]*sz_bz[iz]*
                bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 0);
        }
    }
    // Gather field on particle Ezp from field on grid ez_arr
    // Gather field on particle Bxp from field on grid bx_arr
    for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
        for (int ix=0; ix<=depos_order; ix++){
            Ezp += sx_ez[ix]*sz_ez[iz]*
                ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 0);
            Bxp += sx_bx[ix]*sz_bx[iz]*
                bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 0);
        }
    }
    // Gather field on particle Byp from field on grid by_arr
    for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
        for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
            Byp += sx_by[ix]*sz_by[iz]*
                by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 0);
        }
    }

#elif defined(WARPX_DIM_RZ)

    amrex::ParticleReal Erp = 0.;
    amrex::ParticleReal Ethetap = 0.;
    amrex::ParticleReal Brp = 0.;
    amrex::ParticleReal Bthetap = 0.;

    // Gather field on particle Ethetap from field on grid ey_arr
    for (int iz=0; iz<=depos_order; iz++){
        for (int ix=0; ix<=depos_order; ix++){
            Ethetap += sx_ey[ix]*sz_ey[iz]*
                ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 0);
        }
    }
    // Gather field on particle Erp from field on grid ex_arr
    // Gather field on particle Bzp from field on grid bz_arr
    for (int iz=0; iz<=depos_order; iz++){
        for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
            Erp += sx_ex[ix]*sz_ex[iz]*
                ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 0);
            Bzp += sx_bz[ix]*sz_bz[iz]*
                bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 0);
        }
    }
    // Gather field on particle Ezp from field on grid ez_arr
    // Gather field on particle Brp from field on grid bx_arr
    for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
        for (int ix=0; ix<=depos_order; ix++){
            Ezp += sx_ez[ix]*sz_ez[iz]*
                ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 0);
            Brp += sx_bx[ix]*sz_bx[iz]*
                bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 0);
        }
    }
    // Gather field on particle Bthetap from field on grid by_arr
    for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
        for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
            Bthetap += sx_by[ix]*sz_by[iz]*
                by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 0);
        }
    }

    amrex::Real costheta;
    amrex::Real sintheta;
    if (rp > 0.) {
        costheta = xp/rp;
        sintheta = yp/rp;
    } else {
        costheta = 1.;
        sintheta = 0.;
    }
    const Complex xy0 = Complex{costheta, -sintheta};
    Complex xy = xy0;

    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {

        // Gather field on particle Ethetap from field on grid ey_arr
        for (int iz=0; iz<=depos_order; iz++){
            for (int ix=0; ix<=depos_order; ix++){
                const amrex::Real dEy = (+ ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 2*imode-1)*xy.real()
                                         - ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 2*imode)*xy.imag());
                Ethetap += sx_ey[ix]*sz_ey[iz]*dEy;
            }
        }
        // Gather field on particle Erp from field on grid ex_arr
        // Gather field on particle Bzp from field on grid bz_arr
        for (int iz=0; iz<=depos_order; iz++){
            for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
                const amrex::Real dEx = (+ ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 2*imode-1)*xy.real()
                                         - ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 2*imode)*xy.imag());
                Erp += sx_ex[ix]*sz_ex[iz]*dEx;
                const amrex::Real dBz = (+ bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 2*imode-1)*xy.real()
                                         - bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 2*imode)*xy.imag());
                Bzp += sx_bz[ix]*sz_bz[iz]*dBz;
            }
        }
        // Gather field on particle Ezp from field on grid ez_arr
        // Gather field on particle Brp from field on grid bx_arr
        for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
            for (int ix=0; ix<=depos_order; ix++){
                const amrex::Real dEz = (+ ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 2*imode-1)*xy.real()
                                         - ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 2*imode)*xy.imag());
                Ezp += sx_ez[ix]*sz_ez[iz]*dEz;
                const amrex::Real dBx = (+ bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 2*imode-1)*xy.real()
                                         - bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 2*imode)*xy.imag());
                Brp += sx_bx[ix]*sz_bx[iz]*dBx;
            }
        }
        // Gather field on particle Bthetap from field on grid by_arr
        for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
            for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
                const amrex::Real dBy = (+ by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 2*imode-1)*xy.real()
                                         - by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 2*imode)*xy.imag());
                Bthetap += sx_by[ix]*sz_by[iz]*dBy;
            }
        }
        xy = xy*xy0;
    }

    // Convert Erp and Ethetap to Ex and Ey
    Exp += costheta*Erp - sintheta*Ethetap;
    Eyp += costheta*Ethetap + sintheta*Erp;
    Bxp += costheta*Brp - sintheta*Bthetap;
    Byp += costheta*Bthetap + sintheta*Brp;

#else // defined(WARPX_DIM_3D)
    // Gather field on particle Exp from field on grid ex_arr
    for (int iz=0; iz<=depos_order; iz++){
        for (int iy=0; iy<=depos_order; iy++){
            for (int ix=0; ix<= depos_order - galerkin_interpolation; ix++){
                Exp += sx_ex[ix]*sy_ex[iy]*sz_ex[iz]*
                    ex_arr(lo.x+j_ex+ix, lo.y+k_ex+iy, lo.z+l_ex+iz);
            }
        }
    }
    // Gather field on particle Eyp from field on grid ey_arr
    for (int iz=0; iz<=depos_order; iz++){
        for (int iy=0; iy<= depos_order - galerkin_interpolation; iy++){
            for (int ix=0; ix<=depos_order; ix++){
                Eyp += sx_ey[ix]*sy_ey[iy]*sz_ey[iz]*
                    ey_arr(lo.x+j_ey+ix, lo.y+k_ey+iy, lo.z+l_ey+iz);
            }
        }
    }
    // Gather field on particle Ezp from field on grid ez_arr
    for (int iz=0; iz<= depos_order - galerkin_interpolation; iz++){
        for (int iy=0; iy<=depos_order; iy++){
            for (int ix=0; ix<=depos_order; ix++){
                Ezp += sx_ez[ix]*sy_ez[iy]*sz_ez[iz]*
                    ez_arr(lo.x+j_ez+ix, lo.y+k_ez+iy, lo.z+l_ez+iz);
            }
        }
    }
    // Gather field on particle Bzp from field on grid bz_arr
    for (int iz=0; iz<=depos_order; iz++){
        for (int iy=0; iy<= depos_order - galerkin_interpolation; iy++){
            for (int ix=0; ix<= depos_order - galerkin_interpolation; ix++){
                Bzp += sx_bz[ix]*sy_bz[iy]*sz_bz[iz]*
                    bz_arr(lo.x+j_bz+ix, lo.y+k_bz+iy, lo.z+l_bz+iz);
            }
        }
    }
    // Gather field on particle Byp from field on grid by_arr
    for (int iz=0; iz<= depos_order - galerkin_interpolation; iz++){
        for (int iy=0; iy<=depos_order; iy++){
            for (int ix=0; ix<= depos_order - galerkin_interpolation; ix++){
                Byp += sx_by[ix]*sy_by[iy]*sz_by[iz]*
                    by_arr(lo.x+j_by+ix, lo.y+k_by+iy, lo.z+l_by+iz);
            }
        }
    }
    // Gather field on particle Bxp from field on grid bx_arr
    for (int iz=0; iz<= depos_order - galerkin_interpolation; iz++){
        for (int iy=0; iy<= depos_order - galerkin_interpolation; iy++){
            for (int ix=0; ix<=depos_order; ix++){
                Bxp += sx_bx[ix]*sy_bx[iy]*sz_bx[iz]*
                    bx_arr(lo.x+j_bx+ix, lo.y+k_bx+iy, lo.z+l_bx+iz);
            }
        }
    }
#endif
}

/**
 * \brief Energy conserving field gather for thread thread_num for the implicit scheme
 *        This uses the same stencil for the gather that is used for Esirkepov current deposition.
 *
 * \tparam depos_order              Particle shape order
 * \param xp_n,yp_n,zp_n                  Particle position coordinates at start of step
 * \param xp_nph,yp_nph,zp_nph            Particle position coordinates at half step
 * \param Exp,Eyp,Ezp                     Electric field on particles.
 * \param Bxp,Byp,Bzp                     Magnetic field on particles.
 * \param Ex_arr,Ey_arr,Ez_arr            Array4 of the electric field, either full array or tile.
 * \param Bx_arr,By_arr,Bz_arr            Array4 of the magnetic field, either full array or tile.
 * \param Ex_type,Ey_type,Ez_type         IndexType of the electric field
 * \param Bx_type,By_type,Bz_type         IndexType of the magnetic field
 * \param dx                        3D cell spacing
 * \param xyzmin                    Physical lower bounds of domain in x, y, z.
 * \param lo                        Index lower bounds of domain.
 * \param n_rz_azimuthal_modes       Number of azimuthal modes when using RZ geometry
 */
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void doGatherShapeNEsirkepovStencilImplicit (
                     [[maybe_unused]] const amrex::ParticleReal xp_n,
                     [[maybe_unused]] const amrex::ParticleReal yp_n,
                     const amrex::ParticleReal zp_n,
                     [[maybe_unused]] const amrex::ParticleReal xp_nph,
                     [[maybe_unused]] const amrex::ParticleReal yp_nph,
                     const amrex::ParticleReal zp_nph,
                     amrex::ParticleReal& Exp,
                     amrex::ParticleReal& Eyp,
                     amrex::ParticleReal& Ezp,
                     amrex::ParticleReal& Bxp,
                     amrex::ParticleReal& Byp,
                     amrex::ParticleReal& Bzp,
                     amrex::Array4<amrex::Real const> const& Ex_arr,
                     amrex::Array4<amrex::Real const> const& Ey_arr,
                     amrex::Array4<amrex::Real const> const& Ez_arr,
                     amrex::Array4<amrex::Real const> const& Bx_arr,
                     amrex::Array4<amrex::Real const> const& By_arr,
                     amrex::Array4<amrex::Real const> const& Bz_arr,
                     [[maybe_unused]] const amrex::IndexType Ex_type,
                     [[maybe_unused]] const amrex::IndexType Ey_type,
                     [[maybe_unused]] const amrex::IndexType Ez_type,
                     [[maybe_unused]] const amrex::IndexType Bx_type,
                     [[maybe_unused]] const amrex::IndexType By_type,
                     [[maybe_unused]] const amrex::IndexType Bz_type,
                     const amrex::GpuArray<amrex::Real, 3>& dx,
                     const amrex::GpuArray<amrex::Real, 3>& xyzmin,
                     const amrex::Dim3& lo,
                     const int n_rz_azimuthal_modes)
{
    using namespace amrex;
#if !defined(WARPX_DIM_RZ)
    ignore_unused(n_rz_azimuthal_modes);
#endif

#if !defined(WARPX_DIM_1D_Z)
    Real const dxi = 1.0_rt / dx[0];
#endif
#if !defined(WARPX_DIM_1D_Z)
    Real const xmin = xyzmin[0];
#endif
#if defined(WARPX_DIM_3D)
    Real const dyi = 1.0_rt / dx[1];
    Real const ymin = xyzmin[1];
#endif
    Real const dzi = 1.0_rt / dx[2];
    Real const zmin = xyzmin[2];

#if !defined(WARPX_DIM_1D_Z)
    Real constexpr one_third = 1.0_rt / 3.0_rt;
    Real constexpr one_sixth = 1.0_rt / 6.0_rt;
#endif

#if !defined(WARPX_DIM_1D_Z)
    ParticleReal xp_np1 = 2._prt*xp_nph - xp_n;
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
    ParticleReal yp_np1 = 2._prt*yp_nph - yp_n;
#endif
    ParticleReal zp_np1 = 2._prt*zp_nph - zp_n;

    // computes current and old position in grid units
#if defined(WARPX_DIM_RZ)
    amrex::Real const xp_new = xp_np1;
    amrex::Real const yp_new = yp_np1;
    amrex::Real const xp_mid = xp_nph;
    amrex::Real const yp_mid = yp_nph;
    amrex::Real const xp_old = xp_n;
    amrex::Real const yp_old = yp_n;
    amrex::Real const rp_new = std::sqrt(xp_new*xp_new + yp_new*yp_new);
    amrex::Real const rp_old = std::sqrt(xp_old*xp_old + yp_old*yp_old);
    amrex::Real const rp_mid = (rp_new + rp_old)/2._rt;
    amrex::Real costheta_mid, sintheta_mid;
    if (rp_mid > 0._rt) {
        costheta_mid = xp_mid/rp_mid;
        sintheta_mid = yp_mid/rp_mid;
    } else {
        costheta_mid = 1._rt;
        sintheta_mid = 0._rt;
    }
    const Complex xy_mid0 = Complex{costheta_mid, sintheta_mid};
    // Keep these double to avoid bug in single precision
    double const x_new = (rp_new - xmin)*dxi;
    double const x_old = (rp_old - xmin)*dxi;
#else
#if !defined(WARPX_DIM_1D_Z)
    // Keep these double to avoid bug in single precision
    double const x_new = (xp_np1 - xmin)*dxi;
    double const x_old = (xp_n - xmin)*dxi;
#endif
#endif
#if defined(WARPX_DIM_3D)
    // Keep these double to avoid bug in single precision
    double const y_new = (yp_np1 - ymin)*dyi;
    double const y_old = (yp_n - ymin)*dyi;
#endif
    // Keep these double to avoid bug in single precision
    double const z_new = (zp_np1 - zmin)*dzi;
    double const z_old = (zp_n - zmin)*dzi;

    // Shape factor arrays
    // Note that there are extra values above and below
    // to possibly hold the factor for the old particle
    // which can be at a different grid location.
    // Keep these double to avoid bug in single precision
#if !defined(WARPX_DIM_1D_Z)
    double sx_E_new[depos_order + 3] = {0.};
    double sx_E_old[depos_order + 3] = {0.};
#endif
#if defined(WARPX_DIM_3D)
    // Keep these double to avoid bug in single precision
    double sy_E_new[depos_order + 3] = {0.};
    double sy_E_old[depos_order + 3] = {0.};
#endif
    // Keep these double to avoid bug in single precision
    double sz_E_new[depos_order + 3] = {0.};
    double sz_E_old[depos_order + 3] = {0.};

#if defined(WARPX_DIM_3D)
    double sx_B_new[depos_order + 3] = {0.};
    double sx_B_old[depos_order + 3] = {0.};
    double sy_B_new[depos_order + 3] = {0.};
    double sy_B_old[depos_order + 3] = {0.};
    double sz_B_new[depos_order + 3] = {0.};
    double sz_B_old[depos_order + 3] = {0.};
#endif

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    // Special shape functions are needed for By which is cell
    // centered in both x and z. One lower order shape function is used.
    double sx_By_new[depos_order + 2] = {0.};
    double sx_By_old[depos_order + 2] = {0.};
    double sz_By_new[depos_order + 2] = {0.};
    double sz_By_old[depos_order + 2] = {0.};
#endif

    // --- Compute shape factors
    // Compute shape factors for position as they are now and at old positions
    // [ijk]_new: leftmost grid point that the particle touches
    Compute_shape_factor< depos_order > compute_shape_factor;
    Compute_shifted_shape_factor< depos_order > compute_shifted_shape_factor;

#if !defined(WARPX_DIM_1D_Z)
    const int i_E_new = compute_shape_factor(sx_E_new+1, x_new);
    const int i_E_old = compute_shifted_shape_factor(sx_E_old, x_old, i_E_new);
#endif
#if defined(WARPX_DIM_3D)
    const int j_E_new = compute_shape_factor(sy_E_new+1, y_new);
    const int j_E_old = compute_shifted_shape_factor(sy_E_old, y_old, j_E_new);
#endif
    const int k_E_new = compute_shape_factor(sz_E_new+1, z_new);
    const int k_E_old = compute_shifted_shape_factor(sz_E_old, z_old, k_E_new);

#if defined(WARPX_DIM_3D)
    const int i_B_new = compute_shape_factor(sx_B_new+1, x_new - 0.5_rt);
    const int i_B_old = compute_shifted_shape_factor(sx_B_old, x_old - 0.5_rt, i_B_new);
    const int j_B_new = compute_shape_factor(sy_B_new+1, y_new - 0.5_rt);
    const int j_B_old = compute_shifted_shape_factor(sy_B_old, y_old - 0.5_rt, j_B_new);
    const int k_B_new = compute_shape_factor(sz_B_new+1, z_new - 0.5_rt);
    const int k_B_old = compute_shifted_shape_factor(sz_B_old, z_old - 0.5_rt, k_B_new);
#endif

    // computes min/max positions of current contributions
#if !defined(WARPX_DIM_1D_Z)
    int dil_E = 1, diu_E = 1;
    if (i_E_old < i_E_new) { dil_E = 0; }
    if (i_E_old > i_E_new) { diu_E = 0; }
#endif
#if defined(WARPX_DIM_3D)
    int djl_E = 1, dju_E = 1;
    if (j_E_old < j_E_new) { djl_E = 0; }
    if (j_E_old > j_E_new) { dju_E = 0; }
#endif
    int dkl_E = 1, dku_E = 1;
    if (k_E_old < k_E_new) { dkl_E = 0; }
    if (k_E_old > k_E_new) { dku_E = 0; }

#if defined(WARPX_DIM_3D)
    int dil_B = 1, diu_B = 1;
    if (i_B_old < i_B_new) { dil_B = 0; }
    if (i_B_old > i_B_new) { diu_B = 0; }
    int djl_B = 1, dju_B = 1;
    if (j_B_old < j_B_new) { djl_B = 0; }
    if (j_B_old > j_B_new) { dju_B = 0; }
    int dkl_B = 1, dku_B = 1;
    if (k_B_old < k_B_new) { dkl_B = 0; }
    if (k_B_old > k_B_new) { dku_B = 0; }
#endif

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    Compute_shape_factor< depos_order-1 > compute_shape_factor_By;
    Compute_shifted_shape_factor< depos_order-1 > compute_shifted_shape_factor_By;
    const int i_By_new = compute_shape_factor_By(sx_By_new+1, x_new - 0.5_rt);
    const int i_By_old = compute_shifted_shape_factor_By(sx_By_old, x_old - 0.5_rt, i_By_new);
    const int k_By_new = compute_shape_factor_By(sz_By_new+1, z_new - 0.5_rt);
    const int k_By_old = compute_shifted_shape_factor_By(sz_By_old, z_old - 0.5_rt, k_By_new);
    int dil_By = 1, diu_By = 1;
    if (i_By_old < i_By_new) { dil_By = 0; }
    if (i_By_old > i_By_new) { diu_By = 0; }
    int dkl_By = 1, dku_By = 1;
    if (k_By_old < k_By_new) { dkl_By = 0; }
    if (k_By_old > k_By_new) { dku_By = 0; }
#endif

#if defined(WARPX_DIM_3D)

    for (int k=dkl_E; k<=depos_order+2-dku_E; k++) {
        for (int j=djl_E; j<=depos_order+2-dju_E; j++) {
            amrex::Real sdzjk = one_third*(sy_E_new[j]*sz_E_new[k] + sy_E_old[j]*sz_E_old[k])
                               +one_sixth*(sy_E_new[j]*sz_E_old[k] + sy_E_old[j]*sz_E_new[k]);
            amrex::Real sdxi = 0._rt;
            for (int i=dil_E; i<=depos_order+1-diu_E; i++) {
                sdxi += (sx_E_old[i] - sx_E_new[i]);
                auto sdxiov = static_cast<amrex::Real>((x_new - x_old) == 0. ? 1. : sdxi/(x_new - x_old));
                Exp += Ex_arr(lo.x+i_E_new-1+i, lo.y+j_E_new-1+j, lo.z+k_E_new-1+k)*sdxiov*sdzjk;
            }
        }
    }
    for (int k=dkl_E; k<=depos_order+2-dku_E; k++) {
        for (int i=dil_E; i<=depos_order+2-diu_E; i++) {
            amrex::Real sdyik = one_third*(sx_E_new[i]*sz_E_new[k] + sx_E_old[i]*sz_E_old[k])
                               +one_sixth*(sx_E_new[i]*sz_E_old[k] + sx_E_old[i]*sz_E_new[k]);
            amrex::Real sdyj = 0._rt;
            for (int j=djl_E; j<=depos_order+1-dju_E; j++) {
                sdyj += (sy_E_old[j] - sy_E_new[j]);
                auto sdyjov = static_cast<amrex::Real>((y_new - y_old) == 0. ? 1. : sdyj/(y_new - y_old));
                Eyp += Ey_arr(lo.x+i_E_new-1+i, lo.y+j_E_new-1+j, lo.z+k_E_new-1+k)*sdyjov*sdyik;
            }
        }
    }
    for (int j=djl_E; j<=depos_order+2-dju_E; j++) {
        for (int i=dil_E; i<=depos_order+2-diu_E; i++) {
            amrex::Real sdzij = one_third*(sx_E_new[i]*sy_E_new[j] + sx_E_old[i]*sy_E_old[j])
                               +one_sixth*(sx_E_new[i]*sy_E_old[j] + sx_E_old[i]*sy_E_new[j]);
            amrex::Real sdzk = 0._rt;
            for (int k=dkl_E; k<=depos_order+1-dku_E; k++) {
                sdzk += (sz_E_old[k] - sz_E_new[k]);
                auto sdzkov = static_cast<amrex::Real>((z_new - z_old) == 0. ? 1. : sdzk/(z_new - z_old));
                Ezp += Ez_arr(lo.x+i_E_new-1+i, lo.y+j_E_new-1+j, lo.z+k_E_new-1+k)*sdzkov*sdzij;
            }
        }
    }
    for (int k=dkl_B; k<=depos_order+2-dku_B; k++) {
        for (int j=djl_B; j<=depos_order+2-dju_B; j++) {
            amrex::Real sdzjk = one_third*(sy_B_new[j]*sz_B_new[k] + sy_B_old[j]*sz_B_old[k])
                               +one_sixth*(sy_B_new[j]*sz_B_old[k] + sy_B_old[j]*sz_B_new[k]);
            amrex::Real sdxi = 0._rt;
            for (int i=dil_B; i<=depos_order+1-diu_B; i++) {
                sdxi += (sx_B_old[i] - sx_B_new[i]);
                auto sdxiov = static_cast<amrex::Real>((x_new - x_old) == 0. ? 1. : sdxi/(x_new - x_old));
                Bxp += Bx_arr(lo.x+i_B_new-1+i, lo.y+j_B_new-1+j, lo.z+k_B_new-1+k)*sdxiov*sdzjk;
            }
        }
    }
    for (int k=dkl_B; k<=depos_order+2-dku_B; k++) {
        for (int i=dil_B; i<=depos_order+2-diu_B; i++) {
            amrex::Real sdyik = one_third*(sx_B_new[i]*sz_B_new[k] + sx_B_old[i]*sz_B_old[k])
                               +one_sixth*(sx_B_new[i]*sz_B_old[k] + sx_B_old[i]*sz_B_new[k]);
            amrex::Real sdyj = 0._rt;
            for (int j=djl_B; j<=depos_order+1-dju_B; j++) {
                sdyj += (sy_B_old[j] - sy_B_new[j]);
                auto sdyjov = static_cast<amrex::Real>((y_new - y_old) == 0. ? 1. : sdyj/(y_new - y_old));
                Byp += By_arr(lo.x+i_B_new-1+i, lo.y+j_B_new-1+j, lo.z+k_B_new-1+k)*sdyjov*sdyik;
            }
        }
    }
    for (int j=djl_B; j<=depos_order+2-dju_B; j++) {
        for (int i=dil_B; i<=depos_order+2-diu_B; i++) {
            amrex::Real sdzij = one_third*(sx_B_new[i]*sy_B_new[j] + sx_B_old[i]*sy_B_old[j])
                               +one_sixth*(sx_B_new[i]*sy_B_old[j] + sx_B_old[i]*sy_B_new[j]);
            amrex::Real sdzk = 0._rt;
            for (int k=dkl_B; k<=depos_order+1-dku_B; k++) {
                sdzk += (sz_B_old[k] - sz_B_new[k]);
                auto sdzkov = static_cast<amrex::Real>((z_new - z_old) == 0. ? 1. : sdzk/(z_new - z_old));
                Bzp += Bz_arr(lo.x+i_B_new-1+i, lo.y+j_B_new-1+j, lo.z+k_E_new-1+k)*sdzkov*sdzij;
            }
        }
    }

#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)

    for (int k=dkl_E; k<=depos_order+2-dku_E; k++) {
        amrex::Real sdzk = 0.5_rt*(sz_E_new[k] + sz_E_old[k]);
        amrex::Real sdxi = 0._rt;
        for (int i=dil_E; i<=depos_order+1-diu_E; i++) {
            sdxi += (sx_E_old[i] - sx_E_new[i]);
            auto sdxiov = static_cast<amrex::Real>((x_new - x_old) == 0. ? 1. : sdxi/(x_new - x_old));
            Exp += Ex_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 0)*sdxiov*sdzk;
            Bzp += Bz_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 0)*sdxiov*sdzk;
        }
    }
    for (int k=dkl_E; k<=depos_order+2-dku_E; k++) {
        for (int i=dil_E; i<=depos_order+2-diu_E; i++) {
            Real const sdyj = (
                one_third*(sx_E_new[i]*sz_E_new[k] + sx_E_old[i]*sz_E_old[k])
               +one_sixth*(sx_E_new[i]*sz_E_old[k] + sx_E_old[i]*sz_E_new[k]));
            Eyp += Ey_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 0)*sdyj;
        }
    }
    for (int i=dil_E; i<=depos_order+2-diu_E; i++) {
        amrex::Real sdxi = 0.5_rt*(sx_E_new[i] + sx_E_old[i]);
        amrex::Real sdzk = 0._rt;
        for (int k=dkl_E; k<=depos_order+1-dku_E; k++) {
            sdzk += (sz_E_old[k] - sz_E_new[k]);
            auto sdzkov = static_cast<amrex::Real>((z_new - z_old) == 0. ? 1. : sdzk/(z_new - z_old));
            Ezp += Ez_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 0)*sdzkov*sdxi;
            Bxp += Bx_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 0)*sdzkov*sdxi;
        }
    }
    for (int k=dkl_By; k<=depos_order+1-dku_By; k++) {
        for (int i=dil_By; i<=depos_order+1-diu_By; i++) {
            Real const sdyj = (
                one_third*(sx_By_new[i]*sz_By_new[k] + sx_By_old[i]*sz_By_old[k])
               +one_sixth*(sx_By_new[i]*sz_By_old[k] + sx_By_old[i]*sz_By_new[k]));
            Byp += By_arr(lo.x+i_By_new-1+i, lo.y+k_By_new-1+k, 0, 0)*sdyj;
        }
    }

#ifdef WARPX_DIM_RZ
    Complex xy_mid = xy_mid0;

    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {

        // Gather field on particle Exp from field on grid ex_arr
        // Gather field on particle Bzp from field on grid bz_arr
        for (int k=dkl_E; k<=depos_order+2-dku_E; k++) {
            amrex::Real sdzk = 0.5_rt*(sz_E_new[k] + sz_E_old[k]);
            amrex::Real sdxi = 0._rt;
            for (int i=dil_E; i<=depos_order+1-diu_E; i++) {
                sdxi += (sx_E_old[i] - sx_E_new[i]);
                auto sdxiov = static_cast<amrex::Real>((x_new - x_old) == 0. ? 1. : sdxi/(x_new - x_old));
                const amrex::Real dEx = (+ Ex_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode-1)*xy_mid.real()
                                         - Ex_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode)*xy_mid.imag());
                const amrex::Real dBz = (+ Bz_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode-1)*xy_mid.real()
                                         - Bz_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode)*xy_mid.imag());
                Exp += dEx*sdxiov*sdzk;
                Bzp += dBz*sdxiov*sdzk;
            }
        }
        // Gather field on particle Eyp from field on grid ey_arr
        for (int k=dkl_E; k<=depos_order+2-dku_E; k++) {
            for (int i=dil_E; i<=depos_order+2-diu_E; i++) {
                Real const sdyj = (
                    one_third*(sx_E_new[i]*sz_E_new[k] + sx_E_old[i]*sz_E_old[k])
                   +one_sixth*(sx_E_new[i]*sz_E_old[k] + sx_E_old[i]*sz_E_new[k]));
                const amrex::Real dEy = (+ Ey_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode-1)*xy_mid.real()
                                         - Ey_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode)*xy_mid.imag());
                Eyp += dEy*sdyj;
            }
        }
        // Gather field on particle Ezp from field on grid ez_arr
        // Gather field on particle Bxp from field on grid bx_arr
        for (int i=dil_E; i<=depos_order+2-diu_E; i++) {
            amrex::Real sdxi = 0.5_rt*(sx_E_new[i] + sx_E_old[i]);
            amrex::Real sdzk = 0._rt;
            for (int k=dkl_E; k<=depos_order+1-dku_E; k++) {
                sdzk += (sz_E_old[k] - sz_E_new[k]);
                auto sdzkov = static_cast<amrex::Real>((z_new - z_old) == 0. ? 1. : sdzk/(z_new - z_old));
                const amrex::Real dEz = (+ Ez_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode-1)*xy_mid.real()
                                         - Ez_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode)*xy_mid.imag());
                const amrex::Real dBx = (+ Bx_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode-1)*xy_mid.real()
                                         - Bx_arr(lo.x+i_E_new-1+i, lo.y+k_E_new-1+k, 0, 2*imode)*xy_mid.imag());
                Ezp += dEz*sdzkov*sdxi;
                Bxp += dBx*sdzkov*sdxi;
            }
        }
        // Gather field on particle Byp from field on grid by_arr
        for (int k=dkl_By; k<=depos_order+1-dku_By; k++) {
            for (int i=dil_By; i<=depos_order+1-diu_By; i++) {
                Real const sdyj = (
                    one_third*(sx_By_new[i]*sz_By_new[k] + sx_By_old[i]*sz_By_old[k])
                   +one_sixth*(sx_By_new[i]*sz_By_old[k] + sx_By_old[i]*sz_By_new[k]));
                const amrex::Real dBy = (+ By_arr(lo.x+i_By_new-1+i, lo.y+k_By_new-1+k, 0, 2*imode-1)*xy_mid.real()
                                         - By_arr(lo.x+i_By_new-1+i, lo.y+k_By_new-1+k, 0, 2*imode)*xy_mid.imag());
                Byp += dBy*sdyj;
            }
        }
        xy_mid = xy_mid*xy_mid0;
    }

    // Convert Exp and Eyp (which are actually Er and Etheta) to Ex and Ey
    const amrex::Real Exp_save = Exp;
    Exp = costheta_mid*Exp - sintheta_mid*Eyp;
    Eyp = costheta_mid*Eyp + sintheta_mid*Exp_save;
    const amrex::Real Bxp_save = Bxp;
    Bxp = costheta_mid*Bxp - sintheta_mid*Byp;
    Byp = costheta_mid*Byp + sintheta_mid*Bxp_save;

#endif

#elif defined(WARPX_DIM_1D_Z)

    for (int k=dkl_E; k<=depos_order+2-dku_E; k++) {
        amrex::Real const sdzk = 0.5_rt*(sz_E_old[k] + sz_E_new[k]);
        Exp += Ex_arr(lo.x+k_E_new-1+k, 0, 0, 0)*sdzk;
        Eyp += Ey_arr(lo.x+k_E_new-1+k, 0, 0, 0)*sdzk;
        Bzp += Bz_arr(lo.x+k_E_new-1+k, 0, 0, 0)*sdzk;
    }
    amrex::Real sdzk = 0._rt;
    for (int k=dkl_E; k<=depos_order+1-dku_E; k++) {
        sdzk += (sz_E_old[k] - sz_E_new[k]);
        auto sdzkov = static_cast<amrex::Real>((z_new - z_old) == 0. ? 1. : sdzk/(z_new - z_old));
        Bxp += Bx_arr(lo.x+k_E_new-1+k, 0, 0, 0)*sdzkov;
        Byp += By_arr(lo.x+k_E_new-1+k, 0, 0, 0)*sdzkov;
        Ezp += Ez_arr(lo.x+k_E_new-1+k, 0, 0, 0)*sdzkov;
    }
#endif
}

/**
 * \brief Field gather for particles
 *
 * \tparam depos_order         deposition order
 * \tparam lower_in_v          lower shape order in parallel direction (Galerkin)
 * \param getPosition          A functor for returning the particle position.
 * \param getExternalEB        A functor for assigning the external E and B fields.
 * \param Exp,Eyp,Ezp          Pointer to array of electric field on particles.
 * \param Bxp,Byp,Bzp          Pointer to array of magnetic field on particles.
 * \param exfab,eyfab,ezfab    Array4 of the electric field, either full array or tile.
 * \param bxfab,byfab,bzfab    Array4 of the magnetic field, either full array or tile.
 * \param np_to_gather         Number of particles for which field is gathered.
 * \param dx                   3D cell size
 * \param xyzmin               Physical lower bounds of domain.
 * \param lo                   Index lower bounds of domain.
 * \param n_rz_azimuthal_modes Number of azimuthal modes when using RZ geometry
 */
template <int depos_order, int lower_in_v>
void doGatherShapeN(const GetParticlePosition<PIdx>& getPosition,
                    const GetExternalEBField& getExternalEB,
                    amrex::ParticleReal * const Exp, amrex::ParticleReal * const Eyp,
                    amrex::ParticleReal * const Ezp, amrex::ParticleReal * const Bxp,
                    amrex::ParticleReal * const Byp, amrex::ParticleReal * const Bzp,
                    amrex::FArrayBox const * const exfab,
                    amrex::FArrayBox const * const eyfab,
                    amrex::FArrayBox const * const ezfab,
                    amrex::FArrayBox const * const bxfab,
                    amrex::FArrayBox const * const byfab,
                    amrex::FArrayBox const * const bzfab,
                    const long np_to_gather,
                    const std::array<amrex::Real, 3>& dx,
                    const std::array<amrex::Real, 3> xyzmin,
                    const amrex::Dim3 lo,
                    const int n_rz_azimuthal_modes)
{

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

    // Loop over particles and gather fields from
    // {e,b}{x,y,z}_arr to {E,B}{xyz}p.
    amrex::ParallelFor(
        np_to_gather,
        [=] AMREX_GPU_DEVICE (long ip) {

            amrex::ParticleReal xp, yp, zp;
            getPosition(ip, xp, yp, zp);
            getExternalEB(ip, Exp[ip], Eyp[ip], Ezp[ip], Bxp[ip], Byp[ip], Bzp[ip]);

            doGatherShapeN<depos_order, lower_in_v>(
                xp, yp, zp, Exp[ip], Eyp[ip], Ezp[ip], Bxp[ip], Byp[ip], Bzp[ip],
                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        }
        );
}

/**
 * \brief Field gather for a single particle
 *
 * \param xp,yp,zp                Particle position coordinates
 * \param Exp,Eyp,Ezp             Electric field on particles.
 * \param Bxp,Byp,Bzp             Magnetic field on particles.
 * \param ex_arr,ey_arr,ez_arr    Array4 of the electric field, either full array or tile.
 * \param bx_arr,by_arr,bz_arr    Array4 of the magnetic field, either full array or tile.
 * \param ex_type,ey_type,ez_type IndexType of the electric field
 * \param bx_type,by_type,bz_type IndexType of the magnetic field
 * \param dx_arr                  3D cell spacing
 * \param xyzmin_arr              Physical lower bounds of domain in x, y, z.
 * \param lo                      Index lower bounds of domain.
 * \param n_rz_azimuthal_modes    Number of azimuthal modes when using RZ geometry
 * \param nox                     order of the particle shape function
 * \param galerkin_interpolation  whether to use lower order in v
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void doGatherShapeN (const amrex::ParticleReal xp,
                     const amrex::ParticleReal yp,
                     const amrex::ParticleReal zp,
                     amrex::ParticleReal& Exp,
                     amrex::ParticleReal& Eyp,
                     amrex::ParticleReal& Ezp,
                     amrex::ParticleReal& Bxp,
                     amrex::ParticleReal& Byp,
                     amrex::ParticleReal& Bzp,
                     amrex::Array4<amrex::Real const> const& ex_arr,
                     amrex::Array4<amrex::Real const> const& ey_arr,
                     amrex::Array4<amrex::Real const> const& ez_arr,
                     amrex::Array4<amrex::Real const> const& bx_arr,
                     amrex::Array4<amrex::Real const> const& by_arr,
                     amrex::Array4<amrex::Real const> const& bz_arr,
                     const amrex::IndexType ex_type,
                     const amrex::IndexType ey_type,
                     const amrex::IndexType ez_type,
                     const amrex::IndexType bx_type,
                     const amrex::IndexType by_type,
                     const amrex::IndexType bz_type,
                     const amrex::GpuArray<amrex::Real, 3>& dx_arr,
                     const amrex::GpuArray<amrex::Real, 3>& xyzmin_arr,
                     const amrex::Dim3& lo,
                     const int n_rz_azimuthal_modes,
                     const int nox,
                     const bool galerkin_interpolation)
{
    if (galerkin_interpolation) {
        if (nox == 1) {
            doGatherShapeN<1,1>(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 2) {
            doGatherShapeN<2,1>(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 3) {
            doGatherShapeN<3,1>(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        }
    } else {
        if (nox == 1) {
            doGatherShapeN<1,0>(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 2) {
            doGatherShapeN<2,0>(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 3) {
            doGatherShapeN<3,0>(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        }
    }
}


/**
 * \brief Field gather for a single particle
 *
 * \param xp_n,yp_n,zp_n          Particle position coordinates at start of step
 * \param xp_nph,yp_nph,zp_nph    Particle position coordinates at half time level (n + half)
 * \param Exp,Eyp,Ezp             Electric field on particles.
 * \param Bxp,Byp,Bzp             Magnetic field on particles.
 * \param ex_arr,ey_arr,ez_arr    Array4 of the electric field, either full array or tile.
 * \param bx_arr,by_arr,bz_arr    Array4 of the magnetic field, either full array or tile.
 * \param ex_type,ey_type,ez_type IndexType of the electric field
 * \param bx_type,by_type,bz_type IndexType of the magnetic field
 * \param dx_arr                  3D cell spacing
 * \param xyzmin_arr              Physical lower bounds of domain in x, y, z.
 * \param lo                      Index lower bounds of domain.
 * \param n_rz_azimuthal_modes    Number of azimuthal modes when using RZ geometry
 * \param nox                     order of the particle shape function
 * \param galerkin_interpolation  whether to use lower order in v
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void doGatherShapeNImplicit (
                     const amrex::ParticleReal xp_n,
                     const amrex::ParticleReal yp_n,
                     const amrex::ParticleReal zp_n,
                     const amrex::ParticleReal xp_nph,
                     const amrex::ParticleReal yp_nph,
                     const amrex::ParticleReal zp_nph,
                     amrex::ParticleReal& Exp,
                     amrex::ParticleReal& Eyp,
                     amrex::ParticleReal& Ezp,
                     amrex::ParticleReal& Bxp,
                     amrex::ParticleReal& Byp,
                     amrex::ParticleReal& Bzp,
                     amrex::Array4<amrex::Real const> const& ex_arr,
                     amrex::Array4<amrex::Real const> const& ey_arr,
                     amrex::Array4<amrex::Real const> const& ez_arr,
                     amrex::Array4<amrex::Real const> const& bx_arr,
                     amrex::Array4<amrex::Real const> const& by_arr,
                     amrex::Array4<amrex::Real const> const& bz_arr,
                     const amrex::IndexType ex_type,
                     const amrex::IndexType ey_type,
                     const amrex::IndexType ez_type,
                     const amrex::IndexType bx_type,
                     const amrex::IndexType by_type,
                     const amrex::IndexType bz_type,
                     const amrex::GpuArray<amrex::Real, 3>& dx_arr,
                     const amrex::GpuArray<amrex::Real, 3>& xyzmin_arr,
                     const amrex::Dim3& lo,
                     const int n_rz_azimuthal_modes,
                     const int nox,
                     const bool galerkin_interpolation)
{
    if (galerkin_interpolation) {
        if (nox == 1) {
            doGatherShapeNEsirkepovStencilImplicit<1>(xp_n, yp_n, zp_n, xp_nph, yp_nph, zp_nph,
                                                      Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                                      ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                                      ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                                      dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 2) {
            doGatherShapeNEsirkepovStencilImplicit<2>(xp_n, yp_n, zp_n, xp_nph, yp_nph, zp_nph,
                                                      Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                                      ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                                      ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                                      dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 3) {
            doGatherShapeNEsirkepovStencilImplicit<3>(xp_n, yp_n, zp_n, xp_nph, yp_nph, zp_nph,
                                                      Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                                      ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                                      ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                                      dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        }
    } else {
        if (nox == 1) {
            doGatherShapeN<1,0>(xp_nph, yp_nph, zp_nph, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 2) {
            doGatherShapeN<2,0>(xp_nph, yp_nph, zp_nph, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        } else if (nox == 3) {
            doGatherShapeN<3,0>(xp_nph, yp_nph, zp_nph, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                                dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        }
    }
}

#endif // FIELDGATHER_H_
