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

#if (AMREX_SPACEDIM == 2 || AMREX_SPACEDIM == 3)
    const amrex::Real dxi = 1.0_rt/dx[0];
#endif
    const amrex::Real dzi = 1.0_rt/dx[2];
#if (AMREX_SPACEDIM == 3)
    const amrex::Real dyi = 1.0_rt/dx[1];
#endif

#if (AMREX_SPACEDIM == 2 || AMREX_SPACEDIM == 3)
    const amrex::Real xmin = xyzmin[0];
#endif
#if (AMREX_SPACEDIM == 3)
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real zmin = xyzmin[2];

    constexpr int zdir = (AMREX_SPACEDIM - 1);
    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // --- Compute shape factors

    Compute_shape_factor< depos_order > const compute_shape_factor;
    Compute_shape_factor<depos_order - galerkin_interpolation > const compute_shape_factor_galerkin;

#if (AMREX_SPACEDIM == 2 || AMREX_SPACEDIM == 3)
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

#if (AMREX_SPACEDIM == 3)
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

#if (AMREX_SPACEDIM == 1)
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

#elif (AMREX_SPACEDIM == 2)
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

#ifdef WARPX_DIM_RZ

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

        // Gather field on particle Eyp from field on grid ey_arr
        for (int iz=0; iz<=depos_order; iz++){
            for (int ix=0; ix<=depos_order; ix++){
                const amrex::Real dEy = (+ ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 2*imode-1)*xy.real()
                                         - ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 2*imode)*xy.imag());
                Eyp += sx_ey[ix]*sz_ey[iz]*dEy;
            }
        }
        // Gather field on particle Exp from field on grid ex_arr
        // Gather field on particle Bzp from field on grid bz_arr
        for (int iz=0; iz<=depos_order; iz++){
            for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
                const amrex::Real dEx = (+ ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 2*imode-1)*xy.real()
                                         - ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 2*imode)*xy.imag());
                Exp += sx_ex[ix]*sz_ex[iz]*dEx;
                const amrex::Real dBz = (+ bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 2*imode-1)*xy.real()
                                         - bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 2*imode)*xy.imag());
                Bzp += sx_bz[ix]*sz_bz[iz]*dBz;
            }
        }
        // Gather field on particle Ezp from field on grid ez_arr
        // Gather field on particle Bxp from field on grid bx_arr
        for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
            for (int ix=0; ix<=depos_order; ix++){
                const amrex::Real dEz = (+ ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 2*imode-1)*xy.real()
                                         - ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 2*imode)*xy.imag());
                Ezp += sx_ez[ix]*sz_ez[iz]*dEz;
                const amrex::Real dBx = (+ bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 2*imode-1)*xy.real()
                                         - bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 2*imode)*xy.imag());
                Bxp += sx_bx[ix]*sz_bx[iz]*dBx;
            }
        }
        // Gather field on particle Byp from field on grid by_arr
        for (int iz=0; iz<=depos_order-galerkin_interpolation; iz++){
            for (int ix=0; ix<=depos_order-galerkin_interpolation; ix++){
                const amrex::Real dBy = (+ by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 2*imode-1)*xy.real()
                                         - by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 2*imode)*xy.imag());
                Byp += sx_by[ix]*sz_by[iz]*dBy;
            }
        }
        xy = xy*xy0;
    }

    // Convert Exp and Eyp (which are actually Er and Etheta) to Ex and Ey
    const amrex::Real Exp_save = Exp;
    Exp = costheta*Exp - sintheta*Eyp;
    Eyp = costheta*Eyp + sintheta*Exp_save;
    const amrex::Real Bxp_save = Bxp;
    Bxp = costheta*Bxp - sintheta*Byp;
    Byp = costheta*Byp + sintheta*Bxp_save;
#endif

#else // (AMREX_SPACEDIM == 3)
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
 * \brief Field gather for particles
 *
 * \tparam depos_order         deposition order
 * \tparam lower_in_v          lower shape order in parallel direction (Galerkin)
 * \param getPosition          A functor for returning the particle position.
 * \param getExternalE         A functor for assigning the external E field.
 * \param getExternalB         A functor for assigning the external B field.
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
void doGatherShapeN(const GetParticlePosition& getPosition,
                    const GetExternalEField& getExternalE, const GetExternalBField& getExternalB,
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
            getExternalE(ip, Exp[ip], Eyp[ip], Ezp[ip]);
            getExternalB(ip, Bxp[ip], Byp[ip], Bzp[ip]);

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

// TODO Move this to more appropriate location?
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_centering_coefficients(amrex::GpuArray<amrex::Real, 8>& cc,
                                    const amrex::Real r, const int m)
{
    using namespace amrex;

    if (m == 4)
    {
        const amrex::Real r2 = r*r;
        const amrex::Real r3 = r*r2;
        cc[0] = - r3/6. + r/6.;
        cc[1] = r3/2. + r2/2. - r;
        cc[2] = - r3/2. - r2 + r/2. + 1.;
        cc[3] = r3/6. + r2/2. + r/3.;
    }
    else if (m == 6)
    {
        const amrex::Real r2 = r*r;
        const amrex::Real r3 = r*r2;
        const amrex::Real r4 = r*r3;
        const amrex::Real r5 = r*r4;
        cc[0] = - r5/120. + r3/24. - r/30.;
        cc[1] = r5/24. + r4/24. - r3*7./24. - r2/24. + r/4.;
        cc[2] = - r5/12. - r4/6. + r3*7./12. + r2*2./3. - r;
        cc[3] = r5/12. + r4/4. - r3*5./12. - r2*5./4. + r/3. + 1.;
        cc[4] = - r5/24. - r4/6. + r3/24. + r2*2./3. + r/2.;
        cc[5] = r5/120. + r4/24. + r3/24. - r2/24. - r/20.;
    }
    else if (m == 8)
    {
        const amrex::Real r2 = r*r;
        const amrex::Real r3 = r*r2;
        const amrex::Real r4 = r*r3;
        const amrex::Real r5 = r*r4;
        const amrex::Real r6 = r*r5;
        const amrex::Real r7 = r*r6;
        cc[0] = - r7/5040. + r5/360. - r3*7./720. + r/140.;
        cc[1] = r7/720. + r6/720. - r5*17./720. - r4/144. + r3*4./45. + r2/180. - r/15.;
        cc[2] = - r7/240. - r6/120. + r5*3./40. + r4/12. - r3*89./240. - r2*3./40. + r*3./10.;
        cc[3] = r7/144. + r6/48. - r5*17./144. - r4*13./48. + r3*11./18. + r2*3./4. - r;
        cc[4] = - r7/144. - r6/36. + r5*7./72. + r4*7./18. - r3*49./144. - r2*49./36. + r/4. + 1.;
        cc[5] = r7/240. + r6/48. - r5*3./80. - r4*13./48. - r3/15. + r2*3./4. + r*3./5.;
        cc[6] = - r7/720. - r6/120. + r5/360. + r4/12. + r3*71./720. - r2*3./40. - r/10.;
        cc[7] = r7/5040. + r6/720. + r5/720. - r4/144. - r3/90. + r2/180. + r/105.;
    }
    else
    {
        amrex::Abort("Centering coefficients implemented only for orders 4,6,8");
    }
}

// TODO Add Doxygen documentation
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void doGatherFiniteCenteringSinglePart (const amrex::ParticleReal xp,
                                        const amrex::ParticleReal yp,
                                        const amrex::ParticleReal zp,
                                        amrex::ParticleReal& Exp,
                                        amrex::ParticleReal& Eyp,
                                        amrex::ParticleReal& Ezp,
                                        amrex::ParticleReal& Bxp,
                                        amrex::ParticleReal& Byp,
                                        amrex::ParticleReal& Bzp,
                                        const amrex::Array4<amrex::Real const>& ex_arr,
                                        const amrex::Array4<amrex::Real const>& ey_arr,
                                        const amrex::Array4<amrex::Real const>& ez_arr,
                                        const amrex::Array4<amrex::Real const>& bx_arr,
                                        const amrex::Array4<amrex::Real const>& by_arr,
                                        const amrex::Array4<amrex::Real const>& bz_arr,
                                        const amrex::GpuArray<amrex::Real,3>& dxyz,
                                        const amrex::GpuArray<amrex::Real,3>& xyz_min,
                                        const amrex::Dim3& lo,
                                        const int nox,
                                        const int noy,
                                        const int noz)
{
    using namespace amrex;

#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(yp, noy);
#endif

    // Pad input arrays with zeros beyond ghost cells
    // for out-of-bound accesses due to large-stencil operations
    const auto ex_arr_zeropad = [ex_arr] (const int ii, const int jj, const int kk) noexcept
    {
        return ex_arr.contains(ii,jj,kk) ? ex_arr(ii,jj,kk) : 0.0_rt;
    };
    const auto ey_arr_zeropad = [ey_arr] (const int ii, const int jj, const int kk) noexcept
    {
        return ey_arr.contains(ii,jj,kk) ? ey_arr(ii,jj,kk) : 0.0_rt;
    };
    const auto ez_arr_zeropad = [ez_arr] (const int ii, const int jj, const int kk) noexcept
    {
        return ez_arr.contains(ii,jj,kk) ? ez_arr(ii,jj,kk) : 0.0_rt;
    };
    const auto bx_arr_zeropad = [bx_arr] (const int ii, const int jj, const int kk) noexcept
    {
        return bx_arr.contains(ii,jj,kk) ? bx_arr(ii,jj,kk) : 0.0_rt;
    };
    const auto by_arr_zeropad = [by_arr] (const int ii, const int jj, const int kk) noexcept
    {
        return by_arr.contains(ii,jj,kk) ? by_arr(ii,jj,kk) : 0.0_rt;
    };
    const auto bz_arr_zeropad = [bz_arr] (const int ii, const int jj, const int kk) noexcept
    {
        return bz_arr.contains(ii,jj,kk) ? bz_arr(ii,jj,kk) : 0.0_rt;
    };

    const amrex::Real x = (xp - xyz_min[0]) / dxyz[0];
#if (AMREX_SPACEDIM == 3)
    const amrex::Real y = (yp - xyz_min[1]) / dxyz[1];
#endif
    const amrex::Real z = (zp - xyz_min[2]) / dxyz[2];

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

    compute_centering_coefficients(sci, ri, noi);
    compute_centering_coefficients(scj, rj, noj);
#if (AMREX_SPACEDIM == 3)
    compute_centering_coefficients(sck, rk, nok);
#endif

    amrex::Real ci = 1.0_rt;
    amrex::Real cj = 1.0_rt;
    amrex::Real ck = 1.0_rt;

    Exp = 0._rt;
    Eyp = 0._rt;
    Ezp = 0._rt;
    Bxp = 0._rt;
    Byp = 0._rt;
    Bzp = 0._rt;

    for (int kk = 0; kk <= nk; kk++)
    {
#if (AMREX_SPACEDIM == 3)
        ck = sck[kk];
#endif
        for (int jj = 0; jj <= nj; jj++)
        {
            cj = scj[jj];

            for (int ii = 0; ii <= ni; ii++)
            {
                ci = sci[ii];

                Exp += ci * cj * ck * ex_arr_zeropad(lo.x+imin+ii,lo.y+jmin+jj,lo.z+kmin+kk);
                Eyp += ci * cj * ck * ey_arr_zeropad(lo.x+imin+ii,lo.y+jmin+jj,lo.z+kmin+kk);
                Ezp += ci * cj * ck * ez_arr_zeropad(lo.x+imin+ii,lo.y+jmin+jj,lo.z+kmin+kk);

                Bxp += ci * cj * ck * bx_arr_zeropad(lo.x+imin+ii,lo.y+jmin+jj,lo.z+kmin+kk);
                Byp += ci * cj * ck * by_arr_zeropad(lo.x+imin+ii,lo.y+jmin+jj,lo.z+kmin+kk);
                Bzp += ci * cj * ck * bz_arr_zeropad(lo.x+imin+ii,lo.y+jmin+jj,lo.z+kmin+kk);
            }
        }
    }
}

// TODO Add Doxygen documentation
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void doGatherFiniteCentering (const amrex::ParticleReal xp,
                              const amrex::ParticleReal yp,
                              const amrex::ParticleReal zp,
                              amrex::ParticleReal& Exp,
                              amrex::ParticleReal& Eyp,
                              amrex::ParticleReal& Ezp,
                              amrex::ParticleReal& Bxp,
                              amrex::ParticleReal& Byp,
                              amrex::ParticleReal& Bzp,
                              const amrex::Array4<amrex::Real const>& ex_arr,
                              const amrex::Array4<amrex::Real const>& ey_arr,
                              const amrex::Array4<amrex::Real const>& ez_arr,
                              const amrex::Array4<amrex::Real const>& bx_arr,
                              const amrex::Array4<amrex::Real const>& by_arr,
                              const amrex::Array4<amrex::Real const>& bz_arr,
                              const amrex::GpuArray<amrex::Real,3>& dxyz,
                              const amrex::GpuArray<amrex::Real,3>& xyz_min,
                              const amrex::Dim3& lo,
                              const int nox,
                              const int noy,
                              const int noz,
                              const int shape_nox,
                              const int shape_noy,
                              const int shape_noz)
{
    using namespace amrex::literals;

#ifdef WARPX_DIM_XZ
    amrex::ignore_unused(shape_noy);
#endif

    // Use only original particle position for linear shape factors
    int nxp = 1;
    int nyp = 1;
    int nzp = 1;

    // Use three particle positions for quadratic and cubic shape factors
    if (shape_nox == 2 || shape_nox == 3) nxp = 3;
#ifdef WARPX_DIM_3D
    if (shape_noy == 2 || shape_noy == 3) nyp = 3;
#endif
    if (shape_noz == 2 || shape_noz == 3) nzp = 3;

    amrex::Real xxp[3];
    amrex::Real yyp[3];
    amrex::Real zzp[3];

    xxp[0] = xp;
    yyp[0] = yp;
    zzp[0] = zp;

    amrex::Real sx[3];
    amrex::Real sy[3];
    amrex::Real sz[3];

    sx[0] = 1._rt;
    sy[0] = 1._rt;
    sz[0] = 1._rt;

    if (shape_nox == 2 || shape_nox == 3)
    {
        // dxyz[0] is dx
        xxp[0] = xp - dxyz[0];
        xxp[1] = xp;
        xxp[2] = xp + dxyz[0];

        // Spline coefficients for a spline centered in xp
        sx[0] = (shape_nox == 2) ? 1._rt / 8._rt : 1._rt / 6._rt;
        sx[1] = (shape_nox == 2) ? 3._rt / 4._rt : 2._rt / 3._rt;
        sx[2] = (shape_nox == 2) ? 1._rt / 8._rt : 1._rt / 6._rt;
    }

#ifdef WARPX_DIM_3D
    if (shape_noy == 2 || shape_noy == 3)
    {
        // dxyz[1] is dy
        yyp[0] = yp - dxyz[1];
        yyp[1] = yp;
        yyp[2] = yp + dxyz[1];

        // Spline coefficients for a spline centered in yp
        sy[0] = (shape_noy == 2) ? 1._rt / 8._rt : 1._rt / 6._rt;
        sy[1] = (shape_noy == 2) ? 3._rt / 4._rt : 2._rt / 3._rt;
        sy[2] = (shape_noy == 2) ? 1._rt / 8._rt : 1._rt / 6._rt;
    }
#endif

    if (shape_noz == 2 || shape_noz == 3)
    {
        // dxyz[2] is dz
        zzp[0] = zp - dxyz[2];
        zzp[1] = zp;
        zzp[2] = zp + dxyz[2];

        // Spline coefficients for a spline centered in zp
        sz[0] = (shape_noz == 2) ? 1._rt / 8._rt : 1._rt / 6._rt;
        sz[1] = (shape_noz == 2) ? 3._rt / 4._rt : 2._rt / 3._rt;
        sz[2] = (shape_noz == 2) ? 1._rt / 8._rt : 1._rt / 6._rt;
    }

    amrex::ParticleReal Exp_tmp, Eyp_tmp, Ezp_tmp;
    amrex::ParticleReal Bxp_tmp, Byp_tmp, Bzp_tmp;

    for (int iz = 0; iz < nzp; iz++)
    {
        for (int iy = 0; iy < nyp; iy++)
        {
            for (int ix = 0; ix < nxp; ix++)
            {
                doGatherFiniteCenteringSinglePart(xxp[ix], yyp[iy], zzp[iz],
                    Exp_tmp, Eyp_tmp, Ezp_tmp, Bxp_tmp, Byp_tmp, Bzp_tmp,
                    ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                    dxyz, xyz_min, lo, nox, noy, noz);

                Exp_tmp *= sx[ix] * sy[iy] * sz[iz];
                Eyp_tmp *= sx[ix] * sy[iy] * sz[iz];
                Ezp_tmp *= sx[ix] * sy[iy] * sz[iz];

                Bxp_tmp *= sx[ix] * sy[iy] * sz[iz];
                Byp_tmp *= sx[ix] * sy[iy] * sz[iz];
                Bzp_tmp *= sx[ix] * sy[iy] * sz[iz];

                Exp += Exp_tmp;
                Eyp += Eyp_tmp;
                Ezp += Ezp_tmp;
                Bxp += Bxp_tmp;
                Byp += Byp_tmp;
                Bzp += Bzp_tmp;
            }
        }
    }
}

#endif // FIELDGATHER_H_
