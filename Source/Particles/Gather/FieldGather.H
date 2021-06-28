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
 * \param xp, yp, zp                Particle position coordinates
 * \param Exp, Eyp, Ezp             Electric field on particles.
 * \param Bxp, Byp, Bzp             Magnetic field on particles.
 * \param ex_arr ey_arr ez_arr      Array4 of the electric field, either full array or tile.
 * \param bx_arr by_arr bz_arr      Array4 of the magnetic field, either full array or tile.
 * \param ex_type, ey_type, ez_type IndexType of the electric field
 * \param bx_type, by_type, bz_type IndexType of the magnetic field
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

#ifndef WARPX_DIM_RZ
    amrex::ignore_unused(n_rz_azimuthal_modes);
#endif

    const amrex::Real dxi = 1.0_rt/dx[0];
    const amrex::Real dzi = 1.0_rt/dx[2];
#if (AMREX_SPACEDIM == 3)
    const amrex::Real dyi = 1.0_rt/dx[1];
#endif

    const amrex::Real xmin = xyzmin[0];
#if (AMREX_SPACEDIM == 3)
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real zmin = xyzmin[2];

    constexpr int zdir = (AMREX_SPACEDIM - 1);
    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // --- Compute shape factors
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
    Compute_shape_factor< depos_order > const compute_shape_factor;
    Compute_shape_factor<depos_order - galerkin_interpolation > const compute_shape_factor_galerkin;
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
#if (AMREX_SPACEDIM == 2)
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
 * /param getPosition          : A functor for returning the particle position.
 * /param getExternalEField    : A functor for assigning the external E field.
 * /param getExternalBField    : A functor for assigning the external B field.
 * \param Exp, Eyp, Ezp        : Pointer to array of electric field on particles.
 * \param Bxp, Byp, Bzp        : Pointer to array of magnetic field on particles.
 * \param exfab eyfab ezfab    : Array4 of the electric field, either full array or tile.
 * \param ezfab bxfab bzfab    : Array4 of the magnetic field, either full array or tile.
 * \param np_to_gather         : Number of particles for which field is gathered.
 * \param dx                   : 3D cell size
 * \param xyzmin               : Physical lower bounds of domain.
 * \param lo                   : Index lower bounds of domain.
 * \param n_rz_azimuthal_modes : Number of azimuthal modes when using RZ geometry
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
 * /param xp, yp, zp                : Particle position coordinates
 * \param Exp, Eyp, Ezp             : Electric field on particles.
 * \param Bxp, Byp, Bzp             : Magnetic field on particles.
 * \param ex_arr ey_arr ez_arr      : Array4 of the electric field, either full array or tile.
 * \param bx_arr by_arr bz_arr      : Array4 of the magnetic field, either full array or tile.
 * \param ex_type, ey_type, ez_type : IndexType of the electric field
 * \param bx_type, by_type, bz_type : IndexType of the magnetic field
 * \param dx                        : 3D cell spacing
 * \param xyzmin                    : Physical lower bounds of domain in x, y, z.
 * \param lo                        : Index lower bounds of domain.
 * \param n_rz_azimuthal_modes      : Number of azimuthal modes when using RZ geometry
 * \param nox                       : order of the particle shape function
 * \param galerkin_interpolation        : whether to use lower order in v
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

#endif // FIELDGATHER_H_
