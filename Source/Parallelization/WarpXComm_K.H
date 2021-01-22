/* Copyright 2019 Axel Huebl, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMM_K_H_
#define WARPX_COMM_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp (int j, int k, int l,
                   amrex::Array4<amrex::Real      > const& arr_aux,
                   amrex::Array4<amrex::Real const> const& arr_fine,
                   amrex::Array4<amrex::Real const> const& arr_coarse,
                   const amrex::IntVect& arr_stag,
                   const amrex::IntVect& rr)
{
    using namespace amrex;

    // NOTE Indices (j,k,l) in the following refer to (x,z,-) in 2D and (x,y,z) in 3D

    // Refinement ratio
    const int rj = rr[0];
    const int rk = rr[1];
    const int rl = (AMREX_SPACEDIM == 2) ? 1 : rr[2];

    // Staggering (0: cell-centered; 1: nodal)
    const int sj = arr_stag[0];
    const int sk = arr_stag[1];
    const int sl = (AMREX_SPACEDIM == 2) ? 0 : arr_stag[2];

    // Number of points used for interpolation from coarse grid to fine grid
    const int nj = (sj == 0) ? 1 : 2;
    const int nk = (sk == 0) ? 1 : 2;
    const int nl = (sl == 0) ? 1 : 2;

    const int jc = amrex::coarsen(j, rj);
    const int kc = amrex::coarsen(k, rk);
    const int lc = amrex::coarsen(l, rl);

    amrex::Real wj;
    amrex::Real wk;
    amrex::Real wl;

    // Interpolate from coarse grid to fine grid using either 1 point with weight 1, if both grids
    // are cell-centered, or 2 points with weights depending on the distance, if both grids are nodal
    amrex::Real res = 0.0_rt;
    for         (int jj = 0; jj < nj; jj++) {
        for     (int kk = 0; kk < nk; kk++) {
            for (int ll = 0; ll < nl; ll++) {
                wj = (sj == 0) ? 1.0_rt : (rj - amrex::Math::abs(j - (jc + jj) * rj))
                                          / static_cast<amrex::Real>(rj);
                wk = (sk == 0) ? 1.0_rt : (rk - amrex::Math::abs(k - (kc + kk) * rk))
                                          / static_cast<amrex::Real>(rk);
                wl = (sl == 0) ? 1.0_rt : (rl - amrex::Math::abs(l - (lc + ll) * rl))
                                          / static_cast<amrex::Real>(rl);
                res += wj * wk * wl * arr_coarse(jc+jj,kc+kk,lc+ll);
            }
        }
    }
    arr_aux(j,k,l) = arr_fine(j,k,l) + res;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bxa,
                               amrex::Array4<amrex::Real const> const& Bxf,
                               amrex::Array4<amrex::Real const> const& Bxc,
                               amrex::Array4<amrex::Real const> const& Bxg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real owx = 1.0_rt-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real owy = 1.0_rt-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Bxg(jg  ,kg  ,0)
        +     owx *  wy * Bxg(jg  ,kg+1,0)
        +      wx * owy * Bxg(jg+1,kg  ,0)
        +      wx *  wy * Bxg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wy = 0.5_rt-wy;  owy = 1.0_rt-wy;
    Real bc = owx * owy * Bxc(jg  ,kg  ,0)
        +     owx *  wy * Bxc(jg  ,kg-1,0)
        +      wx * owy * Bxc(jg+1,kg  ,0)
        +      wx *  wy * Bxc(jg+1,kg-1,0);

    // interp from fine staggered to fine nodal
    Real bf = 0.5_rt*(Bxf(j,k-1,0) + Bxf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0_rt : 0.5_rt;
    Real owz = 1.0_rt-wz;

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * owz * Bxg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bxg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Bxg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Bxg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Bxg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Bxg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Bxg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Bxg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wy = 0.5_rt-wy;  owy = 1.0_rt-wy;
    wz = 0.5_rt-wz;  owz = 1.0_rt-wz;
    Real bc = owx * owy * owz * Bxc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bxc(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Bxc(jg  ,kg-1,lg  )
        +      wx *  wy * owz * Bxc(jg+1,kg-1,lg  )
        +     owx * owy *  wz * Bxc(jg  ,kg  ,lg-1)
        +      wx * owy *  wz * Bxc(jg+1,kg  ,lg-1)
        +     owx *  wy *  wz * Bxc(jg  ,kg-1,lg-1)
        +      wx *  wy *  wz * Bxc(jg+1,kg-1,lg-1);

    // interp from fine stagged to fine nodal
    Real bf = 0.25_rt*(Bxf(j,k-1,l-1) + Bxf(j,k,l-1) + Bxf(j,k-1,l) + Bxf(j,k,l));
#endif

    Bxa(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bya,
                               amrex::Array4<amrex::Real const> const& Byf,
                               amrex::Array4<amrex::Real const> const& Byc,
                               amrex::Array4<amrex::Real const> const& Byg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real owx = 1.0_rt-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real owy = 1.0_rt-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Byg(jg  ,kg  ,0)
        +     owx *  wy * Byg(jg  ,kg+1,0)
        +      wx * owy * Byg(jg+1,kg  ,0)
        +      wx *  wy * Byg(jg+1,kg+1,0);

    // interp from coarse stagged (cell-centered for By) to fine nodal
    wx = 0.5_rt-wx;  owx = 1.0_rt-wx;
    wy = 0.5_rt-wy;  owy = 1.0_rt-wy;
    Real bc = owx * owy * Byc(jg  ,kg  ,0)
        +     owx *  wy * Byc(jg  ,kg-1,0)
        +      wx * owy * Byc(jg-1,kg  ,0)
        +      wx *  wy * Byc(jg-1,kg-1,0);

    // interp form fine stagger (cell-centered for By) to fine nodal
    Real bf = 0.25_rt*(Byf(j,k,0) + Byf(j-1,k,0) + Byf(j,k-1,0) + Byf(j-1,k-1,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0_rt : 0.5_rt;
    Real owz = 1.0_rt-wz;

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * owz * Byg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Byg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Byg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Byg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Byg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Byg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Byg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Byg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wx = 0.5_rt-wx;  owx = 1.0_rt-wx;
    wz = 0.5_rt-wz;  owz = 1.0_rt-wz;
    Real bc = owx * owy * owz * Byc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Byc(jg-1,kg  ,lg  )
        +     owx *  wy * owz * Byc(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Byc(jg-1,kg+1,lg  )
        +     owx * owy *  wz * Byc(jg  ,kg  ,lg-1)
        +      wx * owy *  wz * Byc(jg-1,kg  ,lg-1)
        +     owx *  wy *  wz * Byc(jg  ,kg+1,lg-1)
        +      wx *  wy *  wz * Byc(jg-1,kg+1,lg-1);

    // interp from fine stagged to fine nodal
    Real bf = 0.25_rt*(Byf(j-1,k,l-1) + Byf(j,k,l-1) + Byf(j-1,k,l) + Byf(j,k,l));

#endif

    Bya(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bza,
                               amrex::Array4<amrex::Real const> const& Bzf,
                               amrex::Array4<amrex::Real const> const& Bzc,
                               amrex::Array4<amrex::Real const> const& Bzg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real owx = 1.0_rt-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real owy = 1.0_rt-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Bzg(jg  ,kg  ,0)
        +     owx *  wy * Bzg(jg  ,kg+1,0)
        +      wx * owy * Bzg(jg+1,kg  ,0)
        +      wx *  wy * Bzg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wx = 0.5_rt-wx;  owx = 1.0_rt-wx;
    Real bc = owx * owy * Bzc(jg  ,kg  ,0)
        +     owx *  wy * Bzc(jg  ,kg+1,0)
        +      wx * owy * Bzc(jg-1,kg  ,0)
        +      wx *  wy * Bzc(jg-1,kg+1,0);

    // interp from fine staggered to fine nodal
    Real bf = 0.5_rt*(Bzf(j-1,k,0) + Bzf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0_rt : 0.5_rt;
    Real owz = 1.0_rt-wz;

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * owz * Bzg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bzg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Bzg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Bzg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Bzg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Bzg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Bzg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Bzg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wx = 0.5_rt-wx;  owx = 1.0_rt-wx;
    wy = 0.5_rt-wy;  owy = 1.0_rt-wy;
    Real bc = owx * owy * owz * Bzc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bzc(jg-1,kg  ,lg  )
        +     owx *  wy * owz * Bzc(jg  ,kg-1,lg  )
        +      wx *  wy * owz * Bzc(jg-1,kg-1,lg  )
        +     owx * owy *  wz * Bzc(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Bzc(jg-1,kg  ,lg+1)
        +     owx *  wy *  wz * Bzc(jg  ,kg-1,lg+1)
        +      wx *  wy *  wz * Bzc(jg-1,kg-1,lg+1);

    // interp from fine stagged to fine nodal
    Real bf = 0.25_rt*(Bzf(j-1,k-1,l) + Bzf(j,k-1,l) + Bzf(j-1,k,l) + Bzf(j,k,l));

#endif

    Bza(j,k,l) = bg + (bf-bc);
}

// With the FDTD Maxwell solver, this is the linear interpolation function used to
// interpolate the Bx field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bxa,
                               amrex::Array4<amrex::Real const> const& Bxf)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    Bxa(j,k,0) = 0.5_rt*(Bxf(j,k-1,0) + Bxf(j,k,0));
    amrex::ignore_unused(l);
#else
    Bxa(j,k,l) = 0.25_rt*(Bxf(j,k-1,l-1) + Bxf(j,k,l-1) + Bxf(j,k-1,l) + Bxf(j,k,l));
#endif
}

// With the PSATD Maxwell solver, this is the arbitrary-order interpolation function used to
// interpolate the Bx field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bxa,
                               amrex::Array4<amrex::Real const> const& Bxf,
                               const int noy,
                               const int noz,
                               amrex::Real const* stencil_coef_y,
                               amrex::Real const* stencil_coef_z)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(noy);
    amrex::ignore_unused(stencil_coef_y);
    amrex::Real res = 0.;
    for (int nz = 1; nz < noz/2+1; nz++) {
        res += 0.5_rt * stencil_coef_z[nz-1] * (Bxf(j, k + nz - 1, l) + Bxf(j, k - nz, l));
    }
    Bxa(j,k,l) = res;
#else
    amrex::Real res = 0.;
    for     (int nz = 1; nz < noz/2+1; nz++) {
        for (int ny = 1; ny < noy/2+1; ny++) {
            res += 0.25_rt * stencil_coef_y[ny-1] * stencil_coef_z[nz-1] *
                   (Bxf(j, k + ny - 1, l + nz - 1) + Bxf(j, k - ny, l + nz - 1)
                  + Bxf(j, k + ny - 1, l - nz    ) + Bxf(j, k - ny, l - nz    ));
        }
    }
    Bxa(j,k,l) = res;
#endif
}

// With the FDTD Maxwell solver, this is the linear interpolation function used to
// interpolate the By field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bya,
                               amrex::Array4<amrex::Real const> const& Byf)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    Bya(j,k,0) = 0.25_rt*(Byf(j,k,0) + Byf(j-1,k,0) + Byf(j,k-1,0) + Byf(j-1,k-1,0));
    amrex::ignore_unused(l);
#else
    Bya(j,k,l) = 0.25_rt*(Byf(j-1,k,l-1) + Byf(j,k,l-1) + Byf(j-1,k,l) + Byf(j,k,l));
#endif
}

// With the PSATD Maxwell solver, this is the arbitrary-order interpolation function used to
// interpolate the By field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bya,
                               amrex::Array4<amrex::Real const> const& Byf,
                               const int nox,
                               const int noz,
                               amrex::Real const* stencil_coef_x,
                               amrex::Real const* stencil_coef_z)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    amrex::Real res = 0.;
    for     (int nz = 1; nz < noz/2+1; nz++) {
        for (int nx = 1; nx < nox/2+1; nx++) {
            res += 0.25_rt * stencil_coef_x[nx-1] * stencil_coef_z[nz-1] *
                   (Byf(j + nx - 1, k + nz - 1, l) + Byf(j - nx, k + nz - 1, l)
                  + Byf(j + nx - 1, k - nz    , l) + Byf(j - nx, k - nz    , l));
        }
    }
    Bya(j,k,l) = res;
#else
    amrex::Real res = 0.;
    for     (int nz = 1; nz < noz/2+1; nz++) {
        for (int nx = 1; nx < nox/2+1; nx++) {
            res += 0.25_rt * stencil_coef_x[nx-1] * stencil_coef_z[nz-1] *
                   (Byf(j + nx - 1, k, l + nz - 1) + Byf(j - nx, k, l + nz - 1)
                  + Byf(j + nx - 1, k, l - nz    ) + Byf(j - nx, k, l - nz    ));
        }
    }
    Bya(j,k,l) = res;
#endif
}

// With the FDTD Maxwell solver, this is the linear interpolation function used to
// interpolate the Bz field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bza,
                               amrex::Array4<amrex::Real const> const& Bzf)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    Bza(j,k,0) = 0.5_rt*(Bzf(j-1,k,0) + Bzf(j,k,0));
    amrex::ignore_unused(l);
#else
    Bza(j,k,l) = 0.25_rt*(Bzf(j-1,k-1,l) + Bzf(j,k-1,l) + Bzf(j-1,k,l) + Bzf(j,k,l));
#endif
}

// With the PSATD Maxwell solver, this is the arbitrary-order interpolation function used to
// interpolate the Bz field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bza,
                               amrex::Array4<amrex::Real const> const& Bzf,
                               const int nox,
                               const int noy,
                               amrex::Real const* stencil_coef_x,
                               amrex::Real const* stencil_coef_y)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(noy);
    amrex::ignore_unused(stencil_coef_y);
    amrex::Real res = 0.;
    for (int nx = 1; nx < nox/2+1; nx++) {
        res += 0.5_rt * stencil_coef_x[nx-1] * (Bzf(j + nx - 1, k, l) + Bzf(j - nx, k, l));
    }
    Bza(j,k,l) = res;
#else
    amrex::Real res = 0.;
    for     (int ny = 1; ny < noy/2+1; ny++) {
        for (int nx = 1; nx < nox/2+1; nx++) {
            res += 0.25_rt * stencil_coef_x[nx-1] * stencil_coef_y[ny-1] *
                   (Bzf(j + nx - 1, k + ny - 1, l) + Bzf(j - nx, k + ny - 1, l)
                  + Bzf(j + nx - 1, k - ny    , l) + Bzf(j - nx, k - ny    , l));
        }
    }
    Bza(j,k,l) = res;
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Exa,
                               amrex::Array4<amrex::Real const> const& Exf,
                               amrex::Array4<amrex::Real const> const& Exc,
                               amrex::Array4<amrex::Real const> const& Exg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real owx = 1.0_rt-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real owy = 1.0_rt-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * Exg(jg  ,kg  ,0)
        +     owx *  wy * Exg(jg  ,kg+1,0)
        +      wx * owy * Exg(jg+1,kg  ,0)
        +      wx *  wy * Exg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wx = 0.5_rt-wx;  owx = 1.0_rt-wx;
    Real ec = owx * owy * Exc(jg  ,kg  ,0)
        +     owx *  wy * Exc(jg  ,kg+1,0)
        +      wx * owy * Exc(jg-1,kg  ,0)
        +      wx *  wy * Exc(jg-1,kg+1,0);

    // interp from fine staggered to fine nodal
    Real ef = 0.5_rt*(Exf(j-1,k,0) + Exf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0_rt-wz;

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * owz * Exg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Exg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Exg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Exg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Exg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Exg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Exg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Exg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wx = 0.5_rt-wx;  owx = 1.0_rt-wx;
    Real ec = owx * owy * owz * Exc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Exc(jg-1,kg  ,lg  )
        +     owx *  wy * owz * Exc(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Exc(jg-1,kg+1,lg  )
        +     owx * owy *  wz * Exc(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Exc(jg-1,kg  ,lg+1)
        +     owx *  wy *  wz * Exc(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Exc(jg-1,kg+1,lg+1);

    // interp from fine staggered to fine nodal
    Real ef = 0.5_rt*(Exf(j-1,k,l) + Exf(j,k,l));

#endif

    Exa(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eya,
                               amrex::Array4<amrex::Real const> const& Eyf,
                               amrex::Array4<amrex::Real const> const& Eyc,
                               amrex::Array4<amrex::Real const> const& Eyg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real owx = 1.0_rt-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real owy = 1.0_rt-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * Eyg(jg  ,kg  ,0)
        +     owx *  wy * Eyg(jg  ,kg+1,0)
        +      wx * owy * Eyg(jg+1,kg  ,0)
        +      wx *  wy * Eyg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal (Eyc is actually nodal)
    Real ec = owx * owy * Eyc(jg  ,kg  ,0)
        +     owx *  wy * Eyc(jg  ,kg+1,0)
        +      wx * owy * Eyc(jg+1,kg  ,0)
        +      wx *  wy * Eyc(jg+1,kg+1,0);

    // interp from fine staggered to fine nodal
    Real ef = Eyf(j,k,0);

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0_rt-wz;

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * owz * Eyg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Eyg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Eyg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Eyg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Eyg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Eyg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Eyg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Eyg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wy = 0.5_rt-wy;  owy = 1.0_rt-wy;
    Real ec = owx * owy * owz * Eyc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Eyc(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Eyc(jg  ,kg-1,lg  )
        +      wx *  wy * owz * Eyc(jg+1,kg-1,lg  )
        +     owx * owy *  wz * Eyc(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Eyc(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Eyc(jg  ,kg-1,lg+1)
        +      wx *  wy *  wz * Eyc(jg+1,kg-1,lg+1);

    // interp from fine staggered to fine nodal
    Real ef = 0.5_rt*(Eyf(j,k-1,l) + Eyf(j,k,l));

#endif

    Eya(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eza,
                               amrex::Array4<amrex::Real const> const& Ezf,
                               amrex::Array4<amrex::Real const> const& Ezc,
                               amrex::Array4<amrex::Real const> const& Ezg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real owx = 1.0_rt-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real owy = 1.0_rt-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * Ezg(jg  ,kg  ,0)
        +     owx *  wy * Ezg(jg  ,kg+1,0)
        +      wx * owy * Ezg(jg+1,kg  ,0)
        +      wx *  wy * Ezg(jg+1,kg+1,0);

    // interp from coarse stagged to fine nodal
    wy = 0.5_rt-wy;  owy = 1.0_rt-wy;
    Real ec = owx * owy * Ezc(jg  ,kg  ,0)
        +     owx *  wy * Ezc(jg  ,kg-1,0)
        +      wx * owy * Ezc(jg+1,kg  ,0)
        +      wx *  wy * Ezc(jg+1,kg-1,0);

    // interp from fine staggered to fine nodal
    Real ef = 0.5_rt*(Ezf(j,k-1,0) + Ezf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0_rt : 0.5_rt;
    Real owz = 1.0_rt-wz;

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * owz * Ezg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Ezg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Ezg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Ezg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Ezg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Ezg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Ezg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Ezg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wz = 0.5_rt-wz;  owz = 1.0_rt-wz;
    Real ec = owx * owy * owz * Ezc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Ezc(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Ezc(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Ezc(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Ezc(jg  ,kg  ,lg-1)
        +      wx * owy *  wz * Ezc(jg+1,kg  ,lg-1)
        +     owx *  wy *  wz * Ezc(jg  ,kg+1,lg-1)
        +      wx *  wy *  wz * Ezc(jg+1,kg+1,lg-1);

    // interp from fine staggered to fine nodal
    Real ef = 0.5_rt*(Ezf(j,k,l-1) + Ezf(j,k,l));

#endif

    Eza(j,k,l) = eg + (ef-ec);
}

// With the FDTD Maxwell solver, this is the linear interpolation function used to
// interpolate the Ex field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Exa,
                               amrex::Array4<amrex::Real const> const& Exf)
{
    using namespace amrex;
    Exa(j,k,l) = 0.5_rt*(Exf(j-1,k,l) + Exf(j,k,l));
}

// With the PSATD Maxwell solver, this is the arbitrary-order interpolation function used to
// interpolate the Ex field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Exa,
                               amrex::Array4<amrex::Real const> const& Exf,
                               const int nox,
                               amrex::Real const* stencil_coef_x)
{
    using namespace amrex;
    amrex::Real res = 0.;
    for (int nx = 1; nx < nox/2+1; nx++) {
        res += 0.5_rt * stencil_coef_x[nx-1] * (Exf(j + nx - 1, k, l) + Exf(j - nx, k, l));
    }
    Exa(j,k,l) = res;
}

// With the FDTD Maxwell solver, this is the linear interpolation function used to
// interpolate the Ey field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eya,
                               amrex::Array4<amrex::Real const> const& Eyf)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    Eya(j,k,0) = Eyf(j,k,0);
    amrex::ignore_unused(l);
#else
    Eya(j,k,l) = 0.5_rt*(Eyf(j,k-1,l) + Eyf(j,k,l));
#endif
}

// With the PSATD Maxwell solver, this is the arbitrary-order interpolation function used to
// interpolate the Ey field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eya,
                               amrex::Array4<amrex::Real const> const& Eyf,
                               const int noy,
                               amrex::Real const* stencil_coef_y)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(noy);
    amrex::ignore_unused(stencil_coef_y);
    Eya(j,k,l) = Eyf(j,k,l);
#else
    amrex::Real res = 0._rt;
    for (int ny = 1; ny < noy/2+1; ny++) {
        res += 0.5_rt * stencil_coef_y[ny-1] * (Eyf(j, k + ny - 1, l) + Eyf(j, k - ny, l));
    }
    Eya(j,k,l) = res;
#endif
}

// With the FDTD Maxwell solver, this is the linear interpolation function used to
// interpolate the Ez field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eza,
                               amrex::Array4<amrex::Real const> const& Ezf)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    Eza(j,k,0) = 0.5_rt*(Ezf(j,k-1,0) + Ezf(j,k,0));
    amrex::ignore_unused(l);
#else
    Eza(j,k,l) = 0.5_rt*(Ezf(j,k,l-1) + Ezf(j,k,l));
#endif
}

// With the PSATD Maxwell solver, this is the arbitrary-order interpolation function used to
// interpolate the Ez field from a staggered grid to a nodal grid, before gathering
// the field from the nodes to the particle positions (momentum-conserving gathering)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eza,
                               amrex::Array4<amrex::Real const> const& Ezf,
                               const int noz,
                               amrex::Real const* stencil_coef_z)
{
    using namespace amrex;
#if (AMREX_SPACEDIM == 2)
    amrex::Real res = 0._rt;
    for (int nz = 1; nz < noz/2+1; nz++) {
        res += 0.5_rt * stencil_coef_z[nz-1] * (Ezf(j, k + nz - 1, l) + Ezf(j, k - nz, l));
    }
    Eza(j,k,l) = res;
#else
    amrex::Real res = 0.;
    for (int nz = 1; nz < noz/2+1; nz++) {
        res += 0.5_rt * stencil_coef_z[nz-1] * (Ezf(j, k, l + nz - 1) + Ezf(j, k, l - nz));
    }
    Eza(j,k,l) = res;
#endif
}

#endif
