/* Copyright 2023 Grant Johnson, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MusclHancock_H_
#define WARPX_MusclHancock_H_

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

// Velocity at the half step
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real V_calc (amrex::Real Q0, amrex::Real Q1, amrex::Real Q2, amrex::Real Q3,
    amrex::Real c, int comp)
{
    // Compute Ux, Uy, Uz:
    auto Ux = Q1/Q0;
    auto Uy = Q2/Q0;
    auto Uz = Q3/Q0;

    // Gamma Calc
    auto gamma = sqrt(1.0 + (Ux*Ux + Uy*Uy + Uz*Uz)/(c*c));

    // Return the proper component
    if (comp == 0) // x
        return Ux/gamma;
    else if (comp == 1) // y
        return Uy/gamma;
    else // (comp == 2) // z
        return Uz/gamma;
}

// mindmod
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real minmod (amrex::Real a, amrex::Real b)
{
    if (a > 0.0 && b > 0.0)
        return std::min(a, b);
    else if (a < 0.0 && b < 0.0)
        return std::max(a, b);
    else
        return 0.0;
}
//maxmod
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real maxmod (amrex::Real a, amrex::Real b)
{
    if (a > 0.0 && b > 0.0)
        return std::max(a, b);
    else if (a < 0.0 && b < 0.0)
        return std::min(a, b);
    else
        return 0.0;
}
// Rusanov Flux
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real flux (amrex::Real Qm, amrex::Real Qp, amrex::Real Vm, amrex::Real Vp)
{
    auto c = std::max( std::abs(Vm) , std::abs(Vp) );
    return 0.5*(Vm*Qm + Vp*Qp) - (0.5*c)*(Qp - Qm);
}
// ave_minmod
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave (amrex::Real a, amrex::Real b)
{
    if (a*b > 0.0)
        return minmod( maxmod(a,b), minmod(2.0*a,2.0*b));
    else
        return 0.0;
}

#endif /*WARPX_MusclHancock_H_*/
