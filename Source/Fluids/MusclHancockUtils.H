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


// Euler push for momentum source (r-direction)
// Note: assumes U normalized by c
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real F_r (amrex::Real r, amrex::Real u_r, amrex::Real u_theta, amrex::Real u_z, amrex::Real dt)
{
    return dt*(-u_theta*u_theta/r)/sqrt(1.0 + u_r*u_r + u_theta*u_theta + u_z*u_z) + u_r;
}

// Euler push for momentum source (theta-direction)
// Note: assumes U normalized by c
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real F_theta (amrex::Real r, amrex::Real u_r, amrex::Real u_theta, amrex::Real u_z, amrex::Real dt)
{
    return dt*(u_theta*u_r/r)/sqrt(1.0 + u_r*u_r + u_theta*u_theta + u_z*u_z) + u_theta;
}

// Velocity at the half step
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real V_calc (amrex::Real Ux, amrex::Real Uy, amrex::Real Uz,
    amrex::Real c, int comp)
{
    // Gamma Calc
    amrex::Real gamma = sqrt(1.0 + (Ux*Ux + Uy*Uy + Uz*Uz)/(c*c));

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
// Min of 3 inputs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real min3 (amrex::Real a, amrex::Real b, amrex::Real c)
{
    return std::min(a, std::min(b, c) );
}
// Max of 3 inputs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real max3 (amrex::Real a, amrex::Real b, amrex::Real c)
{
    return std::max(a, std::max(b, c) );
}
// mindmod
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real minmod3 (amrex::Real a, amrex::Real b , amrex::Real c)
{
    if (a > 0.0 && b > 0.0 && c > 0.0)
        return min3(a,b,c);
    else if (a < 0.0 && b < 0.0 && c < 0.0)
        return max3(a,b,c);
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
    amrex::Real c = std::max( std::abs(Vm) , std::abs(Vp) );
    return 0.5*(Vm*Qm + Vp*Qp) - (0.5*c)*(Qp - Qm);
}
// ave_epsilon
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave_epsilon (amrex::Real a, amrex::Real b, amrex::Real dx)
{
    amrex::Real epsilon_sq = 10.0*dx*dx*dx;
    return ((b*b + epsilon_sq)*a + (a*a + epsilon_sq)*b)/(a*a + b*b + 2.0*epsilon_sq);
}
// ave_minmod high diffusivity
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave_simple (amrex::Real a, amrex::Real b)
{
    return (a+b)/2.0;
}
// ave_minmod high diffusivity
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave_high_diff (amrex::Real a, amrex::Real b)
{
    if (a*b > 0.0)
        return minmod3( (a+b)/2.0, a, b );
    else
        return 0.0;
}
// ave_minmod high diffusivity
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave_adjustable_diff (amrex::Real a, amrex::Real b)
{
    amrex::Real sigma = 2.0*0.732050807568877;
    if (a*b > 0.0)
        return minmod3( (a+b)/2.0, sigma*a, sigma*b );
    else
        return 0.0;
}
// ave_minmod Low diffusivity
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave (amrex::Real a, amrex::Real b)
{
    if (a*b > 0.0)
        return minmod3( (a+b)/2.0, 2.0*a, 2.0*b );
    else
        return 0.0;
}
// ave_superbee
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave_superbee (amrex::Real a, amrex::Real b)
{
    if (a*b > 0.0)
        return minmod( maxmod(a,b), minmod(2.0*a,2.0*b));
    else
        return 0.0;
}
// stage2 slope limiting
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ave_stage2 (amrex::Real dQ, amrex::Real a, amrex::Real b, amrex::Real c)
{
    // sigma = sqrt(3) -1
    amrex::Real sigma = 0.732050807568877;
    amrex::Real dq_min = 2.0*std::min( b - min3(a,b,c), max3(a,b,c) - b);
    return ( std::abs(dQ)/dQ ) * std::min( std::abs(dQ) , sigma*std::abs(dq_min) );
}
// Compute the momentum
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void isolate_momentum (const amrex::Array4<amrex::Real>& N_arr, const amrex::Array4<amrex::Real>& NUx_arr,
const amrex::Array4<amrex::Real>& NUy_arr, const amrex::Array4<amrex::Real>& NUz_arr, int i, int j, int k,
amrex::Real& Ux_p_or_m, amrex::Real& Uy_p_or_m, amrex::Real& Uz_p_or_m)
{
    if (N_arr(i,j,k) > 0.0) {
        Ux_p_or_m = (NUx_arr(i, j, k) / N_arr(i,j,k));
        Uy_p_or_m = (NUy_arr(i, j, k) / N_arr(i,j,k));
        Uz_p_or_m = (NUz_arr(i, j, k) / N_arr(i,j,k));
    }
}
// Compute the zero edges
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_U_edges_to_zero (const amrex::Array4<amrex::Real>& U_edge, int i, int j, int k, amrex::Box box)
{
    if ( box.contains(i,j,k) ) {
        U_edge(i,j,k,0) = 0.0;
        U_edge(i,j,k,1) = 0.0;
        U_edge(i,j,k,2) = 0.0;
        U_edge(i,j,k,3) = 0.0;
    }
}
// Compute the zero edges (minus)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_U_edges_minus (const amrex::Array4<amrex::Real>& U_edge, int i, int j, int k, amrex::Box box,
amrex::Real U_tilde0, amrex::Real U_tilde1, amrex::Real U_tilde2, amrex::Real U_tilde3,
amrex::Real dU0x, amrex::Real dU1x, amrex::Real dU2x, amrex::Real dU3x)
{
    if ( box.contains(i,j,k) ) {
        U_edge(i,j,k,0) = U_tilde0 + dU0x/2.0;
        U_edge(i,j,k,1) = U_tilde1 + dU1x/2.0;
        U_edge(i,j,k,2) = U_tilde2 + dU2x/2.0;
        U_edge(i,j,k,3) = U_tilde3 + dU3x/2.0;
    }
}
// Compute the zero edges (plus)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_U_edges_plus (const amrex::Array4<amrex::Real>& U_edge, int i, int j, int k, amrex::Box box,
amrex::Real U_tilde0, amrex::Real U_tilde1, amrex::Real U_tilde2, amrex::Real U_tilde3,
amrex::Real dU0x, amrex::Real dU1x, amrex::Real dU2x, amrex::Real dU3x)
{
    if ( box.contains(i,j,k) ) {
        U_edge(i,j,k,0) = U_tilde0 - dU0x/2.0;
        U_edge(i,j,k,1) = U_tilde1 - dU1x/2.0;
        U_edge(i,j,k,2) = U_tilde2 - dU2x/2.0;
        U_edge(i,j,k,3) = U_tilde3 - dU3x/2.0;
    }
}
// Positivity Limiter
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void positivity_limiter (const amrex::Array4<amrex::Real>& U_edge_plus,
const amrex::Array4<amrex::Real>& U_edge_minus, const amrex::Array4<amrex::Real>& N_arr,
int i, int j, int k, amrex::Box box, amrex::Real Ux, amrex::Real Uy, amrex::Real Uz,
int comp)
{

    // comp -> x, y, z -> 0, 1, 2
    int ip, jp, kp;
    //plus_index_offsets(i, j, k, ip, jp, kp, comp);

    // Find the correct offsets
#if defined(WARPX_DIM_3D)
    if (comp == 0) { //x
        ip =  i - 1; jp =  j; kp = k;
    } else if (comp == 1){ //y
        ip =  i; jp = j - 1; kp = k;
    } else if (comp == 2){ //z
        ip =  i; jp =  j; kp = k - 1;
    }
 #elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_XZ)
    if (comp == 0) { //x
        ip =  i - 1; jp =  j; kp = k;
    } else if (comp == 2){ //z
        ip =  i; jp =  j - 1; kp = k;
    }
#else
    if (comp == 2) { //z
        ip =  i - 1; jp =  j; kp = k;
    }
#endif

    // Set the edges to zero. If one edge in a cell is zero, we must self-consistently
    // set the slope to zero (hence why we have the three cases, the first is when
    // both points exist, and the second two are are edge cases)
    if (( box.contains(i,j,k) ) && ( box.contains(ip,jp,kp) )) {
        if ((U_edge_minus(i,j,k,0) < 0.0) || (U_edge_plus(ip,jp,kp,0) < 0.0)) {
            U_edge_minus(i,j,k,0) = N_arr(i,j,k);
            U_edge_minus(i,j,k,1) = Ux;
            U_edge_minus(i,j,k,2) = Uy;
            U_edge_minus(i,j,k,3) = Uz;
            U_edge_plus(ip,jp,kp,0) = N_arr(i,j,k);
            U_edge_plus(ip,jp,kp,1) = Ux;
            U_edge_plus(ip,jp,kp,2) = Uy;
            U_edge_plus(ip,jp,kp,3) = Uz;
        }
    } else if (( box.contains(i,j,k) ) && ( box.contains(ip,jp,kp) != 1)) {
        if (U_edge_minus(i,j,k,0) < 0.0) {
            U_edge_minus(i,j,k,0) = N_arr(i,j,k);
            U_edge_minus(i,j,k,1) = Ux;
            U_edge_minus(i,j,k,2) = Uy;
            U_edge_minus(i,j,k,3) = Uz;
        }
    } else if (( box.contains(i,j,k) != 1 ) && ( box.contains(ip,jp,kp) )) {
        if (U_edge_plus(ip,jp,kp,0) < 0.0){
            U_edge_plus(ip,jp,kp,0) = N_arr(i,j,k);
            U_edge_plus(ip,jp,kp,1) = Ux;
            U_edge_plus(ip,jp,kp,2) = Uy;
            U_edge_plus(ip,jp,kp,3) = Uz;
        }
    }
}

#endif /*WARPX_MusclHancock_H_*/
