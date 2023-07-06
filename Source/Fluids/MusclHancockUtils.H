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
// halfstep Q
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Q_step(int i, int j, int k,
    amrex::Real cx_half, amrex::Real cy_half, amrex::Real cz_half,
    amrex::Real Ux, amrex::Real Uy, amrex::Real Uz, amrex::Real a,
    amrex::Real clight, amrex::Real gamma, amrex::Array4<amrex::Real> const &N_arr, 
    amrex::Array4<amrex::Real> const &NUx_arr, amrex::Array4<amrex::Real> const &NUy_arr,
    amrex::Array4<amrex::Real> const &NUz_arr, 
    amrex::Real Q_minus_x[], amrex::Real Q_minus_y[], amrex::Real Q_minus_z[],
    amrex::Real Q_plus_x[], amrex::Real Q_plus_y[], amrex::Real Q_plus_z[])
{
    // Compute the Flux-Jacobian Elements in x
    auto A00x = (Ux*(Uz*Uz)+Ux*(Uy*Uy)+(Ux*Ux*Ux))/a;
    auto A01x = ((clight*clight)+(Uz*Uz)+(Uy*Uy))/a;
    auto A02x = -(Ux*Uy)/a;
    auto A03x = -(Ux*Uz)/a;

    auto A10x = -(Ux*Ux)/(gamma*gamma*gamma);
    auto A11x = (2.0*Ux*(clight*clight)+2.0*Ux*(Uz*Uz)+2.0*Ux*(Uy*Uy)+(Ux*Ux*Ux))/a;
    auto A12x = -((Ux*Ux)*Uy)/a;
    auto A13x = -((Ux*Ux)*Uz)/a;

    auto A20x = -(Ux*Uy)/(gamma*gamma*gamma);
    auto A21x = (Uy*(clight*clight)+Uy*(Uz*Uz)+(Uy*Uy*Uy))/a;
    auto A22x = (Ux*(clight*clight)+Ux*(Uz*Uz)+(Ux*Ux*Ux))/a;
    auto A23x = -(Ux*Uy*Uz)/a;

    auto A30x = -(Ux*Uz)/(gamma*gamma*gamma);
    auto A31x = (Uz*(clight*clight)+(Uz*Uz*Uz)+(Uy*Uy)*Uz)/a;
    auto A32x = -(Ux*Uy*Uz)/a;
    auto A33x = (Ux*(clight*clight)+Ux*(Uy*Uy)+(Ux*Ux*Ux))/a;


    // Compute the Flux-Jacobian Elements in y
    auto A00y = (Uy*(Uz*Uz)+(Uy*Uy*Uy)+(Ux*Ux)*Uy)/a;
    auto A01y = -(Ux*Uy)/a;
    auto A02y = ((clight*clight)+(Uz*Uz)+(Ux*Ux))/a;
    auto A03y = -(Uy*Uz)/a;

    auto A10y = -(Ux*Uy)/(gamma*gamma*gamma);
    auto A11y = (Uy*(clight*clight)+Uy*(Uz*Uz)+(Uy*Uy*Uy))/a;
    auto A12y = (Ux*(clight*clight)+Ux*(Uz*Uz)+(Ux*Ux*Ux))/a;
    auto A13y = -(Ux*Uy*Uz)/a;

    auto A20y = -(Uy*Uy)/(gamma*gamma*gamma);
    auto A21y = -(Ux*(Uy*Uy))/a;
    auto A22y = (2.0*Uy*(clight*clight)+2.0*Uy*(Uz*Uz)+(Uy*Uy*Uy)+2.0*(Ux*Ux)*Uy)/a;
    auto A23y = -((Uy*Uy)*Uz)/a;

    auto A30y = -(Uy*Uz)/(gamma*gamma*gamma);
    auto A31y = -(Ux*Uy*Uz)/a;
    auto A32y = (Uz*(clight*clight)+(Uz*Uz*Uz)+(Ux*Ux)*Uz)/a;
    auto A33y = (Uy*(clight*clight)+(Uy*Uy*Uy)+(Ux*Ux)*Uy)/a;


    // Compute the Flux-Jacobian Elements in z
    auto A00z = ((Uz*Uz*Uz)+((Uy*Uy)+(Ux*Ux))*Uz)/a;
    auto A01z = -(Ux*Uz)/a;
    auto A02z = -(Uy*Uz)/a;
    auto A03z = ((clight*clight)+(Uy*Uy)+(Ux*Ux))/a;

    auto A10z = -(Ux*Uz)/(gamma*gamma*gamma);
    auto A11z = (Uz*(clight*clight)+(Uz*Uz*Uz)+(Uy*Uy)*Uz)/a;
    auto A12z = -(Ux*Uy*Uz)/a;
    auto A13z = (Ux*(clight*clight)+Ux*(Uy*Uy)+(Ux*Ux*Ux))/a;

    auto A20z = -(Uy*Uz)/(gamma*gamma*gamma);
    auto A21z = -(Ux*Uy*Uz)/a;
    auto A22z = (Uz*(clight*clight)+(Uz*Uz*Uz)+(Ux*Ux)*Uz)/a;
    auto A23z = (Uy*(clight*clight)+(Uy*Uy*Uy)+(Ux*Ux)*Uy)/a;

    auto A30z = -(Uz*Uz)/(gamma*gamma*gamma);
    auto A31z = -(Ux*(Uz*Uz))/a;
    auto A32z = -(Uy*(Uz*Uz))/a;
    auto A33z = (2.0*Uz*(clight*clight)+(Uz*Uz*Uz)+(2.0*(Uy*Uy)+2.0*(Ux*Ux))*Uz)/a;

    // Compute the cell slopes x
    auto dQ0x = ave( N_arr(i,j,k) - N_arr(i-1,j,k) , N_arr(i+1,j,k) - N_arr(i,j,k) );
    auto dQ1x = ave( NUx_arr(i,j,k) - NUx_arr(i-1,j,k) , NUx_arr(i+1,j,k) - NUx_arr(i,j,k) );
    auto dQ2x = ave( NUy_arr(i,j,k) - NUy_arr(i-1,j,k) , NUy_arr(i+1,j,k) - NUy_arr(i,j,k) );
    auto dQ3x = ave( NUz_arr(i,j,k) - NUz_arr(i-1,j,k) , NUz_arr(i+1,j,k) - NUz_arr(i,j,k) );

    // Compute the cell slopes y
    auto dQ0y = ave( N_arr(i,j,k) - N_arr(i,j-1,k) , N_arr(i,j+1,k) - N_arr(i,j,k) );
    auto dQ1y = ave( NUx_arr(i,j,k) - NUx_arr(i,j-1,k) , NUx_arr(i,j+1,k) - NUx_arr(i,j,k) );
    auto dQ2y = ave( NUy_arr(i,j,k) - NUy_arr(i,j-1,k) , NUy_arr(i,j+1,k) - NUy_arr(i,j,k) );
    auto dQ3y = ave( NUz_arr(i,j,k) - NUz_arr(i,j-1,k) , NUz_arr(i,j+1,k) - NUz_arr(i,j,k) );

    // Compute the cell slopes z
    auto dQ0z = ave( N_arr(i,j,k) - N_arr(i,j,k-1) , N_arr(i,j,k+1) - N_arr(i,j,k) );
    auto dQ1z = ave( NUx_arr(i,j,k) - NUx_arr(i,j,k-1) , NUx_arr(i,j,k+1) - NUx_arr(i,j,k) );
    auto dQ2z = ave( NUy_arr(i,j,k) - NUy_arr(i,j,k-1) , NUy_arr(i,j,k+1) - NUy_arr(i,j,k) );
    auto dQ3z = ave( NUz_arr(i,j,k) - NUz_arr(i,j,k-1) , NUz_arr(i,j,k+1) - NUz_arr(i,j,k) );

    // Compute Q ([ N, NU]) at the halfsteps (Q_tidle) using the slopes (dQ)
    auto AdQ0x = A00x*dQ0x + A01x*dQ1x + A02x*dQ2x + A03x*dQ3x;
    auto AdQ1x = A10x*dQ0x + A11x*dQ1x + A12x*dQ2x + A13x*dQ3x;
    auto AdQ2x = A20x*dQ0x + A21x*dQ1x + A22x*dQ2x + A23x*dQ3x;
    auto AdQ3x = A30x*dQ0x + A31x*dQ1x + A32x*dQ2x + A33x*dQ3x;
    auto AdQ0y = A00y*dQ0y + A01y*dQ1y + A02y*dQ2y + A03y*dQ3y;
    auto AdQ1y = A10y*dQ0y + A11y*dQ1y + A12y*dQ2y + A13y*dQ3y;
    auto AdQ2y = A20y*dQ0y + A21y*dQ1y + A22y*dQ2y + A23y*dQ3y;
    auto AdQ3y = A30y*dQ0y + A31y*dQ1y + A32y*dQ2y + A33y*dQ3y;
    auto AdQ0z = A00z*dQ0z + A01z*dQ1z + A02z*dQ2z + A03z*dQ3z;
    auto AdQ1z = A10z*dQ0z + A11z*dQ1z + A12z*dQ2z + A13z*dQ3z;
    auto AdQ2z = A20z*dQ0z + A21z*dQ1z + A22z*dQ2z + A23z*dQ3z;
    auto AdQ3z = A30z*dQ0z + A31z*dQ1z + A32z*dQ2z + A33z*dQ3z;
    auto Q_tilde0 = N_arr(i,j,k)   - cx_half*AdQ0x - cy_half*AdQ0y - cz_half*AdQ0z;
    auto Q_tilde1 = NUx_arr(i,j,k) - cx_half*AdQ1x - cy_half*AdQ1y - cz_half*AdQ1z;
    auto Q_tilde2 = NUy_arr(i,j,k) - cx_half*AdQ2x - cy_half*AdQ2y - cz_half*AdQ2z;
    auto Q_tilde3 = NUz_arr(i,j,k) - cx_half*AdQ3x - cy_half*AdQ3y - cz_half*AdQ3z;

    // Predict Q at the cell edges (x)
    Q_minus_x[0] = Q_tilde0 + dQ0x/2.0;
    Q_minus_x[1] = Q_tilde1 + dQ1x/2.0;
    Q_minus_x[2] = Q_tilde2 + dQ2x/2.0;
    Q_minus_x[3] = Q_tilde3 + dQ3x/2.0;
    Q_plus_x[0] = Q_tilde0 - dQ0x/2.0;
    Q_plus_x[1] = Q_tilde1 - dQ1x/2.0;
    Q_plus_x[2] = Q_tilde2 - dQ2x/2.0;
    Q_plus_x[3] = Q_tilde3 - dQ3x/2.0;

    // Predict Q at the cell edges (y)
    Q_minus_y[0] = Q_tilde0 + dQ0y/2.0;
    Q_minus_y[1] = Q_tilde1 + dQ1y/2.0;
    Q_minus_y[2] = Q_tilde2 + dQ2y/2.0;
    Q_minus_y[3] = Q_tilde3 + dQ3y/2.0;
    Q_plus_y[0] = Q_tilde0 - dQ0y/2.0;
    Q_plus_y[1] = Q_tilde1 - dQ1y/2.0;
    Q_plus_y[2] = Q_tilde2 - dQ2y/2.0;
    Q_plus_y[3] = Q_tilde3 - dQ3y/2.0;

    // Predict Q at the cell edges (z)
    Q_minus_z[0] = Q_tilde0 + dQ0z/2.0;
    Q_minus_z[1] = Q_tilde1 + dQ1z/2.0;
    Q_minus_z[2] = Q_tilde2 + dQ2z/2.0;
    Q_minus_z[3] = Q_tilde3 + dQ3z/2.0;
    Q_plus_z[0] = Q_tilde0 - dQ0z/2.0;
    Q_plus_z[1] = Q_tilde1 - dQ1z/2.0;
    Q_plus_z[2] = Q_tilde2 - dQ2z/2.0;
    Q_plus_z[3] = Q_tilde3 - dQ3z/2.0;
}


#endif /*WARPX_MusclHancock_H_*/
