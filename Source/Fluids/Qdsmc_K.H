/* Copyright 2024 Marco Acciarri (Helion Energy Inc.)
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_Qdsmc_K_h
#define WARPX_Qdsmc_K_h

#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>

#include <cmath>

// check: assume nodal, cell centered ?
// check: handling different dimensions ?
// only linear interpolation for now
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void do_deposit_scalar(amrex::Array4<amrex::Real> const& scalar_field,
                    const amrex::Real xp, const amrex::Real yp, const amrex::Real zp,
                    const amrex::Real ip, const amrex::Real jp, const amrex::Real kp,
                    const amrex::XDim3 & xyzmin,
                    const amrex::XDim3 & dinv,
                    const amrex::Real valuep)
{
    // calculate di, dj, dk for particle in cell
    int i = ip;
    int j = jp;
    int k = kp;

    amrex::Real dip = (xp - xyzmin.x)*dinv.x - i;
    amrex::Real djp = (yp - xyzmin.y)*dinv.y - j;
    amrex::Real dkp = (zp - xyzmin.z)*dinv.z - k;

    amrex:: Real sx[] = {1.-dip, dip};
    amrex:: Real sy[] = {1.-djp, djp};
    amrex:: Real sz[] = {1.-dkp, dkp};

    // do triple for loop, moving over i,j,k
    for         (int ioff = 0; ioff < 2; ++ioff){
        for     (int joff = 0; joff < 2; ++joff){
            for (int koff = 0; koff < 2; ++koff){
                amrex::Gpu::Atomic::AddNoRet(&scalar_field(i+ioff,j+joff,k+koff), valuep*sx[ioff]*sy[joff]*sz[koff]);
            }
        }
    }
}

#endif //WARPX_Qdsmc_K_h
