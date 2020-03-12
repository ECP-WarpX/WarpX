/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Jean-Luc Vay, Ligia Diana Amorim
 * Luca Fedeli, Mathieu Lobet, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_F_H_
#define WARPX_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus

#if AMREX_SPACEDIM==2
#define WARPX_ARLIM_ANYD(x) std::array<int,3>{(x)[0], 0, (x)[1]}.data()
#else
#define WARPX_ARLIM_ANYD(x) x
#endif

#define WARPX_TO_FORTRAN_BOX(x) WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())
#define WARPX_TO_FORTRAN_ANYD(x) (x).dataPtr(), WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())
#define WARPX_TO_FORTRAN_N_ANYD(x,n) (x).dataPtr(n), WARPX_ARLIM_ANYD((x).loVect()), WARPX_ARLIM_ANYD((x).hiVect())

#endif

#if (AMREX_SPACEDIM == 3)

#define WRPX_ZERO_OUT_BNDRY              warpx_zero_out_bndry_3d
#define WRPX_BUILD_MASK                  warpx_build_mask_3d
#define WRPX_INTERPOLATE_CIC             warpx_interpolate_cic_3d
#define WRPX_INTERPOLATE_CIC_TWO_LEVELS  warpx_interpolate_cic_two_levels_3d
#define WRPX_PUSH_LEAPFROG               warpx_push_leapfrog_3d

#elif (AMREX_SPACEDIM == 2)

#define WRPX_ZERO_OUT_BNDRY              warpx_zero_out_bndry_2d
#define WRPX_BUILD_MASK                  warpx_build_mask_2d
#define WRPX_INTERPOLATE_CIC             warpx_interpolate_cic_2d
#define WRPX_INTERPOLATE_CIC_TWO_LEVELS  warpx_interpolate_cic_two_levels_2d
#define WRPX_PUSH_LEAPFROG               warpx_push_leapfrog_2d

#endif

#ifdef __cplusplus
extern "C"
{
#endif

///
/// These functions are used in electrostatic mode.
///
    void WRPX_ZERO_OUT_BNDRY(const int* lo, const int* hi,
                             amrex::Real* input_data, amrex::Real* bndry_data,
                             const int* mask);

    void WRPX_BUILD_MASK(const int* lo, const int* hi,
                         const int* tmp_mask, int* mask, const int* ncells);

    void WRPX_INTERPOLATE_CIC_TWO_LEVELS(const amrex::ParticleReal* particles, int ns, int np,
                                         amrex::Real* Ex_p, amrex::Real* Ey_p,
#if (AMREX_SPACEDIM == 3)
                                         amrex::Real* Ez_p,
#endif
                                         const amrex::Real* Ex, const amrex::Real* Ey,
#if (AMREX_SPACEDIM == 3)
                                         const amrex::Real* Ez,
#endif
                                         const int* lo, const int* hi, const amrex::Real* dx,
                                         const amrex::Real* cEx, const amrex::Real* cEy,
#if (AMREX_SPACEDIM == 3)
                                         const amrex::Real* cEz,
#endif
                                         const int* mask,
                                         const int* clo, const int* chi, const amrex::Real* cdx,
                                         const amrex::Real* plo, const int* ng, const int* lev);

    void WRPX_INTERPOLATE_CIC(const amrex::ParticleReal* particles, int ns, int np,
                              amrex::Real* Ex_p, amrex::Real* Ey_p,
#if (AMREX_SPACEDIM == 3)
                              amrex::Real* Ez_p,
#endif
                              const amrex::Real* Ex, const amrex::Real* Ey,
#if (AMREX_SPACEDIM == 3)
                              const amrex::Real* Ez,
#endif
                              const int* lo, const int* hi,
                              const amrex::Real* plo, const amrex::Real* dx,
                              const int* ng);

    void WRPX_PUSH_LEAPFROG(amrex::ParticleReal* particles, int ns, int np,
                            amrex::ParticleReal* vx_p, amrex::ParticleReal* vy_p,
#if (AMREX_SPACEDIM == 3)
                            amrex::ParticleReal* vz_p,
#endif
                            const amrex::Real* Ex_p, const amrex::Real* Ey_p,
#if (AMREX_SPACEDIM == 3)
                            const amrex::Real* Ez_p,
#endif
                            const amrex::Real* charge, const amrex::Real* mass, const amrex::Real* dt,
                            const amrex::Real* prob_lo, const amrex::Real* prob_hi);

    void WRPX_PUSH_LEAPFROG_POSITIONS(amrex::ParticleReal* particles, int ns, int np,
                                      amrex::ParticleReal* vx_p, amrex::ParticleReal* vy_p,
#if (AMREX_SPACEDIM == 3)
                                      amrex::ParticleReal* vz_p,
#endif
                                      const amrex::Real* dt, const amrex::Real* prob_lo,
                                      const amrex::Real* prob_hi);

#ifdef __cplusplus
}
#endif

#endif //WARPX_F_H_
