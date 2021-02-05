/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_amrex_qed_wrapper_commons_h_
#define WARPX_amrex_qed_wrapper_commons_h_

/**
 * This header contains some common #define directives
 * used by the QED library wrappers and related
 * components.
 */

#include <AMReX_AmrCore.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

/**
 * PICSAR uses PXRMP_GPU to decorate methods which should be
 * compiled for GPU. The user has to set it to the right value
 * (AMREX_GPU_DEVICE in this case).
 */
#define PXRMP_WITH_GPU
#define PXRMP_GPU_QUALIFIER AMREX_GPU_HOST_DEVICE
//_________________________

/**
 * PICSAR uses internally some specifiers analogous to
 * AMREX_RESTRICT and AMREX_FORCE_INLINE. These definitions
 * set the aformentioned specifiers to AMREX_RESTRICT and
 * AMREX_FORCE_INLINE.
 */
#define PXRMP_RESTRICT AMREX_RESTRICT
#define PXRMP_FORCE_INLINE AMREX_FORCE_INLINE
//_________________________

#endif //WARPX_amrex_qed_wrapper_commons_h_
