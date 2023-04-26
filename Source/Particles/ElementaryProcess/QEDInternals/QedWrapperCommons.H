/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_amrex_qed_wrapper_commons_h_
#define WARPX_amrex_qed_wrapper_commons_h_
/*
 * This header contain some common #define directives
 * used by the QED library wrappers and related
 * components. It also defines a thin wrapper around
 * amrex::Gpu::DeviceVector.
 */

#include <AMReX_AmrCore.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>


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


/**
* Thin Wrapper around amrex::Gpu::DeviceVector<Real> to make it usable
* as the building block of the lookup tables when managed memory is not used.
*
* @tparam Real the floating point type to be used
*/
template <typename Real>
class AmrexDeviceVectorWrapper : public amrex::Gpu::DeviceVector<Real>
{
    using AGDV = amrex::Gpu::DeviceVector<Real>;
    public:

    /**
    * All the arguments passed to AmrexDeviceVectorWrapper constructor
    * are forwarded to amrex::Gpu::DeviceVector constructor.
    *
    * @tparam Args the constructor arguments
    */
    template<typename... Args>
    AmrexDeviceVectorWrapper(Args&&... args) : AGDV(std::forward<Args>(args)...)
    {}

    /**
    * This function may be called in some places inside picsar_table classes.
    * It forces a copy of the CPU data to the GPU
    *
    */
    void pxr_sync()
    {
        amrex::Gpu::synchronize();
    }
};

#endif //WARPX_amrex_qed_wrapper_commons_h_
