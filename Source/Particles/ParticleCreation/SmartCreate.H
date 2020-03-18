/* Copyright 2019-2020 Luca Fedeli, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef SMART_CREATE_H_
#define SMART_CREATE_H_

#include "DefaultInitialization.H"

#include <AMReX_AmrCore.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>

struct SmartCreate
{
    using amrex;

    int m_num_create_real;
    const int* m_comps_r;

    int m_num_create_int;
    const int* m_m_comps_i;

    const InitializationPolicy* m_policy_real;
    const InitializationPolicy* m_policy_int;

    template <typename PartData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (
        PartData& prt, const int i_prt,
        const amrex::Real x = 0_rt,
#if (AMREX_SPACEDIM == 3)
        const amrex::Real y = 0_rt,
#endif
        const amrex::Real z = 0_rt) const noexcept
    {
        prt.m_aos[i_dst].pos(0) = x;
#if (AMREX_SPACEDIM == 3)
        prt.m_aos[i_dst].pos(1) = y;
#endif
        prt.m_aos[i_dst].pos(2) = z;

         // initialize the real components
         for (int j = 0; j < DstData::NAR; ++j)
             dst.m_rdata[j][i_dst] = initializeRealValue(m_policy_real[j]);
         for (int j = 0; j < dst.m_num_runtime_real; ++j)
             dst.m_runtime_rdata[j][i_dst] = initializeRealValue(m_policy_real[j+DstData::NAR]);

         // initialize the int components
         for (int j = 0; j < DstData::NAI; ++j)
             dst.m_idata[j][i_dst] = initializeIntValue(m_policy_int[j]);
         for (int j = 0; j < dst.m_num_runtime_int; ++j)
             dst.m_runtime_idata[j][i_dst] = initializeIntValue(m_policy_int[j+DstData::NAI]);
     }
 };

#endif //SMART_CREATE_H_
