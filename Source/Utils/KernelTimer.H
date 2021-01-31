/* Copyright 2019-2020 Michael Rowan, Axel Huebl, Kevin Gott
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef KERNELTIMER_H_
#define KERNELTIMER_H_

#ifdef AMREX_USE_GPU
#include <limits.h>
/**
 * \brief Defines a timer object to be used on GPU; measured summed thread time.
 */
class
KernelTimer
{
public:
    /** Constructor.
     * \param[in] do_timing Controls whether timer is active.
     * \param[in,out] cost Pointer to cost of corresponding box.
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    KernelTimer(const bool do_timing, long long int* cost) :
                m_do_timing(do_timing) {
        if (do_timing && cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            // To store present thread's time for deposition
            __shared__ long long int cost_shared;
            cost_shared = 0;
            m_cost_shared = &cost_shared;
            m_cost = cost;

            // Start the timer
            m_wt = clock64();

#elif defined(AMREX_USE_DPCPP)
            // Todo: DPC++
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( do_timing == false,
                    "KernelTimer not supported for SYCL." );
#endif
        }
    }

    //** Destructor. */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    ~KernelTimer() {
        if (m_do_timing && m_cost) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            // Thread sum
            m_wt = clock64() - m_wt;
            //amrex::Gpu::Atomic::Add( m_cost_shared, amrex::Real(m_wt));

            // Kernel time
            amrex::Gpu::Atomic::Max( m_cost_shared, m_wt);

            // Make sure everyone finished current deposit
            __syncthreads();

            // Copy tmp smem arrays back to global
            if (threadIdx.x == 0) {
                // Thread sum
                //amrex::Gpu::Atomic::Add( m_cost, *m_cost_shared);

                // Kernel time
                // Cycles
                amrex::Gpu::Atomic::Max( m_cost, *m_cost_shared);
            }
            cudaDeviceSynchronize();

            if (threadIdx.x == 0 && blockIdx.x == 0) {
                // Frequency
                int clock_rate, dev_id;
#ifdef AMREX_USE_HIP
                hipGetDevice(&dev_id);
                hipDeviceGetAttribute(&clock_rate, hipDevAttrClockRate, dev_id);
#elif defined(AMREX_USE_CUDA)
                cudaGetDevice(&dev_id); // kHz
                cudaDeviceGetAttribute(&clock_rate, cudaDevAttrClockRate, dev_id);
#endif
                // shader clock frequency/2 --> primary clock frequency
                amrex::Real time = amrex::Real(*m_cost)/amrex::Real(0.5*clock_rate); // msec
                time *= 1.0e-3;       // sec
                printf("I measured %lld cycles for the dt.            \n", *m_cost);
                printf("I measured %f sec for the current deposition. \n", time);
                printf("I measured %d kHz for the clock rate.         \n", clock_rate);
            }
#endif
        }
    }

private:
    /** Stores whether kernel timer is active. */
    bool m_do_timing;

    /** Location in which to accumulate costs from all blocks.
     *  (after threads accumulate costs in shared memory).
     */
    long long int* m_cost;

    /** Store the time difference (cost) from a single thread. */
    long long int m_wt;

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    /** Location in which which to accumulate from each thread
     *  (this will be declared shared in constructor).
     */
    long long int* m_cost_shared;
#endif
};

#endif
#endif
