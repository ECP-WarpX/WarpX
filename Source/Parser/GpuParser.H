/* Copyright 2019-2020 Maxence Thevenet, Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_GPU_PARSER_H_
#define WARPX_GPU_PARSER_H_

#include "Parser/WarpXParser.H"

#include <AMReX_Gpu.H>
#include <AMReX_Array.H>
#include <AMReX_TypeTraits.H>
#include <AMReX.H>

// When compiled for CPU, wrap WarpXParser and enable threading.
// When compiled for GPU, store one copy of the parser in
// device memory for __device__ code, and one copy of the parser
// in host memory for __host__ code. This way, the parser can be
// efficiently called from both host and device.
template <int N>
class GpuParser
{
public:
    GpuParser (WarpXParser const& wp);

    GpuParser (GpuParser<N> const&) = delete;
    GpuParser (GpuParser<N> &&) = delete;
    void operator= (GpuParser<N> const&) = delete;
    void operator= (GpuParser<N> &&) = delete;

    void clear ();

    template <typename... Ts>
    AMREX_GPU_HOST_DEVICE
    std::enable_if_t<sizeof...(Ts) == N
                     and amrex::Same<amrex::Real,Ts...>::value,
                     amrex::Real>
    operator() (Ts... var) const noexcept
    {
#ifdef AMREX_USE_GPU
#if AMREX_DEVICE_COMPILE
// WarpX compiled for GPU, function compiled for __device__
        amrex::GpuArray<amrex::Real,N> l_var{var...};
        return wp_ast_eval<0>(m_gpu_parser_ast, l_var.data());
#else
// WarpX compiled for GPU, function compiled for __host__
        amrex::ignore_unused(var...);
        return wp_ast_eval<0>(m_cpu_parser->ast, nullptr);
#endif

#else
// WarpX compiled for CPU
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        m_var[tid] = amrex::GpuArray<amrex::Real,N>{var...};
        return wp_ast_eval<0>(m_parser[tid]->ast, nullptr);
#endif
    }

    void init_gpu_parser (WarpXParser const& wp); // public for CUDA

protected:

#ifdef AMREX_USE_GPU
    // Copy of the parser running on __device__
    struct wp_node* m_gpu_parser_ast;
    // Copy of the parser running on __host__
    struct wp_parser* m_cpu_parser;
    mutable amrex::GpuArray<amrex::Real,N> m_var;
#else
    // Only one parser
    struct wp_parser** m_parser;
    mutable amrex::GpuArray<amrex::Real,N>* m_var;
    int nthreads;
#endif
};

template <int N>
GpuParser<N>::GpuParser (WarpXParser const& wp)
{
    AMREX_ALWAYS_ASSERT(wp.depth() <= WARPX_PARSER_DEPTH);

#ifdef AMREX_USE_GPU

    struct wp_parser* a_wp = wp.m_parser;

    // Initialize CPU parser:
    m_cpu_parser = wp_parser_dup(a_wp);
    for (int i = 0; i < N; ++i) {
        wp_parser_regvar(m_cpu_parser, wp.m_varnames[i].c_str(), &m_var[i]);
    }

    // Initialize GPU parser
    init_gpu_parser(wp);

#else // not defined AMREX_USE_GPU

#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else // _OPENMP
    nthreads = 1;
#endif // _OPENMP

    m_parser = ::new struct wp_parser*[nthreads];
    m_var = ::new amrex::GpuArray<amrex::Real,N>[nthreads];

    for (int tid = 0; tid < nthreads; ++tid)
    {
#ifdef _OPENMP
        m_parser[tid] = wp_parser_dup(wp.m_parser[tid]);
        for (int i = 0; i < N; ++i) {
            wp_parser_regvar(m_parser[tid], wp.m_varnames[tid][i].c_str(), &(m_var[tid][i]));
        }
#else // _OPENMP
        m_parser[tid] = wp_parser_dup(wp.m_parser);
        for (int i = 0; i < N; ++i) {
            wp_parser_regvar(m_parser[tid], wp.m_varnames[i].c_str(), &(m_var[tid][i]));
        }
#endif // _OPENMP
    }

#endif // AMREX_USE_GPU
}

template <int N>
void GpuParser<N>::init_gpu_parser (WarpXParser const& wp)
{
#ifdef AMREX_USE_GPU

    // We create a temporary Parser on CPU for memcpy.  We cannot use
    // m_cpu_parser for this because the variables in m_cpu_parser are
    // registered for CPU use.
    struct wp_parser* cpu_tmp = wp_parser_dup(m_cpu_parser);
    for (int i = 0; i < N; ++i) {
        wp_parser_regvar_gpu(cpu_tmp, wp.m_varnames[i].c_str(), i);
    }

    m_gpu_parser_ast = (struct wp_node*)
        amrex::The_Arena()->alloc(cpu_tmp->sz_mempool);
    amrex::Gpu::htod_memcpy_async(m_gpu_parser_ast, cpu_tmp->ast, cpu_tmp->sz_mempool);

    auto dp = m_gpu_parser_ast;
    char* droot = (char*)dp;
    char* croot = (char*)(cpu_tmp->ast);
    amrex::single_task([=] AMREX_GPU_DEVICE () noexcept
    {
        wp_ast_update_device_ptr<0>(dp, droot, croot);
    });

    amrex::Gpu::synchronize();

    wp_parser_delete(cpu_tmp);
#endif
    amrex::ignore_unused(wp);
}

template <int N>
void
GpuParser<N>::clear ()
{
#ifdef AMREX_USE_GPU
    amrex::The_Arena()->free(m_gpu_parser_ast);
    wp_parser_delete(m_cpu_parser);
#else
    for (int tid = 0; tid < nthreads; ++tid)
    {
        wp_parser_delete(m_parser[tid]);
    }
    ::delete[] m_parser;
    ::delete[] m_var;
#endif
}

#endif
