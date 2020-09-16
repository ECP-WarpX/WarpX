/* Copyright 2020 Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARSER_WRAPPER_H_
#define WARPX_PARSER_WRAPPER_H_

#include "Parser/WarpXParser.H"
#include "Parser/GpuParser.H"

#include <AMReX_Gpu.H>

#include <memory>

/**
 * \brief
 *  The HostDeviceParser struct is non-owning and is suitable for being
 *  value captured by device lamba.
 */

template <int N>
struct HostDeviceParser
{
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
        return (*m_gpu_parser)(var...);
#endif
#else
// WarpX compiled for CPU
        return (*m_gpu_parser)(var...);
#endif
    }

#ifdef AMREX_USE_GPU
    struct wp_node * m_gpu_parser_ast = nullptr;
#endif
    GpuParser<N> const* m_gpu_parser = nullptr;
};

/**
 * \brief
 *  The ParserWrapper struct is constructed to safely use the GpuParser class
 *  that can essentially be though of as a raw pointer. The GpuParser does
 *  not have an explicit destructor and the AddPlasma subroutines handle the GpuParser
 *  in a safe way. The ParserWrapper struct is used to avoid memory leak
 *  in the EB parser functions.
 */
template <int N>
struct ParserWrapper
    : public GpuParser<N>
{
    using GpuParser<N>::GpuParser;

    ParserWrapper (ParserWrapper<N> const&) = delete;
    void operator= (ParserWrapper<N> const&) = delete;

    ~ParserWrapper() { GpuParser<N>::clear(); }

    void operator() () const = delete; // Hide GpuParser<N>::operator()

    HostDeviceParser<N> getParser () const {
#ifdef AMREX_USE_GPU
        return HostDeviceParser<N>{this->m_gpu_parser_ast, static_cast<GpuParser<N> const*>(this)};
#else
        return HostDeviceParser<N>{static_cast<GpuParser<N> const*>(this)};
#endif
    }
};

template <int N>
HostDeviceParser<N> getParser (ParserWrapper<N> const* parser_wrapper)
{
    return parser_wrapper ? parser_wrapper->getParser() : HostDeviceParser<N>{};
}

template <int N>
HostDeviceParser<N> getParser (std::unique_ptr<ParserWrapper<N> > const& parser_wrapper)
{
    return parser_wrapper ? parser_wrapper->getParser() : HostDeviceParser<N>{};
}

#endif
