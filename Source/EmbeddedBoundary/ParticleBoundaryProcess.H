/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARYPROCESS_H_
#define PARTICLEBOUNDARYPROCESS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_Random.H>

namespace ParticleBoundaryProcess {

struct NoOp {
    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const PData& /*ptd*/, int /*i*/,
                     const amrex::RealVect& /*pos*/, const amrex::RealVect& /*normal*/,
                     amrex::RandomEngine const& /*engine*/) const noexcept
    {}
};

struct Absorb {
    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const PData& ptd, int i,
                     const amrex::RealVect& /*pos*/, const amrex::RealVect& /*normal*/,
                     amrex::RandomEngine const& /*engine*/) const noexcept
    {
        auto& p = ptd.m_aos[i];
        p.id() = -p.id();
    }
};
}

#endif
