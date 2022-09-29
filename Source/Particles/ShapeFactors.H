/* Copyright 2019-2021 Maxence Thevenet, Michael Rowan, Luca Fedeli, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef SHAPEFACTORS_H_
#define SHAPEFACTORS_H_

#include <AMReX.H>
#include <AMReX_GpuQualifiers.H>


/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specializations are defined for orders 0 to 3 (using "if constexpr").
 *  Shape factor functors may be evaluated with double arguments
 *  in current deposition to ensure that current deposited by
 *  particles that move only a small distance is still resolved.
 *  Without this safeguard, single and double precision versions
 *  can give disagreeing results in the time evolution for some
 *  problem setups.
 */
template <int depos_order>
struct Compute_shape_factor
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(
        T* const sx,
        T xmid) const
    {
        if constexpr (depos_order == 0){
            const auto j = static_cast<int>(xmid + T(0.5));
            sx[0] = T(1.0);
            return j;
        }
        else if constexpr (depos_order == 1){
            const auto j = static_cast<int>(xmid);
            const T xint = xmid - T(j);
            sx[0] = T(1.0) - xint;
            sx[1] = xint;
            return j;
        }
        else if constexpr (depos_order == 2){
            const auto j = static_cast<int>(xmid + T(0.5));
            const T xint = xmid - T(j);
            sx[0] = T(0.5)*(T(0.5) - xint)*(T(0.5) - xint);
            sx[1] = T(0.75) - xint*xint;
            sx[2] = T(0.5)*(T(0.5) + xint)*(T(0.5) + xint);
            // index of the leftmost cell where particle deposits
            return j-1;
        }
        else if constexpr (depos_order == 3){
            const auto j = static_cast<int>(xmid);
            const T xint = xmid - T(j);
            sx[0] = (T(1.0))/(T(6.0))*(T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - xint);
            sx[1] = (T(2.0))/(T(3.0)) - xint*xint*(T(1.0) - xint/(T(2.0)));
            sx[2] = (T(2.0))/(T(3.0)) - (T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - T(0.5)*(T(1.0) - xint));
            sx[3] = (T(1.0))/(T(6.0))*xint*xint*xint;
            // index of the leftmost cell where particle deposits
            return j-1;
        }
        else{
            amrex::Abort("Unknown particle shape selected in Compute_shape_factor");
            amrex::ignore_unused(sx, xmid);
        }
        return 0;
    }
};

template <int depos_order>
struct Compute_shape_factor_uni
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
            int i_min,
            T* const sx,
            T xmid) const
    {
        if constexpr (depos_order == 0){
            sx[0] = T(1.0);
        }
        else if constexpr (depos_order == 1){
            const T xint = xmid - T(i_min);
            sx[0] = T(1.0) - xint;
            sx[1] = xint;
        }
        else if constexpr (depos_order == 2){
            const T xint = xmid - T(i_min);
            sx[0] = T(0.5)*(T(0.5) - xint)*(T(0.5) - xint);
            sx[1] = T(0.75) - xint*xint;
            sx[2] = T(0.5)*(T(0.5) + xint)*(T(0.5) + xint);
        }
        else if constexpr (depos_order == 3){
            const T xint = xmid - T(i_min);
            sx[0] = (T(1.0))/(T(6.0))*(T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - xint);
            sx[1] = (T(2.0))/(T(3.0)) - xint*xint*(T(1.0) - xint/(T(2.0)));
            sx[2] = (T(2.0))/(T(3.0)) - (T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - T(0.5)*(T(1.0) - xint));
            sx[3] = (T(1.0))/(T(6.0))*xint*xint*xint;
        }
        else{
            amrex::Abort("Unknown particle shape selected in Compute_shape_factor");
            amrex::ignore_unused(i_min, sx, xmid);
        }
    }
};




/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specializations are defined below for orders 1, 2 and 3 (using "if constexpr").
 */
template <int depos_order>
struct Compute_shifted_shape_factor
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(
        T* const sx,
        const T x_old,
        const int i_new) const
    {
        if constexpr (depos_order == 1){
            const auto i = static_cast<int>(x_old);
            const T xint = x_old - T(i);

            // we avoid dynamic indexing here, to enable vectorization (CPU) and
            // reduce register/lmem footprint (GPU).
            // The below is equivalent to this:
            const int i_shift = i - i_new;
            //sx[1+i_shift] = T(1.0) - xint;
            //sx[2+i_shift] = xint;
            sx[1] = T(1.0) - xint;
            sx[2] = xint;

            AMREX_ASSERT_WITH_MESSAGE(
                (i_shift == 0) || (i_shift == -1) || (i_shift == 1),
                "Compute_shifted_shape_factor: shift out of range");

            const bool downshift = i_shift == -1;
            sx[0] = downshift ? sx[1] : T(0.0);
            sx[1] = downshift ? sx[2] : sx[1];
            sx[2] = downshift ? sx[3] : sx[2];
            sx[3] = downshift ? T(0.0) : sx[3];

            const bool upshift = i_shift == 1;
            sx[3] = upshift ? sx[2] : T(0.0);
            sx[2] = upshift ? sx[1] : sx[2];
            sx[1] = upshift ? sx[0] : sx[1];
            sx[0] = upshift ? T(0.0) : sx[0];

            return i;
        }
        else if constexpr (depos_order == 2){
            const auto i = static_cast<int>(x_old + T(0.5));
            const T xint = x_old - T(i);

            // we avoid dynamic indexing here, to enable vectorization (CPU) and
            // reduce register/lmem footprint (GPU).
            // The below is equivalent to this:
            const int i_shift = i - (i_new + 1);
            //sx[1+i_shift] = T(0.5)*(T(0.5) - xint)*(T(0.5) - xint);
            //sx[2+i_shift] = T(0.75) - xint*xint;
            //sx[3+i_shift] = T(0.5)*(T(0.5) + xint)*(T(0.5) + xint);

            AMREX_ASSERT_WITH_MESSAGE(
                (i_shift == 0) || (i_shift == -1) || (i_shift == 1),
                "Compute_shifted_shape_factor: shift out of range");

            sx[1] = T(0.5)*(T(0.5) - xint)*(T(0.5) - xint);
            sx[2] = T(0.75) - xint*xint;
            sx[3] = T(0.5)*(T(0.5) + xint)*(T(0.5) + xint);

            const bool downshift = i_shift == -1;
            sx[0] = downshift ? sx[1] : T(0.0);
            sx[1] = downshift ? sx[2] : sx[1];
            sx[2] = downshift ? sx[3] : sx[2];
            sx[3] = downshift ? sx[4] : sx[3];
            sx[4] = downshift ? T(0.0) : sx[4];

            const bool upshift = i_shift == 1;
            sx[4] = upshift ? sx[3] : T(0.0);
            sx[3] = upshift ? sx[2] : sx[3];
            sx[2] = upshift ? sx[1] : sx[2];
            sx[1] = upshift ? sx[0] : sx[1];
            sx[0] = upshift ? T(0.0) : sx[0];

            // index of the leftmost cell where particle deposits
            return i - 1;
        }
        else if constexpr (depos_order == 3){
            const auto i = static_cast<int>(x_old);
            const T xint = x_old - i;

            // we avoid dynamic indexing here, to enable vectorization (CPU) and
            // reduce register/lmem footprint (GPU).
            // The below is equivalent to this:
            const int i_shift = i - (i_new + 1);
            //sx[1+i_shift] = (T(1.0))/(T(6.0))*(T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - xint);
            //sx[2+i_shift] = (T(2.0))/(T(3.0)) - xint*xint*(T(1.0) - xint/(T(2.0)));
            //sx[3+i_shift] = (T(2.0))/(T(3.0)) - (T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - T(0.5)*(T(1.0) - xint));
            //sx[4+i_shift] = (T(1.0))/(T(6.0))*xint*xint*xint;

            AMREX_ASSERT_WITH_MESSAGE(
                (i_shift == 0) || (i_shift == -1) || (i_shift == 1),
                "Compute_shifted_shape_factor: shift out of range");

            sx[1] = (T(1.0))/(T(6.0))*(T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - xint);
            sx[2] = (T(2.0))/(T(3.0)) - xint*xint*(T(1.0) - xint/(T(2.0)));
            sx[3] = (T(2.0))/(T(3.0)) - (T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - T(0.5)*(T(1.0) - xint));
            sx[4] = (T(1.0))/(T(6.0))*xint*xint*xint;

            const bool downshift = i_shift == -1;
            sx[0] = downshift ? sx[1] : T(0.0);
            sx[1] = downshift ? sx[2] : sx[1];
            sx[2] = downshift ? sx[3] : sx[2];
            sx[3] = downshift ? sx[4] : sx[3];
            sx[4] = downshift ? sx[5] : sx[4];
            sx[5] = downshift ? T(0.0) : sx[5];

            const bool upshift = i_shift == 1;
            sx[5] = upshift ? sx[4] : T(0.0);
            sx[4] = upshift ? sx[3] : sx[4];
            sx[3] = upshift ? sx[2] : sx[3];
            sx[2] = upshift ? sx[1] : sx[2];
            sx[1] = upshift ? sx[0] : sx[1];
            sx[0] = upshift ? T(0.0) : sx[0];

            // index of the leftmost cell where particle deposits
            return i - 1;
        }
        else{
            amrex::Abort("Unknown particle shape selected in Compute_shifted_shape_factor");
            amrex::ignore_unused(sx, x_old, i_new);
        }
        return 0;
    }
};

#endif // SHAPEFACTORS_H_
