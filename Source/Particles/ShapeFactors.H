/* Copyright 2019 Maxence Thevenet, Michael Rowan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef SHAPEFACTORS_H_
#define SHAPEFACTORS_H_

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialized templates are defined below for orders 0 to 3.
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
    int operator()(T* const /*sx*/, T /*xint*/) const { return 0; }
};

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 0
 */
template <>
struct Compute_shape_factor< 0 >
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, T xmid) const
    {
        const auto j = static_cast<int>(xmid + T(0.5));
        sx[0] = T(1.0);
        return j;
    }
};

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 1
 */
template <>
struct Compute_shape_factor< 1 >
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, T xmid) const
    {
        const auto j = static_cast<int>(xmid);
        const T xint = xmid - T(j);
        sx[0] = T(1.0) - xint;
        sx[1] = xint;
        return j;
    }
};

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 2
 */
template <>
struct Compute_shape_factor< 2 >
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, T xmid) const
    {
        const auto j = static_cast<int>(xmid + T(0.5));
        const T xint = xmid - T(j);
        sx[0] = T(0.5)*(T(0.5) - xint)*(T(0.5) - xint);
        sx[1] = T(0.75) - xint*xint;
        sx[2] = T(0.5)*(T(0.5) + xint)*(T(0.5) + xint);
        // index of the leftmost cell where particle deposits
        return j-1;
    }
};

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 3
 */
template <>
struct Compute_shape_factor< 3 >
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, T xmid) const
    {
        const auto j = static_cast<int>(xmid);
        const T xint = xmid - T(j);
        sx[0] = (T(1.0))/(T(6.0))*(T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - xint);
        sx[1] = (T(2.0))/(T(3.0)) - xint*xint*(T(1.0) - xint/(T(2.0)));
        sx[2] = (T(2.0))/(T(3.0)) - (T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - T(0.5)*(T(1.0) - xint));
        sx[3] = (T(1.0))/(T(6.0))*xint*xint*xint;
        // index of the leftmost cell where particle deposits
        return j-1;
    }
};

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialized templates are defined below for orders 1, 2 and 3.
 */
template <int depos_order>
struct Compute_shifted_shape_factor
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, const T x_old, const int i_new) const;
};

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 1
 */
template <>
struct Compute_shifted_shape_factor< 1 >
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, const T x_old, const int i_new) const
    {
        const auto i = static_cast<int>(x_old);
        const int i_shift = i - i_new;
        const T xint = x_old - T(i);
        sx[1+i_shift] = T(1.0) - xint;
        sx[2+i_shift] = xint;
        return i;
    }
};

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 2
 */
template <>
struct Compute_shifted_shape_factor< 2 >
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, const T x_old, const int i_new) const
    {
        const auto i = static_cast<int>(x_old + T(0.5));
        const int i_shift = i - (i_new + 1);
        const T xint = x_old - T(i);
        sx[1+i_shift] = T(0.5)*(T(0.5) - xint)*(T(0.5) - xint);
        sx[2+i_shift] = T(0.75) - xint*xint;
        sx[3+i_shift] = T(0.5)*(T(0.5) + xint)*(T(0.5) + xint);
        // index of the leftmost cell where particle deposits
        return i - 1;
    }
};

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 3
 */
template <>
struct Compute_shifted_shape_factor< 3 >
{
    template< typename T >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator()(T* const sx, const T x_old, const int i_new) const
    {
        const auto i = static_cast<int>(x_old);
        const int i_shift = i - (i_new + 1);
        const T xint = x_old - i;
        sx[1+i_shift] = (T(1.0))/(T(6.0))*(T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - xint);
        sx[2+i_shift] = (T(2.0))/(T(3.0)) - xint*xint*(T(1.0) - xint/(T(2.0)));
        sx[3+i_shift] = (T(2.0))/(T(3.0)) - (T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - T(0.5)*(T(1.0) - xint));
        sx[4+i_shift] = (T(1.0))/(T(6.0))*xint*xint*xint;
        // index of the leftmost cell where particle deposits
        return i - 1;
    }
};

#endif // SHAPEFACTORS_H_
