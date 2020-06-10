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
 *  The `double` version of shape factor functions can be used
 *  in current deposition to ensure that current deposited by
 *  particles that move only a small distance is still resolved.
 *  Without this, single and double precision version of the code
 *  can give disagreeing results in the time evolution for some
 *  problem setups.
 */
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor(amrex::Real* const sx, amrex::Real xint)
{
    return 0;
};

template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor_dbl(double* const sx, double xint)
{
    return 0;
};

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 0
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <0> (amrex::Real* const sx, amrex::Real xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid+0.5_rt);
    sx[0] = 1.0_rt;
    return j;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor_dbl <0> (double* const sx, double xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid+0.5);
    sx[0] = 1.0;
    return j;
}

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 1
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <1> (amrex::Real* const sx, amrex::Real xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid);
    const amrex::Real xint = xmid-j;
    sx[0] = 1.0_rt - xint;
    sx[1] = xint;
    return j;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor_dbl <1> (double* const sx, double xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid);
    const double xint = xmid-j;
    sx[0] = 1.0 - xint;
    sx[1] = xint;
    return j;
}

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 2
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <2> (amrex::Real* const sx, amrex::Real xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid+0.5_rt);
    const amrex::Real xint = xmid-j;
    sx[0] = 0.5_rt*(0.5_rt-xint)*(0.5_rt-xint);
    sx[1] = 0.75_rt-xint*xint;
    sx[2] = 0.5_rt*(0.5_rt+xint)*(0.5_rt+xint);
    // index of the leftmost cell where particle deposits
    return j-1;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor_dbl <2> (double* const sx, double xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid+0.5);
    const double xint = xmid-j;
    sx[0] = 0.5*(0.5-xint)*(0.5-xint);
    sx[1] = 0.75-xint*xint;
    sx[2] = 0.5*(0.5+xint)*(0.5+xint);
    // index of the leftmost cell where particle deposits
    return j-1;
}

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 3
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <3> (amrex::Real* const sx, amrex::Real xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid);
    const amrex::Real xint = xmid-j;
    sx[0] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[1] = 2.0/3.0-xint*xint*(1.0-xint/2.0);
    sx[2] = 2.0/3.0-(1.0-xint)*(1.0-xint)*(1.0-0.5*(1.0-xint));
    sx[3] = 1.0/6.0*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return j-1;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor_dbl <3> (double* const sx, double xmid){

    using namespace amrex;

    const auto j = static_cast<int>(xmid);
    const double xint = xmid-j;
    sx[0] = 1.0_rt/6.0_rt*(1.0_rt-xint)*(1.0_rt-xint)*(1.0_rt-xint);
    sx[1] = 2.0_rt/3.0_rt-xint*xint*(1.0_rt-xint/2.0_rt);
    sx[2] = 2.0_rt/3.0_rt-(1.0_rt-xint)*(1.0_rt-xint)*(1.0_rt-0.5_rt*(1.0_rt-xint));
    sx[3] = 1.0_rt/6.0_rt*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return j-1;
}

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialized templates are defined below for orders 1, 2 and 3.
 */
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor (amrex::Real* const sx,
                                  const amrex::Real x_old,
                                  const int i_new);

template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor_dbl (double* const sx,
                                      const double x_old,
                                      const int i_new);

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 1
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <1> (amrex::Real* const sx,
                                      const amrex::Real x_old,
                                      const int i_new){

    using namespace amrex;

    const auto i = static_cast<int>(x_old);
    const int i_shift = i - i_new;
    const amrex::Real xint = x_old - i;
    sx[1+i_shift] = 1.0_rt - xint;
    sx[2+i_shift] = xint;
    return i;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor_dbl <1> (double* const sx,
                                          double x_old,
                                          const int i_new){

    using namespace amrex;

    const auto i = static_cast<int>(x_old);
    const int i_shift = i - i_new;
    const double xint = x_old - i;
    sx[1+i_shift] = 1.0 - xint;
    sx[2+i_shift] = xint;
    return i;
}

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 2
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <2> (amrex::Real* const sx,
                                      const amrex::Real x_old,
                                      const int i_new){

    using namespace amrex;

    const auto i = static_cast<int>(x_old+0.5_rt);
    const int i_shift = i - (i_new + 1);
    const amrex::Real xint = x_old - i;
    sx[1+i_shift] = 0.5_rt*(0.5_rt-xint)*(0.5_rt-xint);
    sx[2+i_shift] = 0.75_rt-xint*xint;
    sx[3+i_shift] = 0.5_rt*(0.5_rt+xint)*(0.5_rt+xint);
    // index of the leftmost cell where particle deposits
    return i-1;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor_dbl <2> (double* const sx,
                                          const double x_old,
                                          const int i_new){

    using namespace amrex;

    const auto i = static_cast<int>(x_old+0.5);
    const int i_shift = i - (i_new + 1);
    const double xint = x_old - i;
    sx[1+i_shift] = 0.5*(0.5-xint)*(0.5-xint);
    sx[2+i_shift] = 0.75-xint*xint;
    sx[3+i_shift] = 0.5*(0.5+xint)*(0.5+xint);
    // index of the leftmost cell where particle deposits
    return i-1;
}

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 3
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <3> (amrex::Real* const sx,
                                      const amrex::Real x_old,
                                      const int i_new){

    using namespace amrex;

    const auto i = static_cast<int>(x_old);
    const int i_shift = i - (i_new + 1);
    const amrex::Real xint = x_old - i;
    sx[1+i_shift] = 1.0_rt/6.0_rt*(1.0_rt-xint)*(1.0_rt-xint)*(1.0_rt-xint);
    sx[2+i_shift] = 2.0_rt/3.0_rt-xint*xint*(1.0_rt-xint/2.0_rt);
    sx[3+i_shift] = 2.0_rt/3.0_rt-(1.0_rt-xint)*(1.0_rt-xint)*(1.0_rt-0.5_rt*(1.0_rt-xint));
    sx[4+i_shift] = 1.0_rt/6.0_rt*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return i-1;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor_dbl <3> (double* const sx,
                                          const double x_old,
                                          const int i_new){

    using namespace amrex;

    const auto i = static_cast<int>(x_old);
    const int i_shift = i - (i_new + 1);
    const double xint = x_old - i;
    sx[1+i_shift] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[2+i_shift] = 2.0/3.0-xint*xint*(1.0-xint/2.0);
    sx[3+i_shift] = 2.0/3.0-(1.0-xint)*(1.0-xint)*(1.0-0.5*(1.0-xint));
    sx[4+i_shift] = 1.0/6.0*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return i-1;
}

#endif // SHAPEFACTORS_H_
