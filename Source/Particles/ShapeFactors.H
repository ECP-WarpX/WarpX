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
    amrex::Real operator()(T xmid, int cell_index) const
    {
        using namespace amrex::literals;
        if constexpr (depos_order == 0){
            return (1.0);
        }
        else if constexpr (depos_order == 1){
            const T xint = xmid - static_cast<int>(xmid);
            if (cell_index==0){
                return (1.0) - xint;
            } else if (cell_index == 1){
                return xint;
            } else {
                return 0;
            }
        }
        else if constexpr (depos_order == 2){
            const T xint = xmid - static_cast<int>(xmid + (0.5_rt));
            if (cell_index==0){
                return 0.5_rt*((0.5_rt) - xint)*((0.5_rt) - xint);
            } else if (cell_index==1) {
                return (0.75_rt) - xint*xint;
            } else if (cell_index == 2){
                return (0.5_rt)*((0.5_rt) + xint)*((0.5_rt) + xint);
            } else {
                return 0;
            }
        }
        else if constexpr (depos_order == 3){
            const T xint = xmid - static_cast<int>(xmid);
            if (cell_index == 0){
                return ((1.0_rt))/((6.0_rt))*((1.0_rt) - xint)*((1.0_rt) - xint)*((1.0_rt) - xint);
            } else if (cell_index == 1) {
                return ((2.0_rt))/((3.0_rt)) - xint*xint*((1.0_rt) - xint/((2.0_rt)));
            } else if (cell_index == 2){
                return ((2.0_rt))/((3.0_rt)) - ((1.0_rt) - xint)*((1.0_rt) - xint)*((1.0_rt) - (0.5_rt)*((1.0_rt) - xint));
            } else if (cell_index == 3){
                return ((1.0_rt))/((6.0_rt))*xint*xint*xint;
            } else {
                return 0;
            }
        }

        else{
            amrex::Abort("Unknown particle shape selected in Compute_shape_factor");
            amrex::ignore_unused(xmid);
            return 0;
        }
        return 0;
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
    amrex::Real operator()(
        const T x_old, int cell_index,
        const int i_new) const
    {

        if constexpr (depos_order == 1){
            const auto i = static_cast<int>(x_old);
            const int i_shift = i - i_new;
            const T xint = x_old - T(i);
            if (cell_index == 1+i_shift){
                return T(1.0) - xint;
            } else if (cell_index == 2+i_shift){
                return xint;
            } else {
                return 0;
            }
        }

        else if constexpr (depos_order == 2){
            const auto i = static_cast<int>(x_old + T(0.5));
            const int i_shift = i - (i_new + 1);
            const T xint = x_old - T(i);
            if (cell_index ==1+i_shift){
                return T(0.5)*(T(0.5) - xint)*(T(0.5) - xint);
            } else if (cell_index == 2+i_shift){
                return T(0.75) - xint*xint;
            } else if (cell_index == 3+i_shift){
                return T(0.5)*(T(0.5) + xint)*(T(0.5) + xint);
            } else {
                return 0;
            }
        }

        else if constexpr (depos_order == 3){
            const auto i = static_cast<int>(x_old);
            const int i_shift = i - (i_new + 1);
            const T xint = x_old - i;
            if (cell_index == 1+i_shift){
                return (T(1.0))/(T(6.0))*(T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - xint);
            } else if (cell_index == 2+i_shift){
                return (T(2.0))/(T(3.0)) - xint*xint*(T(1.0) - xint/(T(2.0)));
            } else if (cell_index == 3+i_shift){
                return (T(2.0))/(T(3.0)) - (T(1.0) - xint)*(T(1.0) - xint)*(T(1.0) - T(0.5)*(T(1.0) - xint));
            } else if (cell_index == 4+i_shift){
                return (T(1.0))/(T(6.0))*xint*xint*xint;
            } else {
                return 0;
            }
            // index of the leftmost cell where particle deposits
        }
        else{
            amrex::Abort("Unknown particle shape selected in Compute_shifted_shape_factor");
            amrex::ignore_unused(i_new);
            return 0;
        }
        return 0;
    }
};

#endif // SHAPEFACTORS_H_
