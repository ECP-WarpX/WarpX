/* Copyright 2023
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_STATIONFUNCTOR_H_
#define WARPX_STATIONFUNCTOR_H_

#include "ComputeDiagFunctor.H"
#include <AMReX_BaseFwd.H>
#include <array>

class
StationFunctor final : public ComputeDiagFunctor
{
public:
    //StationFunctor ( const amrex::MultiFab* const mf_src, amrex::Real zlocation,
    StationFunctor (const std::array<const amrex::MultiFab* const, 6> arr_mf_src, amrex::Real zlocation,
                    const int lev, amrex::IntVect crse_ratio, const int ncomp = 1);
    virtual void operator () (amrex::MultiFab & mf_dst, const int dcomp, const int i_buffer) const override;
    void PrepareFunctorData (int i_station, bool slice_in_domain, amrex::Real zlocation, amrex::Box buffer_box,
                             int k_index, const int max_box_size, const int buffer_full) override;
    std::unique_ptr<amrex::MultiFab> AllocateSlice(amrex::Vector<int>& slice_to_full_ba_map) const;
private:
    std::array<const amrex::MultiFab * const, 6> m_arr_mf_src;
    //const amrex::MultiFab * const m_mf_src = nullptr;
    int const m_lev;
    amrex::Real m_z_location;
    int m_k_index;
    int m_max_box_size;
    amrex::Box m_buffer_box;
    int m_slice_in_domain = 0;
};





#endif
