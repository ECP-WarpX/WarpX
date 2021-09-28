/* Copyright 2019-2021 Maxence Thevenet, Weiqun Zhang, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <AMReX_Dim3.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#ifndef WARPX_FILTER_H_
#define WARPX_FILTER_H_

class Filter
{
public:
    Filter () = default;

    /* Apply stencil on MultiFab (2D/3D).
     *
     * \param dstmf Destination MultiFab
     * \param srcmf source MultiFab
     * \param[in] lev mesh refinement level
     * \param scomp first component of srcmf on which the filter is applied
     * \param dcomp first component of dstmf on which the filter is applied
     * \param ncomp Number of components on which the filter is applied.
     */
    void ApplyStencil (amrex::MultiFab& dstmf,
                       const amrex::MultiFab& srcmf, const int lev, int scomp=0,
                       int dcomp=0, int ncomp=10000);

    /** Apply stencil on a FArrayBox.
     *
     * \param dstfab Destination FArrayBox
     * \param srcfab source FArrayBox
     * \param tbx Grown box on which srcfab is defined.
     * \param scomp first component of srcfab on which the filter is applied
     * \param dcomp first component of dstfab on which the filter is applied
     * \param ncomp Number of components on which the filter is applied.
     */
    void ApplyStencil (amrex::FArrayBox& dstfab,
                       const amrex::FArrayBox& srcfab, const amrex::Box& tbx,
                       int scomp=0, int dcomp=0, int ncomp=10000);

    /** Apply Filter
     *
     * This also handles guard cells (padded to zero).
     *
     * @param tbx Grown box on which srcfab is defined
     * @param src Source
     * @param dst Destination
     * @param scomp first component of srcfab on which the filter is applied
     * @param dcomp first component of dstfab on which the filter is applied
     * @param ncomp Number of components on which the filter is applied
     */
    void DoFilter (const amrex::Box& tbx,
                   amrex::Array4<amrex::Real const> const& src,
                   amrex::Array4<amrex::Real      > const& dst,
                   int scomp, int dcomp, int ncomp);

    //! In 2D, stencil_length_each_dir = {length(stencil_x), length(stencil_z)}
    amrex::IntVect stencil_length_each_dir;

protected:
    /** Stencil along each direction.
     *  in 2D, stencil_y is not initialized. */
    amrex::Gpu::DeviceVector<amrex::Real> stencil_x, stencil_y, stencil_z;

    /** Length of each stencil.
     *  In 2D, slen = {length(stencil_x), length(stencil_z), 1} */
    amrex::Dim3 slen;

private:

};
#endif // #ifndef WARPX_FILTER_H_
