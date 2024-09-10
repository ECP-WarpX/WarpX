/* This file is part of Warpx.
 *
 * Authors: Avigdor Veksler
 * License: BSD-3-Clause-LBNL
*/
#include "JdispFunctor.H"

#include "WarpX.H"
#include "FieldSolver/Fields.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

using namespace amrex;
using namespace warpx::fields;

JdispFunctor::JdispFunctor (int dir, int lev,
        amrex::IntVect crse_ratio, bool convertRZmodes2cartesian, int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_dir(dir), m_lev(lev),
    m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{ }

void
JdispFunctor::operator() (amrex::MultiFab& mf_dst, int dcomp, const int /*i_buffer*/) const
{
    using ablastr::fields::Direction;

    auto& warpx = WarpX::GetInstance();
    auto* hybrid_pic_model = warpx.get_pointer_HybridPICModel();

    /** pointer to total simulation current (J) multifab */
    amrex::MultiFab* mf_j = warpx.m_fields.get("current_fp",Direction{m_dir},m_lev);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(hybrid_pic_model,
        "Displacement current diagnostic is only implemented for the HybridPICModel.");
    AMREX_ASSUME(hybrid_pic_model != nullptr);

     /** pointer to current calculated from Ampere's Law (Jamp) multifab */
    amrex::MultiFab* mf_curlB = hybrid_pic_model->get_pointer_current_fp_ampere(m_lev, m_dir);;

    //if (!hybrid_pic_model) {
        // To finish this implementation, we need to implement a method to
        // calculate (∇ x B).

        // Skeleton for future implementation for solvers other than HybridPIC.
        // Get curlB multifab

        // Divide curlB multifab by mu0 to get units of current
        // mf_curlB->mult(1.0/PhysConsts::mu0)
    //}

    // A Jdisp multifab is generated to hold displacement current.
    amrex::MultiFab Jdisp( mf_j->boxArray(), mf_j->DistributionMap(), 1, mf_j->nGrowVect() );
    Jdisp.setVal(0);

    // J_displacement = curl x B / mu0 - J
    amrex::MultiFab::LinComb(
        Jdisp, 1, *mf_curlB, 0,
        -1, *mf_j, 0, 0, 1, Jdisp.nGrowVect()
    );

    if (hybrid_pic_model) {
        // Subtract the interpolated j_external value from j_displacement.
        /** pointer to external currents (Jext) multifab */
        amrex::MultiFab* mf_j_external = hybrid_pic_model->get_pointer_current_fp_external(m_lev, m_dir);

        // Index type required for interpolating Jext from their respective
        // staggering (nodal) to the Jx_displacement, Jy_displacement, Jz_displacement
        // locations. The staggering of J_displacement is the same as the
        // staggering for J, so we use J_stag as the interpolation map.
        // For interp to work below, the indices of the undefined dimensions
        // must match. We set them as (1,1,1).
        amrex::GpuArray<int, 3> Jext_IndexType = {1, 1, 1};
        amrex::GpuArray<int, 3> J_IndexType = {1, 1, 1};
        amrex::IntVect Jext_stag = mf_j_external->ixType().toIntVect();
        amrex::IntVect J_stag = mf_j->ixType().toIntVect();

        // Index types for the dimensions simulated are overwritten.
        for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Jext_IndexType[idim] = Jext_stag[idim];
            J_IndexType[idim] = J_stag[idim];
        }

        // Parameters for `interp` that maps from Jext to J.
        // The "coarsening is just 1 i.e. no coarsening"
        amrex::GpuArray<int, 3> const& coarsen = {1, 1, 1};

        // Loop through the grids, and over the tiles within each grid to
        // subtract the interpolated Jext from J_displacement.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(Jdisp, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            Array4<Real> const& Jdisp_arr = Jdisp.array(mfi);
            Array4<Real const> const& Jext = mf_j_external->const_array(mfi);

            //  Loop over cells and update the Jdisp MultiFab
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k){
                // Interpolate Jext to the staggering of J
                auto const jext_interp = ablastr::coarsen::sample::Interp(Jext, Jext_IndexType, J_IndexType, coarsen, i, j, k, 0);
                Jdisp_arr(i, j, k, 0) -= jext_interp;
            });
        }
    }

    InterpolateMFForDiag(mf_dst, Jdisp, dcomp, warpx.DistributionMap(m_lev),
                         m_convertRZmodes2cartesian);
}
