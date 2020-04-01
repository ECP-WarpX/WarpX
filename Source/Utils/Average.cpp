#include "Average.H"

using namespace amrex;

void
Average::ToCellCenter ( MultiFab& mf_out,
                        const MultiFab& mf_in,
                        const int dcomp,
                        const int ngrow,
                        const int scomp,
                        const int ncomp )
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // Loop over boxes (or tiles if not on GPU)
    for (MFIter mfi( mf_out, TilingIfNotGPU() ); mfi.isValid(); ++mfi)
    {
        const Box bx = mfi.growntilebox( ngrow );
        Array4<Real> const& mf_out_arr = mf_out.array( mfi );
        Array4<Real const> const& mf_in_arr = mf_in.const_array( mfi );
        const IntVect stag = mf_in.boxArray().ixType().ixType();
        ParallelFor( bx, ncomp,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k, int n )
                     {
                         mf_out_arr(i,j,k,n+dcomp) = Average::ToCellCenter( mf_in_arr, stag, i, j, k, n+scomp );
                     } );
    }
}

void
Average::CoarsenLoop ( MultiFab& mf_cp,
                       const MultiFab& mf_fp,
                       const int dcomp,
                       const int scomp,
                       const int ncomp,
                       const IntVect ratio )
{
    const IntVect stag = mf_fp.boxArray().ixType().ixType();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // Loop over boxes (or tiles if not on GPU)
    for (MFIter mfi( mf_cp, TilingIfNotGPU() ); mfi.isValid(); ++mfi)
    {
        // NOTE: tiles defined at the coarse level
        const Box& bx = mfi.tilebox();
        Array4<Real> const& mf_cp_arr = mf_cp.array( mfi );
        Array4<Real const> const& mf_fp_arr = mf_fp.const_array( mfi );
        ParallelFor( bx, ncomp,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k, int n )
                     {
                         mf_cp_arr(i,j,k,n+dcomp) = Average::CoarsenKernel( mf_fp_arr, stag, i, j, k, n+scomp, ratio );
                     } );
    }
}

void
Average::Coarsen ( MultiFab& mf_cp,
                   const MultiFab& mf_fp,
                   const int dcomp,
                   const int scomp,
                   const int ncomp,
                   const IntVect ratio )
{
    AMREX_ASSERT( mf_cp.nComp() == mf_fp.nComp() );

    // Coarsen() fine data
    BoxArray coarsened_mf_fp_ba = mf_fp.boxArray();
    coarsened_mf_fp_ba.coarsen( ratio );

    if (coarsened_mf_fp_ba == mf_cp.boxArray() and mf_fp.DistributionMap() == mf_cp.DistributionMap())
        Average::CoarsenLoop( mf_cp, mf_fp, dcomp, scomp, ncomp, ratio );
    else
    // We copy from component scomp of the fine FArrayBox into component 0 of the
    // coarse FArrayBox because the coarse FArrayBox is a temporary FArrayBox
    // starting at component 0 and is not part of the actual coarse MultiFab mf_cp
    {
        MultiFab coarsened_mf_fp( coarsened_mf_fp_ba, mf_fp.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory() );
        Average::CoarsenLoop( coarsened_mf_fp, mf_fp, 0, scomp, ncomp, ratio );
        mf_cp.copy( coarsened_mf_fp, 0, scomp, ncomp );
    }
}
