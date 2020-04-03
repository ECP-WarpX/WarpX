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
Average::CoarsenAndInterpolateLoop ( MultiFab& mf_cp,
                                     const MultiFab& mf_fp,
                                     const int dcomp,
                                     const int scomp,
                                     const int ncomp,
                                     const IntVect ratio )
{
    // Staggerings of input fine MultiFab and output coarse MultiFab
    const IntVect stag_fp = mf_fp.boxArray().ixType().ixType();
    const IntVect stag_cp = mf_cp.boxArray().ixType().ixType();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( mf_fp.nGrowVect() >= IntVect(stag_cp-stag_fp),
        "input fine MultiFab does not have enough guard cells for this interpolation" );

    // Auxiliary integer arrays (always 3D unlike IntVect objects)
    int sf[3], sc[3], cr[3];

    sf[0] = stag_fp[0];
    sf[1] = stag_fp[1];
#if   (AMREX_SPACEDIM == 2)
    sf[2] = 0;
#elif (AMREX_SPACEDIM == 3)
    sf[2] = stag_fp[2];
#endif

    sc[0] = stag_cp[0];
    sc[1] = stag_cp[1];
#if   (AMREX_SPACEDIM == 2)
    sc[2] = 0;
#elif (AMREX_SPACEDIM == 3)
    sc[2] = stag_cp[2];
#endif

    cr[0] = ratio[0];
    cr[1] = ratio[1];
#if   (AMREX_SPACEDIM == 2)
    cr[2] = 1;
#elif (AMREX_SPACEDIM == 3)
    cr[2] = ratio[2];
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // Loop over boxes (or tiles if not on GPU)
    for (MFIter mfi( mf_cp, TilingIfNotGPU() ); mfi.isValid(); ++mfi)
    {
        // Tiles defined at the coarse level
        const Box& bx = mfi.tilebox();
        Array4<Real> const& mf_cp_arr = mf_cp.array( mfi );
        Array4<Real const> const& mf_fp_arr = mf_fp.const_array( mfi );
        ParallelFor( bx, ncomp,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k, int n )
                     {
                         mf_cp_arr(i,j,k,n+dcomp) = Average::CoarsenAndInterpolateKernel(
                             mf_fp_arr, sf, sc, cr, i, j, k, n+scomp );
                     } );
    }
}

void
Average::CoarsenAndInterpolate ( MultiFab& mf_cp,
                                 const MultiFab& mf_fp,
                                 const int dcomp,
                                 const int scomp,
                                 const int ncomp,
                                 const IntVect ratio )
{
    BL_PROFILE( "Average::CoarsenAndInterpolate" );

    AMREX_D_TERM( AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ratio[0] != 0 and ( ratio[0] == 1 or ratio[0]%2 == 0 ),
                      "coarsening ratio must be an integer power of 2" );,
                  AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ratio[1] != 0 and ( ratio[1] == 1 or ratio[1]%2 == 0 ),
                      "coarsening ratio must be an integer power of 2" );,
                  AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ratio[2] != 0 and ( ratio[2] == 1 or ratio[2]%2 == 0 ),
                      "coarsening ratio must be an integer power of 2" ); );

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( mf_cp.nComp() == mf_fp.nComp(),
        "input fine MultiFab and output coarse MultiFab must have the same number of components" );

    // Coarsen() fine data
    BoxArray coarsened_mf_fp_ba = mf_fp.boxArray();
    coarsened_mf_fp_ba.coarsen( ratio );

    if (coarsened_mf_fp_ba == mf_cp.boxArray() and mf_fp.DistributionMap() == mf_cp.DistributionMap())
        Average::CoarsenAndInterpolateLoop( mf_cp, mf_fp, dcomp, scomp, ncomp, ratio );
    else
    // Copy from component scomp of the fine FArrayBox into component 0 of the coarse
    // FArrayBox because the coarse FArrayBox is a temporary FArrayBox starting at
    // component 0 and is not part of the actual coarse MultiFab mf_cp
    {
        MultiFab coarsened_mf_fp( coarsened_mf_fp_ba, mf_fp.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory() );
        Average::CoarsenAndInterpolateLoop( coarsened_mf_fp, mf_fp, 0, scomp, ncomp, ratio );
        mf_cp.copy( coarsened_mf_fp, 0, scomp, ncomp );
    }
}
