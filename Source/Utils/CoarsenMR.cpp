#include "CoarsenMR.H"

using namespace amrex;

void
CoarsenMR::Loop ( MultiFab& mf_dst,
                  const MultiFab& mf_src,
                  const int ncomp,
                  const IntVect ngrow,
                  const IntVect crse_ratio )
{
    // Staggering of source fine MultiFab and destination coarse MultiFab
    const IntVect stag_src = mf_src.boxArray().ixType().toIntVect();
    const IntVect stag_dst = mf_dst.boxArray().ixType().toIntVect();

    // Auxiliary integer arrays (always 3D)
    GpuArray<int,3> sf; // staggering of source fine MultiFab
    GpuArray<int,3> sc; // staggering of destination coarse MultiFab
    GpuArray<int,3> cr; // coarsening ratio

    sf[0] = stag_src[0];
    sf[1] = stag_src[1];
#if   (AMREX_SPACEDIM == 2)
    sf[2] = 0;
#elif (AMREX_SPACEDIM == 3)
    sf[2] = stag_src[2];
#endif

    sc[0] = stag_dst[0];
    sc[1] = stag_dst[1];
#if   (AMREX_SPACEDIM == 2)
    sc[2] = 0;
#elif (AMREX_SPACEDIM == 3)
    sc[2] = stag_dst[2];
#endif

    cr[0] = crse_ratio[0];
    cr[1] = crse_ratio[1];
#if   (AMREX_SPACEDIM == 2)
    cr[2] = 1;
#elif (AMREX_SPACEDIM == 3)
    cr[2] = crse_ratio[2];
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    // Loop over boxes (or tiles if not on GPU)
    for (MFIter mfi( mf_dst, TilingIfNotGPU() ); mfi.isValid(); ++mfi)
    {
        // Tiles defined at the coarse level
        const Box& bx = mfi.growntilebox( ngrow );
        Array4<Real> const& arr_dst = mf_dst.array( mfi );
        Array4<Real const> const& arr_src = mf_src.const_array( mfi );
        ParallelFor( bx, ncomp,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k, int n )
                     {
                         arr_dst(i,j,k,n) = CoarsenMR::Interp(
                             arr_src, sf, sc, cr, i, j, k, n );
                     } );
    }
}

void
CoarsenMR::Coarsen ( MultiFab& mf_dst,
                     const MultiFab& mf_src,
                     const IntVect crse_ratio )
{
    BL_PROFILE("CoarsenMR::Coarsen()");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( mf_src.ixType() == mf_dst.ixType(),
        "source MultiFab and destination MultiFab have different IndexType" );

    // Number of guard cells to fill on coarse patch and number of components
    const IntVect ngrow = ( mf_src.nGrowVect() + 1 ) / crse_ratio;
    const int ncomp = mf_src.nComp();

    CoarsenMR::Loop( mf_dst, mf_src, ncomp, ngrow, crse_ratio );
}
