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
Average::CoarsenAndInterpolateLoop ( MultiFab& mf_dst,
                                     const MultiFab& mf_src,
                                     const int dcomp,
                                     const int scomp,
                                     const int ncomp,
                                     const IntVect crse_ratio )
{
    // Staggerings of input fine MultiFab and output coarse MultiFab
    const IntVect stag_src = mf_src.boxArray().ixType().toIntVect();
    const IntVect stag_dst = mf_dst.boxArray().ixType().toIntVect();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( mf_src.nGrowVect() >= IntVect(stag_dst-stag_src),
        "input MultiFab does not have enough guard cells for this interpolation" );

    // Auxiliary integer arrays (always 3D): staggering of source fine MultiFab (sf),
    // staggering of destination coarse MultiFab (sc), and coarsening ratio (cr)
    int sf[3], sc[3], cr[3];

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
        const Box& bx = mfi.tilebox();
        Array4<Real> const& arr_dst = mf_dst.array( mfi );
        Array4<Real const> const& arr_src = mf_src.const_array( mfi );
        ParallelFor( bx, ncomp,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k, int n )
                     {
                         arr_dst(i,j,k,n+dcomp) = Average::CoarsenAndInterpolateKernel(
                             arr_src, sf, sc, cr, i, j, k, n+scomp );
                     } );
    }
}

void
Average::CoarsenAndInterpolate ( MultiFab& mf_dst,
                                 const MultiFab& mf_src,
                                 const int dcomp,
                                 const int scomp,
                                 const int ncomp,
                                 const IntVect crse_ratio )
{
    BL_PROFILE( "Average::CoarsenAndInterpolate" );

    // Convert BoxArray of source MultiFab to staggering of destination MultiFab and coarsen it (if possible)
    BoxArray ba_tmp = amrex::convert( mf_src.boxArray(), mf_dst.ixType().toIntVect() );
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ba_tmp.coarsenable( crse_ratio ),
        "source MultiFab converted to staggering of destination MultiFab is not coarsenable" );
    ba_tmp.coarsen( crse_ratio );

    if ( ba_tmp == mf_dst.boxArray() and mf_src.DistributionMap() == mf_dst.DistributionMap() )
        Average::CoarsenAndInterpolateLoop( mf_dst, mf_src, dcomp, scomp, ncomp, crse_ratio );
    else
    {
        // Cannot coarsen into MultiFab with different BoxArray or DistributionMapping:
        // 1) create temporary MultiFab on coarsened version of source BoxArray with same DistributionMapping
        MultiFab mf_tmp( ba_tmp, mf_src.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory() );
        // 2) interpolate from mf_src to mf_tmp (start writing into component 0)
        Average::CoarsenAndInterpolateLoop( mf_tmp, mf_src, 0, scomp, ncomp, crse_ratio );
        // 3) copy from mf_tmp to mf_dst (with different BoxArray or DistributionMapping)
        mf_dst.copy( mf_tmp, 0, dcomp, ncomp );
    }
}
