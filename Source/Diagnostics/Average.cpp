#include "Average.H"

using namespace amrex;

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real Average::ToCellCenter ( Array4<Real const> const& mf_in_arr,
                             const IntVect stag,
                             AMREX_D_DECL( int i, int j, int k ),
                             const int comp)
{
    AMREX_D_TERM( int sx = stag[0];,
                  int sy = stag[1];,
                  int sz = stag[2]; );
#if ( AMREX_SPACEDIM == 2 )
    return 0.25_rt * ( mf_in_arr(i   ,j   ,0,comp)
                       + mf_in_arr(i+sx,j   ,0,comp)
                       + mf_in_arr(i   ,j+sy,0,comp)
                       + mf_in_arr(i+sx,j+sy,0,comp) );
#elif ( AMREX_SPACEDIM == 3 )
    return 0.125_rt * ( mf_in_arr(i   ,j   ,k   ,comp)
                        + mf_in_arr(i+sx,j   ,k   ,comp)
                        + mf_in_arr(i   ,j+sy,k   ,comp)
                        + mf_in_arr(i   ,j   ,k+sz,comp)
                        + mf_in_arr(i+sx,j+sy,k   ,comp)
                        + mf_in_arr(i   ,j+sy,k+sz,comp)
                        + mf_in_arr(i+sx,j   ,k+sz,comp)
                        + mf_in_arr(i+sx,j+sy,k+sz,comp) );
#endif
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void Average::ToCellCenter ( Box const& bx,
                             Array4<Real> const& mf_out_arr,
                             Array4<Real const> const& mf_in_arr,
                             const IntVect stag,
                             int i,
                             int j,
                             int k,
                             const int dcomp,
                             const int scomp,
                             const int ncomp )
{
    for (int n=0; n < ncomp; ++n) {
#if ( AMREX_SPACEDIM == 2 )
        mf_out_arr(i,j,0,n+dcomp) = Average::ToCellCenter( mf_in_arr, stag, i, j, n+scomp );
#elif ( AMREX_SPACEDIM == 3 )
        mf_out_arr(i,j,k,n+dcomp) = Average::ToCellCenter( mf_in_arr, stag, i, j, k, n+scomp );
#endif
    }
}

void
Average::ToCellCenter ( MultiFab& mf_out,
                        const MultiFab& mf_in,
                        const int dcomp,
                        const int ngrow,
                        const int scomp,
                        const int ncomp )
{
    //TODO: asserts?
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
        ParallelFor( bx,
                     [=] AMREX_GPU_DEVICE( int i, int j, int k )
                     {
                         Average::ToCellCenter( bx, mf_out_arr, mf_in_arr, stag,
                                                i, j, k, dcomp, scomp, ncomp );
                     } );
    }
}
