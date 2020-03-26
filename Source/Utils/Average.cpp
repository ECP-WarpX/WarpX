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
Average::FineToCoarse ( const MultiFab& S_fine,
                        MultiFab& S_crse,
                        int scomp,
                        int ncomp,
                        const IntVect& ratio )
{
    BL_PROFILE("amrex::average_down");
    AMREX_ASSERT(S_crse.nComp() == S_fine.nComp());
    AMREX_ASSERT((S_crse.is_cell_centered() && S_fine.is_cell_centered()) ||
                 (S_crse.is_nodal()         && S_fine.is_nodal()));

    bool is_cell_centered = S_crse.is_cell_centered();

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

    if (crse_S_fine_BA == S_crse.boxArray() and S_fine.DistributionMap() == S_crse.DistributionMap())
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& bx = mfi.tilebox();
            Array4<Real> const& crsearr = S_crse.array(mfi);
            Array4<Real const> const& finearr = S_fine.const_array(mfi);

            if (is_cell_centered) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown(tbx,crsearr,finearr,scomp,scomp,ncomp,ratio);
                });
            } else {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown_nodes(tbx,crsearr,finearr,scomp,scomp,ncomp,ratio);
                });
            }
        }
    }
    else
    {
        MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, 0, MFInfo(), FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& bx = mfi.tilebox();
            Array4<Real> const& crsearr = crse_S_fine.array(mfi);
            Array4<Real const> const& finearr = S_fine.const_array(mfi);

            //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
            //        because the crse fab is a temporary which was made starting at comp 0, it is
            //        not part of the actual crse multifab which came in.

            if (is_cell_centered) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown(tbx,crsearr,finearr,0,scomp,ncomp,ratio);
                });
            } else {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown_nodes(tbx,crsearr,finearr,0,scomp,ncomp,ratio);
                });
            }
        }

        S_crse.copy(crse_S_fine,0,scomp,ncomp);
    }
}
