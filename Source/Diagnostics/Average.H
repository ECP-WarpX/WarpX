#ifndef WARPX_AVERAGE_H_
#define WARPX_AVERAGE_H_

#include "WarpX.H"
using namespace amrex;

namespace Average{

    Real ToCellCenter ( Array4<Real const> const& mf_in_arr,
                        const IntVect stag,
                        AMREX_D_DECL( int i, int j, int k ),
                        const int comp=0 );

    void ToCellCenter ( Box const& bx,
                        Array4<Real> const& mf_out_arr,
                        Array4<Real const> const& mf_in_arr,
                        const IntVect stag,
                        int i,
                        int j,
                        int k,
                        const int dcomp,
                        const int scomp=0,
                        const int ncomp=1 );

    void ToCellCenter ( MultiFab& mf_out,
                        const MultiFab& mf_in,
                        const int dcomp,
                        const int ngrow,
                        const int scomp=0,
                        const int ncomp=1 );
}

#endif // WARPX_AVERAGE_H_
