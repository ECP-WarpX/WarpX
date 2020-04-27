#ifndef WARPX_COARSEN_MR_H_
#define WARPX_COARSEN_MR_H_

#include "WarpX.H"

using namespace amrex;

namespace CoarsenMR{

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Real Interp ( Array4<Real const> const& arr_src,
                  int const* const AMREX_RESTRICT sf,
                  int const* const AMREX_RESTRICT sc,
                  int const* const AMREX_RESTRICT cr,
                  const int i,
                  const int j,
                  const int k,
                  const int comp )
    {
        // This is efficiently used as a method to add neutral elements beyond guard cells
        auto const arr_src_safe = [arr_src]
            AMREX_GPU_DEVICE (int const i, int const j, int const k, int const n) noexcept
            {
                return arr_src.contains(i, j, k) ? arr_src(i, j, k, n) :  0.0_rt;
            };

        // Indices of source coarse array
        const int ic[3] = { i, j, k };

        // Number of points to loop over and starting indices of fine array in each direction
        int np[3], idx_min[3];

        // Compute number of points to loop over in each direction
        for ( int l = 0; l < 3; ++l ) {
            if ( cr[l] == 1 ) np[l] = 1; // no coarsening
            else              np[l] = cr[l]*(1-sf[l])*(1-sc[l]) // if cell-centered
                                      +2*(cr[l]-1)*sf[l]*sc[l]+1*sf[l]*sc[l]; // if nodal
        }

        // Compute starting indices of fine array in each direction
        for ( int l = 0; l < 3; ++l ) {
            if ( cr[l] == 1 ) idx_min[l] = ic[l]; // no coarsening
            else              idx_min[l] = ic[l]*cr[l]*(1-sf[l])*(1-sc[l]) // if cell-centered
                                           +(ic[l]*cr[l]-cr[l]+1)*sf[l]*sc[l]; // if nodal
        }

        // Auxiliary integer variables
        const int numx = np[0];
        const int numy = np[1];
        const int numz = np[2];
        const int imin = idx_min[0];
        const int jmin = idx_min[1];
        const int kmin = idx_min[2];
        const int sfx  = sf[0];
        const int sfy  = sf[1];
        const int sfz  = sf[2];
        const int scx  = sc[0];
        const int scy  = sc[1];
        const int scz  = sc[2];
        const int crx  = cr[0];
        const int cry  = cr[1];
        const int crz  = cr[2];
        int  ii, jj, kk;
        Real wx, wy, wz;

        // Interpolate over points computed above
        Real c = 0.0_rt;
        for (int kref = 0; kref < numz; ++kref) {
            kk = kmin+kref;
            for (int jref = 0; jref < numy; ++jref) {
                jj = jmin+jref;
                for (int iref = 0; iref < numx; ++iref) {
                    ii = imin+iref;
                    wx = (1.0_rt/static_cast<Real>(numx))*(1-sfx)*(1-scx) // if cell-centered
                         +((abs(crx-abs(ii-i*crx)))/static_cast<Real>(crx*crx))*sfx*scx; // if nodal
                    wy = (1.0_rt/static_cast<Real>(numy))*(1-sfy)*(1-scy) // if cell-centered
                         +((abs(cry-abs(jj-j*cry)))/static_cast<Real>(cry*cry))*sfy*scy; // if nodal
                    wz = (1.0_rt/static_cast<Real>(numz))*(1-sfz)*(1-scz) // if cell-centered
                         +((abs(crz-abs(kk-k*crz)))/static_cast<Real>(crz*crz))*sfz*scz; // if nodal
                    c += wx*wy*wz*arr_src_safe(ii,jj,kk,comp);
                }
            }
        }
        return c;
    };

    void Loop ( MultiFab& mf_dst,
                const MultiFab& mf_src,
                const int ncomp,
                const IntVect ngrow,
                const IntVect crse_ratio=IntVect(1) );

    void Coarsen ( MultiFab& mf_dst,
                   const MultiFab& mf_src,
                   const IntVect crse_ratio=IntVect(1) );
}

#endif // WARPX_COARSEN_MR_H_
