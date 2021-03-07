#ifndef WARPX_PEC_KERNELS_H_
#define WARPX_PEC_KERNELS_H_

#include "WarpX.H"
#include "Utils/WarpXAlgorithmSelection.H"


namespace PEC {


    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int is_lo_PEC (int dir) {
        return int( WarpX::field_boundary_lo[dir] == FieldBoundaryType::PEC);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int is_hi_PEC (int dir) {
        return int( WarpX::field_boundary_hi[dir] == FieldBoundaryType::PEC);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void ZeroTangentialEfield ( const int icomp, const amrex::IntVect & dom_lo,
                                const amrex::IntVect & dom_hi,
                                const amrex::IntVect & ijk_vec,
                                amrex::Real& Efield, const amrex::IntVect& stag )
    {
        int setfield = 0;
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            // index of dimension
            const int dir = ( ( (icomp+idim)<AMREX_SPACEDIM)? (icomp+idim) : ( (icomp+idim) - AMREX_SPACEDIM));
            // loop over sides , iside = 0 (lo), iside = 1 (hi)
            for (int iside = 0; iside < 2; ++iside) {
                const int ibndry = dir + iside*AMREX_SPACEDIM;
                const int bndry_cell = ( (ibndry < AMREX_SPACEDIM) ?
                    (ijk_vec[dir] == dom_lo[dir]) : (ijk_vec[dir] == dom_hi[dir] + stag[dir]));
                const int is_tangent_to_PEC =
                    ( ( ibndry < AMREX_SPACEDIM) ? is_lo_PEC(dir) : is_hi_PEC(dir) );
                if ( bndry_cell == 1 && is_tangent_to_PEC==1) setfield = 1;
            }
        }
        if (setfield == 1) Efield = 0.;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void ZeroNormalBfield (const int icomp, const amrex::IntVect & dom_lo,
                           const amrex::IntVect & dom_hi,
                           const amrex::IntVect & ijk_vec,
                           amrex::Real& Bfield, const amrex::IntVect & stag)
    {
        int setfield = 0;
        const int dir = icomp;
        for (int iside = 0; iside < 2; ++iside) {
            const int ibndry = dir + iside * AMREX_SPACEDIM;
            const int bndry_cell = ( ibndry < AMREX_SPACEDIM ) ?
                ( ijk_vec[dir] == dom_lo[dir]) : ( ijk_vec[dir] == dom_hi[dir]+stag[dir]);
            const int is_normal_to_PEC = (
                ( ibndry < AMREX_SPACEDIM) ? is_lo_PEC(dir) : is_hi_PEC(dir) );
            if (bndry_cell == 1 && is_normal_to_PEC == 1 ) setfield = 1;
        }
        if (setfield == 1) Bfield = 0.;
    }


    bool isAnyBoundaryPEC();
    void ApplyPECtoEfield ( std::array<std::unique_ptr<amrex::MultiFab>, 3>& Efield,
                            const int lev, PatchType patch_type);
    void ApplyPECtoBfield ( std::array<std::unique_ptr<amrex::MultiFab>, 3>& Bfield,
                            const int lev, PatchType patch_type);
}

#endif // WarpX_PEC_KERNELS_H_
