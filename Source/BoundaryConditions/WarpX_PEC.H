#ifndef WARPX_PEC_KERNELS_H_
#define WARPX_PEC_KERNELS_H_

#include "WarpX.H"
#include "Utils/WarpXAlgorithmSelection.H"


namespace PEC {
using namespace amrex;
    /**
     * \brief Determines if the field boundary condition stored in fboundary is PEC
     *        in direction, dir, is PEC.
     *
     * \param[in] fboundary  Value containing boundary type
     * \param[in] dir        direction
     *
     * \param[out] 1 if the boundary type is PEC else 0
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int is_boundary_PEC (amrex::GpuArray<int, 3> const& fboundary, int dir) {
        return int( fboundary[dir] == FieldBoundaryType::PEC);
    }

    /**
     * \brief Sets the electric field value tangential to the PEC boundary to zero. The
     *        field value of the guard cells outside the domain boundary are set equal
     *        and opposite to the field in the valid cells at their mirrored locations.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *        For corner cells with mixed boundaries, the mirror location could be outside
     *        valid region, while still ensuring PEC condition is maintained across the
     *        PEC boundary.
     *
     * \param[in] icomp        component of the Efield being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] Efield       field data to be updated if (ijk) is at the boundary or a guard cell
     * \param[in] stag         staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetTangentialEfield (const int icomp, const amrex::IntVect & dom_lo,
                                const amrex::IntVect &dom_hi,
                                const amrex::IntVect &ijk_vec, const int n,
                                amrex::Array4<amrex::Real> const& Efield,
                                const amrex::IntVect& stag,
                                amrex::GpuArray<int, 3> const& fbndry_lo,
                                amrex::GpuArray<int, 3> const& fbndry_hi )
    {
        // Fields in guard cells set equal and opposite to cells in the mirror locations
        // across the PEC boundary. Here we just initialize it.
        amrex::IntVect ijk_mirror = ijk_vec;
        int OnPECBoundary = 0;
        int GuardCell = 0;
        amrex::Real sign = 1._rt; 
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
        const int sdim = ( (icomp == 1) ? 0 : 1 );
#else
        const int sdim = 1;               
#endif
        // Loop over dimensions that are tangential to the Efield component, icomp.
        // For 2D : for icomp==1, (Ey in XZ, Etheta in RZ), both X-Z and R-Z dimensions
        //          are tangential. Therefore, we set sdim to 0 to loop over both
        //          dimensions for this component.
        // Otherwise : sdim = 1, to ensure we only loop over tangential dimensions.
        for (int idim = sdim; idim < AMREX_SPACEDIM; ++idim) {
            // dir is index of direction that is tangential to the component, icomp.
            const int dir = ( ( (icomp+idim)<AMREX_SPACEDIM)
                          ? (icomp+idim) : ( (icomp+idim) - AMREX_SPACEDIM));
            // loop over sides, iside = 0 (lo), iside = 1 (hi)
            for (int iside = 0; iside < 2; ++iside) {
                const int is_tangent_to_PEC = ( ( iside == 0 )
                                            ? is_boundary_PEC(fbndry_lo, dir)
                                            : is_boundary_PEC(fbndry_hi, dir) );
                if (is_tangent_to_PEC == 1) {
                    // guard cell outside the domain by "ig" number cells, in direction, dir
                    const int ig = ( (iside == 0)
                                 ? (dom_lo[dir] - ijk_vec[dir])
                                 : (ijk_vec[dir] - (dom_hi[dir] + stag[dir]) ) );

                    if ( ig == 0) OnPECBoundary = 1;
                    if ( ig > 0 ) {
                        GuardCell = 1;
                        // Mirror location inside the domain by "ig" number of cells
                        // across PEC boundary in direction, dir, and side, iside.
                        ijk_mirror[dir] = ( ( iside == 0)
                                        ? (dom_lo[dir] + ig)
                                        : (dom_hi[dir] + stag[dir] - ig));
                        // The tangential Efield in guard cells is equal and opposite to the
                        // value in the respective mirror location across the PEC boundary.
                        // For corner cases, this sign change will occur for each valid
                        // PEC boundary, and is stored in the variable "sign".
                        sign *= -1._rt;
                    }
                }
            }
        }
        if (GuardCell == 1) Efield(ijk_vec,n) = sign * Efield(ijk_mirror,n);
        // if ijk_vec is on a PEC boundary in any direction, set Etangential to 0.
        if (OnPECBoundary == 1) Efield(ijk_vec,n) = 0._rt;
    }

    /**
     * \brief Sets the magnetic field value normal to the PEC boundary to zero. The
     *        field value of the guard cells outside the domain boundary are set equal
     *        and opposite to the field in the valid cells at their mirrored locations.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *
     * \param[in] icomp        component of the Bfield being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] Bfield       field data to be updated if (ijk) is at the boundary
                               or a guard cell
     * \param[in] stag         staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetNormalBfield (const int icomp, const amrex::IntVect & dom_lo,
                           const amrex::IntVect & dom_hi,
                           const amrex::IntVect & ijk_vec, const int n,
                           amrex::Array4<amrex::Real> const& Bfield,
                           const amrex::IntVect & stag,
                           amrex::GpuArray<int, 3> const& fbndry_lo,
                           amrex::GpuArray<int, 3> const& fbndry_hi )
    {
        amrex::IntVect ijk_mirror = ijk_vec;
        int OnPECBoundary = 0;
        int GuardCell = 0;
        amrex::Real sign = 1._rt;
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
        if (icomp == 1)
            amrex::Abort("This component of the Bfield is not normal to any boundary !!! \n");
        const int dir = ( (icomp==2) ? (AMREX_SPACEDIM-1) : icomp );
#else
        const int dir = icomp;
#endif
        // Loop over sides, iside = 0 (lo), iside = 1 (hi)
        for (int iside = 0; iside < 2; ++iside) {
            const int is_normal_to_PEC = ( ( iside == 0 )
                                       ? is_boundary_PEC(fbndry_lo, dir)
                                       : is_boundary_PEC(fbndry_hi, dir) );
            if (is_normal_to_PEC == 1) {
                // guard cell outside the domain by "ig" number cells, in direction, dir
                const int ig = ( (iside == 0)
                               ? (dom_lo[dir] - ijk_vec[dir])
                               : (ijk_vec[dir] - (dom_hi[dir] + stag[dir]) ) );
                
                if ( ig == 0 ) OnPECBoundary = 1;

                if ( ig > 0 ) {
                    GuardCell = 1;
                    // Mirror location inside the domain by "ig" number of cells
                    // across PEC boundary in direction, dir, and side, iside
                    ijk_mirror[dir] = ( (iside == 0)
                                   ? (dom_lo[dir] + ig)
                                   : (dom_hi[dir] + stag[dir] - ig));
                    // The normal Bfield in guard cells is equal and opposite to the
                    // value in the respective mirror location across the PEC boundary.
                    sign *= -1._rt;
                }
            }
        }
        if (GuardCell == 1) Bfield(ijk_vec,n) = sign * Bfield(ijk_mirror,n);
        // if ijk_vec is on a PEC boundary in any direction, set Bnormal to 0.
        if (OnPECBoundary == 1) Bfield(ijk_vec,n) = 0._rt;
    }

    /** Returns 1 if any domain boundary is set to PEC, else returns 0.*/
    bool isAnyBoundaryPEC();
    /**
     * \brief Sets the tangential electric field at the PEC boundary to zero.
     *        The guard cell values are set equal and opposite to the valid cell
     *        field value at the respective mirror locations.
     *
     * \param[in,out] Efield     Boundary values of tangential Efield are set to zero.
     * \param[in]     lev        level of the Multifab
     * \param[in]     patch_type coarse or fine
     */
    void ApplyPECtoEfield ( std::array<std::unique_ptr<amrex::MultiFab>, 3>& Efield,
                            const int lev, PatchType patch_type);
    /**
     * \brief Sets the normal component of the magnetic field at the PEC boundary to zero.
     *        The guard cell values are set equal and opposite to the valid cell
     *        field value at the respective mirror locations.
     *
     * \param[in,out] Bfield     Boundary values of normal Bfield are set to zero.
     * \param[in]     lev        level of the Multifab
     * \param[in]     patch_type coarse or fine
     */
    void ApplyPECtoBfield ( std::array<std::unique_ptr<amrex::MultiFab>, 3>& Bfield,
                            const int lev, PatchType patch_type);
}

#endif // WarpX_PEC_KERNELS_H_
