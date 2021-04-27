#ifndef WARPX_PEC_KERNELS_H_
#define WARPX_PEC_KERNELS_H_

#include "WarpX.H"
#include "Utils/WarpXAlgorithmSelection.H"


namespace PEC {

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
     *
     * \param[in] icomp        component of the Efield being updated
                               (0=x, 1=y, 2=z in 3D) and (0=x, 1=z in 2D)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] Efield       field data to be updated if (ijk) is at the boundary or a guard cell
     * \param[in] stag         staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetTangentialEfield (const int icomp, const amrex::IntVect & dom_lo,
                                const amrex::IntVect &dom_hi,
                                const amrex::IntVect &ijk_vec,
                                amrex::Array4<amrex::Real> const& Efield,
                                const amrex::IntVect& stag,
                                amrex::GpuArray<int, 3> const& fbndry_lo,
                                amrex::GpuArray<int, 3> const& fbndry_hi )
    {
        // Fields in guard cells set equal and opposite to valid cells
        amrex::IntVect ijk_valid = ijk_vec;
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            // index of dimension
            const int dir = ( ( (icomp+idim)<AMREX_SPACEDIM)? (icomp+idim) : ( (icomp+idim) - AMREX_SPACEDIM));
            // loop over sides , iside = 0 (lo), iside = 1 (hi)
            for (int iside = 0; iside < 2; ++iside) {
                const int is_tangent_to_PEC = ( ( iside == 0 )
                                            ? is_boundary_PEC(fbndry_lo, dir)
                                            : is_boundary_PEC(fbndry_hi, dir) );
                if (is_tangent_to_PEC == 1) {
                    // guard cell outside the domain by "ig" number of cells
                    const int ig = ( (iside == 0)
                                 ? (dom_lo[dir] - ijk_vec[dir])
                                 : (ijk_vec[dir] - (dom_hi[dir]) + stag[dir]) );

                    if ( ig == 0) Efield(ijk_vec) = 0.;

                    if ( ig > 0 ) {

                        // valid cell inside the domain by "ig" number of cells
                        ijk_valid[dir] = ( ( iside == 0)
                                       ? (dom_lo[dir] + ig)
                                       : (dom_hi[dir] + stag[dir] - ig));
                        // setting the guard cell Efield tangential to PEC boundary
                        Efield(ijk_vec) = - Efield(ijk_valid);
                    }
                }
            }
        }
    }

    /**
     * \brief Sets the magnetic field value normal to the PEC boundary to zero. The
     *        field value of the guard cells outside the domain boundary are set equal
     *        and opposite to the field in the valid cells at their mirrored locations.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *
     * \param[in] icomp        component of the Efield being updated
                               (0=x, 1=y, 2=z in 3D) and (0=x, 1=z in 2D)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] Bfield       field data to be updated if (ijk) is at the boundary
                               or a guard cell
     * \param[in] stag         staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetNormalBfield (const int icomp, const amrex::IntVect & dom_lo,
                           const amrex::IntVect & dom_hi,
                           const amrex::IntVect & ijk_vec,
                           amrex::Array4<amrex::Real> const& Bfield,
                           const amrex::IntVect & stag,
                           amrex::GpuArray<int, 3> const& fbndry_lo,
                           amrex::GpuArray<int, 3> const& fbndry_hi )
    {
        amrex::IntVect ijk_valid = ijk_vec;
        const int dir = icomp;
        for (int iside = 0; iside < 2; ++iside) {
            const int is_normal_to_PEC = ( ( iside == 0 )
                                       ? is_boundary_PEC(fbndry_lo, dir)
                                       : is_boundary_PEC(fbndry_hi, dir) );
            if (is_normal_to_PEC == 1) {
                // guard cell outside the domain by "ig" number cells
                const int ig = ( (iside == 0)
                               ? (dom_lo[dir] - ijk_vec[dir])
                               : (ijk_vec[dir] - (dom_hi[dir]) + stag[dir]) );

                if ( ig == 0 ) Bfield(ijk_vec) = 0.;

                if ( ig > 0 ) {

                    // valid cell inside the domain by "ig" number of cells
                    ijk_valid[dir] = ( (iside == 0)
                                   ? (dom_lo[dir] + ig)
                                   : (dom_hi[dir] + stag[dir] - ig));

                    // setting the guard cell Efield tangential to PEC boundary
                    Bfield(ijk_vec) = - Bfield(ijk_valid);
                }
            }
        }
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
