#ifndef WARPX_PERFECTCONDUCTIVITY_KERNELS_H_
#define WARPX_PERFECTCONDUCTIVITY_KERNELS_H_

#include "WarpX.H"
#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <memory>

namespace PerfectConductivity {
using namespace amrex;

    /**
     * \brief Sets field value tangential to the boundary with boundary_type to zero.
     *        This is the PEC condition for the electric field, and the
     *        PMC condition for the magnetic field.
     *        The tangential field components in the guard cells outside the
     *        domain boundary are set equal and opposite to the field in the valid cells
     *        at their mirrored locations. The normal field components in the guard cells
     *        are set equal to the field in the valid cells at their mirrored locations.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *        For corner cells with mixed boundaries, the mirror location could be outside
     *        valid region, while still ensuring the condition is maintained across the
     *        the boundary, and the necessary sign change is accounted for depending on
     *        if the component, icomp, is tangential or normal to the the boundary.
     *
     *        For 3D :
     *            x component is tangential to the y-boundary and z-boundary
     *            y component is tangential to the x-boundary and z-boundary
     *            z component is tangential to the x-boundary and y-boundary
     *            x component is normal to the x-boundary
     *            y component is normal to the y-boundary
     *            z component is normal to the z-boundary
     *            where, x-boundary is the yz-plane at x=xmin and x=xmax
     *                   y-boundary is the xz-plane at y=ymin and y=ymax
     *                   z-boundary is the xy-plane at z=zmin and z=zmax
     *
     *        For 2D : WarpX uses X-Z as the two dimensions
     *            x component is tangential to the z-boundary
     *            y component is tangential to the x-boundary and z-boundary
     *            z component is tangential to the x-boundary
     *            x component is normal to the x-boundary
     *            y component is not normal to any boundary (Only xz dimensions in 2D)
     *            z component is normal to the z-boundary
     *            where, x-boundary is along the line z at x=xmin and x=xmax
     *                   z-boundary is along the line x at z=zmin and z=zmax
     *
     *        For 1D : WarpX uses Z as the only dimension
     *            x component is tangential to the z-boundary
     *            y component is tangential to the z-boundary
     *            z component is not tangential to the z-boundary
     *            x component is not normal to any boundary (Only z dimension in 1D)
     *            y component is not normal to any boundary (Only z dimension in 1D)
     *            z component is normal to the z-boundary
     *            where, z-boundary is a point at z=zmin and z=zmax
     *
     *        For RZ : WarpX uses R-Z as the two dimensions
     *            r component is tangential to the z-boundary
     *            theta_component is tangential to the r-boundary and z-boundary
     *            z component is tangential to the r-boundary
     *            r component is normal to the r-boundary
     *            theta_component is not normal to any boundary (on RZ dimensions are modeled)
     *            z component is normal to the z-boundary
     *            where, r-boundary is along the line z at r=rmin and r=rmax
     *                   z-boundary is along the line r at z=zmin and z=zmax
     *
     *
     * \param[in] icomp        component of the field being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of field Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] field        field data to be updated if (ijk) is at the boundary or a guard cell
     * \param[in] is_nodal     staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     * \param[in] boundary_type The targeted boundary type (either PEC or PMC)
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void ZeroTangentialField (const int icomp, const amrex::IntVect & dom_lo,
                              const amrex::IntVect &dom_hi,
                              const amrex::IntVect &ijk_vec, const int n,
                              amrex::Array4<amrex::Real> const& field,
                              const amrex::IntVect& is_nodal,
                              amrex::GpuArray<FieldBoundaryType, 3> const& fbndry_lo,
                              amrex::GpuArray<FieldBoundaryType, 3> const& fbndry_hi,
                              FieldBoundaryType const boundary_type)
    {
        // Tangential field componentes in guard cells are set equal and opposite to cells
        // in the mirror locations across the boundary, whereas normal field
        // components are set equal to values in the mirror locations across the
        // boundary. Here we just initialize it.
        amrex::IntVect ijk_mirror = ijk_vec;
        bool OnBoundary = false;
        bool GuardCell = false;
        amrex::Real sign = 1._rt;
        // Loop over all the dimensions
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Loop over sides, iside = 0 (lo), iside = 1 (hi)
            // Loop over sides, iside = 0 (lo), iside = 1 (hi)
            for (int iside = 0; iside < 2; ++iside) {
                const bool isBoundary = ( (iside == 0)
                                        ? fbndry_lo[idim] == boundary_type
                                        : fbndry_hi[idim] == boundary_type);
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                // For 2D : for icomp==1, (Ey in XZ, Etheta in RZ),
                //          icomp=1 is tangential to both x and z boundaries
                //          The logic below ensures that the flags are set right for 2D
                const bool is_tangent = ( (icomp == AMREX_SPACEDIM*idim) ? false : true );
#elif (defined WARPX_DIM_1D_Z)
                // For 1D : icomp=0 and icomp=1 (Ex and Ey are tangential to the z boundary)
                //          The logic below ensures that the flags are set right for 1D
                const bool is_tangent = ( ( icomp == idim+2) ? false : true );
#else
                const bool is_tangent = ( ( icomp == idim) ? false : true );
#endif
                if (isBoundary == true) {
                    // guard cell outside the domain by "ig" number cells, in direction, idim
                    const int ig = ( (iside == 0)
                                 ? (dom_lo[idim] - ijk_vec[idim])
                                 : (ijk_vec[idim] - (dom_hi[idim] + is_nodal[idim]) ) );
                    if (ig == 0) {
                        if (is_tangent == true and is_nodal[idim] == 1) {
                            OnBoundary = true;
                        }
                    } else if (ig > 0) {
                        // Find mirror location across the boundary
                        ijk_mirror[idim] = ( ( iside == 0)
                                        ? (dom_lo[idim] + ig)
                                        : (dom_hi[idim] + is_nodal[idim] - ig));
                        GuardCell = true;
                        // tangential components are inverted across the boundary
                        if (is_tangent == true) sign *= -1._rt;
                    }
                } // is boundary
            } // loop over iside
        } // loop over dimensions
        if (OnBoundary == true) {
            // if ijk_vec is on a the boundary in any direction, set tangential field to 0.
            field(ijk_vec,n) = 0._rt;
        } else if (GuardCell == true) {
            field(ijk_vec,n) = sign * field(ijk_mirror,n);
        }
    }

    /**
     * \brief Sets field value normal to the boundary with boundary_type to zero.
     *        This is the PMC condition for the electric field, and the
     *        PEC condition for the magnetic field.
     *        The tangential (and normal) field value of the guard cells outside the
     *        domain boundary are set equal (and opposite) to the respective field components
     *        in the valid cells at their mirrored locations.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *
     *        For 3D :
     *            x component is tangential to the y-boundary and z-boundary
     *            y component is tangential to the x-boundary and z-boundary
     *            z component is tangential to the x-boundary and y-boundary
     *            x component is normal to the x-boundary
     *            y component is normal to the y-boundary
     *            z component is normal to the z-boundary
     *            where, x-boundary is the yz-plane at x=xmin and x=xmax
     *                   y-boundary is the xz-plane at y=ymin and y=ymax
     *                   z-boundary is the xy-plane at z=zmin and z=zmax
     *
     *        For 2D : WarpX uses X-Z as the two dimensions
     *            x component is tangential to the z-boundary
     *            y component is tangential to the x-boundary and z-boundary
     *            z component is tangential to the x-boundary
     *            x component is normal to the x-boundary
     *            y component is not normal to any boundary (Only xz dimensions in 2D)
     *            z component is normal to the z-boundary
     *            where, x-boundary is along the line z at x=xmin and x=xmax
     *                   z-boundary is along the line x at z=zmin and z=zmax
     *
     *        For 1D : WarpX uses Z as the only dimension
     *            x component is tangential to the z-boundary
     *            y component is tangential to the z-boundary
     *            z component is not tangential to the z-boundary
     *            x component is not normal to any boundary (Only z dimension in 1D)
     *            y component is not normal to any boundary (Only z dimension in 1D)
     *            z component is normal to the z-boundary
     *            where, z-boundary is a point at z=zmin and z=zmax
     *
     *        For RZ : WarpX uses R-Z as the two dimensions
     *            r component is tangential to the z-boundary
     *            theta_component is tangential to the r-boundary and z-boundary
     *            z component is tangential to the r-boundary
     *            r component is normal to the r-boundary
     *            theta_component is not normal to any boundary (on RZ dimensions are modeled)
     *            z component is normal to the z-boundary
     *            where, r-boundary is along the line z at r=rmin and r=rmax
     *                   z-boundary is along the line r at z=zmin and z=zmax
     *
     *
     * \param[in] icomp        component of the field being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of field Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] field        field data to be updated if (ijk) is at the boundary
                               or a guard cell
     * \param[in] is_nodal     staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     * \param[in] boundary_type The targeted boundary type (either PEC or PMC)
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void ZeroNormalField (const int icomp, const amrex::IntVect & dom_lo,
                          const amrex::IntVect & dom_hi,
                          const amrex::IntVect & ijk_vec, const int n,
                          amrex::Array4<amrex::Real> const& field,
                          const amrex::IntVect & is_nodal,
                          amrex::GpuArray<FieldBoundaryType, 3> const& fbndry_lo,
                          amrex::GpuArray<FieldBoundaryType, 3> const& fbndry_hi,
                          FieldBoundaryType const boundary_type)
    {
        amrex::IntVect ijk_mirror = ijk_vec;
        bool OnBoundary = false;
        bool GuardCell = false;
        amrex::Real sign = 1._rt;
        // Loop over all dimensions
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Loop over sides, iside = 0 (lo), iside = 1 (hi)
            for (int iside = 0; iside < 2; ++iside) {
                const bool isBoundary = ( (iside == 0 )
                                        ? fbndry_lo[idim] == boundary_type
                                        : fbndry_hi[idim] == boundary_type );
                if (isBoundary == true) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                    // For 2D : for icomp==1, (By in XZ, Btheta in RZ),
                    //          icomp=1 is not normal to x or z boundary
                    //          The logic below ensures that the flags are set right for 2D
                    const bool is_normal = ( (icomp == AMREX_SPACEDIM*idim) ? true : false );
#elif (defined WARPX_DIM_1D_Z)
                    // For 1D : icomp=0 and icomp=1 (Bx and By are not normal to the z boundary)
                    //          The logic below ensures that the flags are set right for 1D
                    const bool is_normal = ( ( icomp == idim+2) ? true : false );
#else
                    const bool is_normal = ( ( icomp == idim) ? true : false );
#endif
                    // guard cell outside the domain by "ig" number cells, in direction, idim
                    const int ig = ( (iside == 0)
                                   ? (dom_lo[idim] - ijk_vec[idim])
                                   : (ijk_vec[idim] - (dom_hi[idim] + is_nodal[idim]) ) );
                    if (ig == 0) {
                        // Only normal component is set to 0
                        if (is_normal == true and is_nodal[idim]==1) {
                            OnBoundary = true;
                        }
                    } else if ( ig > 0) {
                        // Mirror location inside the domain by "ig" number of cells
                        // across the boundary in direction, idim, and side, iside
                        ijk_mirror[idim] = ( (iside == 0)
                                       ? (dom_lo[idim] + ig)
                                       : (dom_hi[idim] + is_nodal[idim] - ig));
                        GuardCell = true;
                        // Sign of the normal component in guard cell is inverted
                        if (is_normal == true) sign *= -1._rt;
                    }
                } // if Boundary
            } // loop over sides
        } // loop of dimensions

        if (OnBoundary == true) {
            // if ijk_vec is on a the boundary in any direction, set Bnormal to 0.
            field(ijk_vec,n) = 0._rt;
        } else if (GuardCell == true) {
            // Bnormal and Btangential is set opposite and equal to the value
            // in the mirror location, respectively.
            field(ijk_vec,n) = sign * field(ijk_mirror,n);
        }
    }

    /** Returns 1 if any domain boundary is set to PEC, else returns 0.*/
    bool isAnyBoundaryPEC();

    /** Returns 1 if any domain boundary is set to PMC, else returns 0.*/
    bool isAnyBoundaryPMC();

    /**
     * \brief Sets the tangential field at boundaries with boundary_type to zero.
     *        The guard cell values are set equal and opposite to the valid cell
     *        field value at the respective mirror locations.
     *
     * \param[in,out] field           Boundary values of tangential field are set to zero
     * \param[in]     lev             level of the Multifab
     * \param[in]     boundary_type   The targeted boundary type (either PEC or PMC)
     * \param[in]     patch_type      coarse or fine
     * \param[in]     split_pml_field whether pml the multifab is the regular field or
     *                                split pml field
     */
    void ApplyTangentialZeroToField (std::array<amrex::MultiFab*, 3> field, const int lev,
                                     FieldBoundaryType boundary_type, PatchType patch_type,
                                     const bool split_pml_field = false);
    /**
     * \brief Sets the normal component of the magnetic field at the PEC boundary to zero.
     *        The guard cell values are set equal and opposite to the valid cell
     *        field value at the respective mirror locations.
     *
     * \param[in,out] field      Boundary values of normal field are set to zero.
     * \param[in]     lev        level of the Multifab
     * \param[in]     boundary_type   The targeted boundary type (either PEC or PMC)
     * \param[in]     patch_type coarse or fine
     * \param[in]     split_pml_field whether pml the multifab is the regular field or
     *                                split pml field
     */
    void ApplyNormalZeroToField (std::array<amrex::MultiFab*, 3> field, const int lev,
                                 FieldBoundaryType boundary_type, PatchType patch_type,
                                 const bool split_pml_field = false);
}

#endif // WARPX_PERFECTCONDUCTIVITY_KERNELS_H_
