#ifndef WARPX_PEC_KERNELS_H_
#define WARPX_PEC_KERNELS_H_

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

namespace PEC {
using namespace amrex;
    /**
     * \brief Determines if the field boundary condition stored in fboundary
     *        in direction, dir, is PEC.
     *
     * \param[in] fboundary  Value containing boundary type
     * \param[in] dir        direction
     *
     * \returns 1 if the boundary type is PEC else 0
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool is_boundary_PEC (amrex::GpuArray<int, 3> const& fboundary, int dir) {
        return ( fboundary[dir] == FieldBoundaryType::PEC );
    }

    /**
     * \brief Determines if the particle boundary condition stored in pboundary
     *        in direction, dir, is reflecting.
     *
     * \param[in] pboundary  Value containing boundary type
     * \param[in] dir        direction
     *
     * \returns 1 if the boundary type is reflecting else 0
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool is_boundary_reflecting (amrex::GpuArray<ParticleBoundaryType, 3> const& pboundary, int dir) {
        return ( pboundary[dir] == ParticleBoundaryType::Reflecting );
    }

    /**
     * \brief Calculates the number of grid points the given index is pass the
     *        domain boundary i.e. a value of +1 means the current cell is
     *        outside of the simulation domain by 1 cell. Note that the high
     *        side domain boundary is between cell dom_hi and dom_hi+1 for cell
     *        centered grids and on cell dom_hi+1 for nodal grid. This is why
     *        (dom_hi[idim] + is_nodal[idim]) is used below.
     *
     * \param[in] dom_lo, dom_hi  Domain boundaries
     * \param[in] ijk_vec         Cell coordinates
     * \param[in] is_nodal        Whether the field of interest is nodal
     * \param[in] idim            Dimension of interest
     * \param[in] iside           0 for low and 1 for high
     *
     * \returns number of grid points to the boundary
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int get_cell_count_to_boundary (const amrex::IntVect& dom_lo,
        const amrex::IntVect& dom_hi, const amrex::IntVect& ijk_vec,
        const amrex::IntVect& is_nodal, const int idim, const int iside)
    {
        return ((iside == 0) ? (dom_lo[idim] - ijk_vec[idim])
                             : (ijk_vec[idim] - (dom_hi[idim] + is_nodal[idim])));
    }

    /**
     * \brief Sets the electric field value tangential to the PEC boundary to zero. The
     *        tangential Efield components in the guard cells outside the
     *        domain boundary are set equal and opposite to the field in the valid cells
     *        at their mirrored locations. The normal Efield components in the guard cells
     *        are set equal to the field in the valid cells at their mirrored locations.
     *        The number or depth of guard cells updated is equal to the shape factor of
     *        particles in each dimension.
     *        For corner cells with mixed boundaries, the mirror location could be outside
     *        valid region, while still ensuring PEC condition is maintained across the
     *        PEC boundary, and the necessary sign change is accounted for depending on
     *        if the component, icomp, is tangential or normal to the PEC boundary.
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
     * \param[in] icomp        component of the Efield being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] Efield       field data to be updated if (ijk) is at the boundary or a guard cell
     * \param[in] is_nodal     staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetEfieldOnPEC (const int icomp, const amrex::IntVect & dom_lo,
                                const amrex::IntVect &dom_hi,
                                const amrex::IntVect &ijk_vec, const int n,
                                amrex::Array4<amrex::Real> const& Efield,
                                const amrex::IntVect& is_nodal,
                                amrex::GpuArray<int, 3> const& fbndry_lo,
                                amrex::GpuArray<int, 3> const& fbndry_hi )
    {
        // Tangential Efield componentes in guard cells set equal and opposite to cells
        // in the mirror locations across the PEC boundary, whereas normal E-field
        // components are set equal to values in the mirror locations across the PEC
        // boundary. Here we just initialize it.
        amrex::IntVect ijk_mirror = ijk_vec;
        bool OnPECBoundary = false;
        bool GuardCell = false;
        amrex::Real sign = 1._rt;
        // Loop over all the dimensions
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Loop over sides, iside = 0 (lo), iside = 1 (hi)
            for (int iside = 0; iside < 2; ++iside) {
                const bool isPECBoundary = ( (iside == 0)
                                        ? is_boundary_PEC(fbndry_lo, idim)
                                        : is_boundary_PEC(fbndry_hi, idim) );
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                // For 2D : for icomp==1, (Ey in XZ, Etheta in RZ),
                //          icomp=1 is tangential to both x and z boundaries
                //          The logic below ensures that the flags are set right for 2D
                const bool is_tangent_to_PEC = ( (icomp == AMREX_SPACEDIM*idim)
                                             ? false : true );
#elif (defined WARPX_DIM_1D_Z)
                // For 1D : icomp=0 and icomp=1 (Ex and Ey are tangential to the z boundary)
                //          The logic below ensures that the flags are set right for 1D
                const bool is_tangent_to_PEC = ( ( icomp == idim+2) ? false : true );
#else
                const bool is_tangent_to_PEC = ( ( icomp == idim) ? false : true );
#endif
                if (isPECBoundary == true) {
                    // grid point ijk_vec is ig number of points pass the
                    // domain boundary in direction, idim
                    const int ig = get_cell_count_to_boundary(
                        dom_lo, dom_hi, ijk_vec, is_nodal, idim, iside);

                    if (ig == 0) {
                        if (is_tangent_to_PEC == true and is_nodal[idim] == 1) {
                            OnPECBoundary = true;
                        }
                    } else if (ig > 0) {
                        // Find mirror location across PEC boundary
                        ijk_mirror[idim] = ( ( iside == 0)
                                        ? (dom_lo[idim] + ig - (1 - is_nodal[idim]))
                                        : (dom_hi[idim] + 1 - ig));
                        GuardCell = true;
                        // tangential components are inverted across PEC boundary
                        if (is_tangent_to_PEC == true) sign *= -1._rt;
                    }
                } // is PEC boundary
            } // loop over iside
        } // loop over dimensions
        if (OnPECBoundary == true) {
            // if ijk_vec is on a PEC boundary in any direction, set Etangential to 0.
            Efield(ijk_vec,n) = 0._rt;
        } else if (GuardCell == true) {
            Efield(ijk_vec,n) = sign * Efield(ijk_mirror,n);
        }
    }


    /**
     * \brief Sets the magnetic field value normal to the PEC boundary to zero. The
     *        tangential (and normal) field value of the guard cells outside the
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
     * \param[in] icomp        component of the Bfield being updated
     *                         (0=x, 1=y, 2=z in Cartesian)
     *                         (0=r, 1=theta, 2=z in RZ)
     * \param[in] dom_lo       index value of the lower domain boundary (cell-centered)
     * \param[in] dom_hi       index value of the higher domain boundary (cell-centered)
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of Efield Array4
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] Bfield       field data to be updated if (ijk) is at the boundary
                               or a guard cell
     * \param[in] is_nodal     staggering of the field data being updated.
     * \param[in] fbndry_lo    Field boundary type at the lower boundaries
     * \param[in] fbndry_hi    Field boundary type at the upper boundaries
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetBfieldOnPEC (const int icomp, const amrex::IntVect & dom_lo,
                           const amrex::IntVect & dom_hi,
                           const amrex::IntVect & ijk_vec, const int n,
                           amrex::Array4<amrex::Real> const& Bfield,
                           const amrex::IntVect & is_nodal,
                           amrex::GpuArray<int, 3> const& fbndry_lo,
                           amrex::GpuArray<int, 3> const& fbndry_hi )
    {
        amrex::IntVect ijk_mirror = ijk_vec;
        bool OnPECBoundary = false;
        bool GuardCell = false;
        amrex::Real sign = 1._rt;
        // Loop over all dimensions
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Loop over sides, iside = 0 (lo), iside = 1 (hi)
            for (int iside = 0; iside < 2; ++iside) {
                const bool isPECBoundary = ( (iside == 0 )
                                        ? is_boundary_PEC(fbndry_lo, idim)
                                        : is_boundary_PEC(fbndry_hi, idim) );
                if (isPECBoundary == true) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                    // For 2D : for icomp==1, (By in XZ, Btheta in RZ),
                    //          icomp=1 is not normal to x or z boundary
                    //          The logic below ensures that the flags are set right for 2D
                    const bool is_normal_to_PEC = ( (icomp == AMREX_SPACEDIM*idim)
                                                ? true : false );
#elif (defined WARPX_DIM_1D_Z)
                    // For 1D : icomp=0 and icomp=1 (Bx and By are not normal to the z boundary)
                    //          The logic below ensures that the flags are set right for 1D
                    const bool is_normal_to_PEC = ( ( icomp == idim+2) ? true : false );
#else
                    const bool is_normal_to_PEC = ( ( icomp == idim) ? true : false );
#endif

                    // grid point ijk_vec is ig number of points pass the
                    // domain boundary in direction, idim
                    const int ig = get_cell_count_to_boundary(
                        dom_lo, dom_hi, ijk_vec, is_nodal, idim, iside);

                    if (ig == 0) {
                        // Only normal component is set to 0
                        if (is_normal_to_PEC == true and is_nodal[idim]==1) {
                            OnPECBoundary = true;
                        }
                    } else if ( ig > 0) {
                        // Mirror location inside the domain by "ig" number of cells
                        // across PEC boundary in direction, idim, and side, iside
                        ijk_mirror[idim] = ( (iside == 0)
                                        ? (dom_lo[idim] + ig - (1 - is_nodal[idim]))
                                        : (dom_hi[idim] + 1 - ig));
                        GuardCell = true;
                        // Sign of the normal component in guard cell is inverted
                        if (is_normal_to_PEC == true) sign *= -1._rt;
                    }
                } // if PEC Boundary
            } // loop over sides
        } // loop of dimensions

        if (OnPECBoundary == true) {
            // if ijk_vec is on a PEC boundary in any direction, set Bnormal to 0.
            Bfield(ijk_vec,n) = 0._rt;
        } else if (GuardCell == true) {
            // Bnormal and Btangential is set opposite and equal to the value
            // in the mirror location, respectively.
            Bfield(ijk_vec,n) = sign * Bfield(ijk_mirror,n);
        }
    }


    /**
     * \brief Sets the rho or J field value in cells close to and on a PEC boundary. The
     *        charge/current density deposited in the guard cells are either reflected
     *        back into the simulation domain (if a reflecting particle
     *        boundary is used), or the opposite charge/current density is deposited
     *        back in the domain to capture the effect of an image charge.
     *        The charge/current density on the PEC boundary is set to 0 while values
     *        in the guard cells are set equal (and opposite) to their mirror
     *        location inside the domain - representing image charges - in the
     *        normal (tangential) direction.
     *
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of the rho Array4
     * \param[in out] field    field data to be updated
     * \param[in] mirrorfac    mirror cell is given by mirrorfac - ijk_vec
     * \param[in] psign        Whether the field value should be flipped across the boundary
     * \param[in] is_pec       Whether the given boundary is PEC
     * \param[in] tangent_to_bndy    Whether a given direction is perpendicular to the boundary
     * \param[in] fabbox       multifab box including ghost cells
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetRhoOrJfieldFromPEC (const int n,
                                const amrex::IntVect & ijk_vec,
                                amrex::Array4<amrex::Real> const& field,
                                amrex::GpuArray<GpuArray<int, 2>, AMREX_SPACEDIM> const& mirrorfac,
                                amrex::GpuArray<GpuArray<amrex::Real, 2>, AMREX_SPACEDIM> const& psign,
                                amrex::GpuArray<GpuArray<bool, 2>, AMREX_SPACEDIM> const& is_pec,
                                amrex::GpuArray<bool, AMREX_SPACEDIM> const& tangent_to_bndy,
                                amrex::Box const& fabbox)
    {
        // The boundary is handled in 2 steps:
        // 1) The cells internal to the domain are updated using the
        //    current deposited in the guard cells
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            for (int iside = 0; iside < 2; ++iside)
            {
                if (!is_pec[idim][iside]) continue;

                // Get the mirror guard cell index
                amrex::IntVect iv_mirror = ijk_vec;
                iv_mirror[idim] = mirrorfac[idim][iside] - ijk_vec[idim];

                // On the PEC boundary the current density is set to 0
                if (ijk_vec == iv_mirror) field(ijk_vec, n) = 0._rt;
                // otherwise update the internal cell if the mirror guard cell exists
                else if (fabbox.contains(iv_mirror))
                {
                    field(ijk_vec,n) += psign[idim][iside] * field(iv_mirror,n);
                }
            }
        }
        // 2) The guard cells are updated with the appropriate image
        //    charge current based on the current in the valid cells
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            for (int iside = 0; iside < 2; ++iside)
            {
                if (!is_pec[idim][iside]) continue;

                amrex::IntVect iv_mirror = ijk_vec;
                iv_mirror[idim] = mirrorfac[idim][iside] - ijk_vec[idim];
                if (ijk_vec != iv_mirror && fabbox.contains(iv_mirror))
                {
                    if (tangent_to_bndy[idim])
                        field(iv_mirror, n) = -field(ijk_vec, n);
                    else
                        field(iv_mirror, n) = field(ijk_vec, n);
                }
            }
        }
    }


    /**
     * \brief This function sets the given field value on a PEC boundary
     *        to enforce a Neumann boundary condition (zero derivative) in the
     *        normal direction.
     *
     * \param[in] n            index of the MultiFab component being updated
     * \param[in] ijk_vec      indices along the x(i), y(j), z(k) of the rho Array4
     * \param[in out] field    field data to be updated
     * \param[in] mirrorfac    mirror cell is given by mirrorfac - ijk_vec
     * \param[in] psign        Whether the field value should be flipped across the boundary
     * \param[in] is_pec       Whether the given boundary is PEC
     * \param[in] tangent_to_bndy    Whether a given direction is perpendicular to the boundary
     * \param[in] fabbox       multifab box including ghost cells
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void SetNeumannOnPEC (const int n,
                          const amrex::IntVect & ijk_vec,
                          amrex::Array4<amrex::Real> const& field,
                          amrex::GpuArray<GpuArray<int, 2>, AMREX_SPACEDIM> const& mirrorfac,
                          amrex::GpuArray<GpuArray<bool, 2>, AMREX_SPACEDIM> const& is_pec,
                          amrex::Box const& fabbox )
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            for (int iside = 0; iside < 2; ++iside)
            {
                if (!is_pec[idim][iside]) continue;

                // Get the mirror guard cell index
                amrex::IntVect iv_mirror = ijk_vec;
                iv_mirror[idim] = mirrorfac[idim][iside] - ijk_vec[idim];

                // On the PEC boundary the field value is set equal to the
                // first value in the domain (nodal fields)
                if (ijk_vec == iv_mirror) {
                    iv_mirror[idim] += (iside == 0) ? 1 : -1;
                    if (fabbox.contains(iv_mirror)) field(ijk_vec, n) = field(iv_mirror, n);
                }
                // otherwise set the mirror guard cell equal to the internal cell value
                else if (fabbox.contains(iv_mirror))
                {
                    field(iv_mirror, n) = field(ijk_vec, n);
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
     * \param[in,out] Efield          Boundary values of tangential Efield are set to zero
     * \param[in]     lev             level of the Multifab
     * \param[in]     patch_type      coarse or fine
     * \param[in]     split_pml_field whether pml the multifab is the regular Efield or
     *                                split pml field
     */
    void ApplyPECtoEfield ( std::array<amrex::MultiFab*, 3> Efield,
                            const int lev, PatchType patch_type,
                            const bool split_pml_field = false);
    /**
     * \brief Sets the normal component of the magnetic field at the PEC boundary to zero.
     *        The guard cell values are set equal and opposite to the valid cell
     *        field value at the respective mirror locations.
     *
     * \param[in,out] Bfield     Boundary values of normal Bfield are set to zero.
     * \param[in]     lev        level of the Multifab
     * \param[in]     patch_type coarse or fine
     */
    void ApplyPECtoBfield ( std::array<amrex::MultiFab*, 3> Bfield,
                            const int lev, PatchType patch_type);

    /**
     * \brief Reflects charge density deposited over the PEC boundary back into
     * the simulation domain.
     *
     * \param[in,out] rho        Multifab containing the charge density
     * \param[in]     lev        level of the Multifab
     * \param[in]     patch_type coarse or fine
     */
    void ApplyPECtoRhofield(amrex::MultiFab* rho, const int lev,
                            PatchType patch_type);

    /**
     * \brief Reflects current density deposited over the PEC boundary back into
     * the simulation domain.
     *
     * \param[in,out] Jx, Jy, Jz Multifabs containing the current density
     * \param[in]     lev        level of the Multifab
     * \param[in]     patch_type coarse or fine
     */
    void ApplyPECtoJfield(amrex::MultiFab* Jx, amrex::MultiFab* Jy,
                          amrex::MultiFab* Jz, const int lev,
                          PatchType patch_type);

    /**
     * \brief Apply the PEC boundary to the electron pressure field.
     *
     * \param[in,out] Pefield    Multifab containing the electron pressure
     * \param[in]     lev        level of the Multifab
     * \param[in]     patch_type coarse or fine
     */
    void ApplyPECtoElectronPressure (amrex::MultiFab* Pefield,
                                     const int lev, PatchType patch_type);
}
#endif // WarpX_PEC_KERNELS_H_
