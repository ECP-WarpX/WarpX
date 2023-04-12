#include "BoundaryConditions/WarpX_PerfectConductivity.H"

#include "WarpX.H"

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

using namespace amrex::literals;

bool
PerfectConductivity::isAnyBoundaryPEC() {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC) return true;
        if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC) return true;
    }
    return false;
}

bool
PerfectConductivity::isAnyBoundaryPMC() {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::PMC) return true;
        if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::PMC) return true;
    }
    return false;
}

void
PerfectConductivity::ApplyTangentialZeroToField (std::array<amrex::MultiFab*, 3> field, const int lev,
                                                 FieldBoundaryType boundary_type, PatchType patch_type,
                                                 const bool split_pml_field)
{
    auto& warpx = WarpX::GetInstance();
    amrex::Box domain_box = warpx.Geom(lev).Domain();
    if (patch_type == PatchType::coarse) {
        amrex::IntVect ref_ratio = ( (lev > 0) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1) );
        domain_box.coarsen(ref_ratio);
    }
    amrex::IntVect domain_lo = domain_box.smallEnd();
    amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_lo;
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    amrex::IntVect fx_nodal = field[0]->ixType().toIntVect();
    amrex::IntVect fy_nodal = field[1]->ixType().toIntVect();
    amrex::IntVect fz_nodal = field[2]->ixType().toIntVect();
    amrex::IntVect ng_fieldgather = warpx.get_ng_fieldgather();
    // For each field multifab, apply the boundary condition to ncomponents
    // If not split field, the condition is applied to the regular field used in Maxwell's eq.
    // If split_pml_field is true, then condition is applied to all the split field components of the tangential field.
    int nComp_x = field[0]->nComp();
    int nComp_y = field[1]->nComp();
    int nComp_z = field[2]->nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*field[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Extract field data
        amrex::Array4<amrex::Real> const& fx = field[0]->array(mfi);
        amrex::Array4<amrex::Real> const& fy = field[1]->array(mfi);
        amrex::Array4<amrex::Real> const& fz = field[2]->array(mfi);

        // Extract tileboxes over which to loop
        // if split field, the box includes nodal flag
        // For field used in Maxwell's update, nodal flag plus cells that particles
        // gather fields from in the guard-cell region are included.
        // Note that for simulations without particles or laser, ng_field_gather is 0
        // and the guard-cell values of the field multifab will not be modified.
        amrex::Box const& tex = (split_pml_field) ? mfi.tilebox(field[0]->ixType().toIntVect())
                                                  : mfi.tilebox(field[0]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tey = (split_pml_field) ? mfi.tilebox(field[1]->ixType().toIntVect())
                                                  : mfi.tilebox(field[1]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tez = (split_pml_field) ? mfi.tilebox(field[2]->ixType().toIntVect())
                                                  : mfi.tilebox(field[2]->ixType().toIntVect(), ng_fieldgather);

        // loop over cells and update fields
        amrex::ParallelFor(
            tex, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
                PerfectConductivity::ZeroTangentialField(icomp, domain_lo, domain_hi, iv, n,
                                         fx, fx_nodal, fbndry_lo, fbndry_hi, boundary_type);
            },
            tey, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
                PerfectConductivity::ZeroTangentialField(icomp, domain_lo, domain_hi, iv, n,
                                         fy, fy_nodal, fbndry_lo, fbndry_hi, boundary_type);
            },
            tez, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
                PerfectConductivity::ZeroTangentialField(icomp, domain_lo, domain_hi, iv, n,
                                         fz, fz_nodal, fbndry_lo, fbndry_hi, boundary_type);
            }
        );

    }
}


void
PerfectConductivity::ApplyNormalZeroToField (std::array<amrex::MultiFab*, 3> field, const int lev,
                                             FieldBoundaryType boundary_type, PatchType patch_type,
                                             const bool split_pml_field)
{
    auto& warpx = WarpX::GetInstance();
    amrex::Box domain_box = warpx.Geom(lev).Domain();
    if (patch_type == PatchType::coarse) {
        amrex::IntVect ref_ratio = ( (lev > 0) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1) );
        domain_box.coarsen(ref_ratio);
    }
    amrex::IntVect domain_lo = domain_box.smallEnd();
    amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_lo;
    amrex::GpuArray<FieldBoundaryType, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    amrex::IntVect fx_nodal = field[0]->ixType().toIntVect();
    amrex::IntVect fy_nodal = field[1]->ixType().toIntVect();
    amrex::IntVect fz_nodal = field[2]->ixType().toIntVect();
    amrex::IntVect ng_fieldgather = warpx.get_ng_fieldgather();
    const int nComp_x = field[0]->nComp();
    const int nComp_y = field[1]->nComp();
    const int nComp_z = field[2]->nComp();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*field[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Extract field data
        amrex::Array4<amrex::Real> const& fx = field[0]->array(mfi);
        amrex::Array4<amrex::Real> const& fy = field[1]->array(mfi);
        amrex::Array4<amrex::Real> const& fz = field[2]->array(mfi);

        // Extract tileboxes over which to loop
        // if split field, the box includes nodal flag
        // For field used in Maxwell's update, nodal flag plus cells that particles
        // gather fields from in the guard-cell region are included.
        // Note that for simulations without particles or laser, ng_field_gather is 0
        // and the guard-cell values of the field multifab will not be modified.
        amrex::Box const& tbx = (split_pml_field) ? mfi.tilebox(field[0]->ixType().toIntVect())
                                                  : mfi.tilebox(field[0]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tby = (split_pml_field) ? mfi.tilebox(field[1]->ixType().toIntVect())
                                                  : mfi.tilebox(field[1]->ixType().toIntVect(), ng_fieldgather);
        amrex::Box const& tbz = (split_pml_field) ? mfi.tilebox(field[2]->ixType().toIntVect())
                                                  : mfi.tilebox(field[2]->ixType().toIntVect(), ng_fieldgather);

        // loop over cells and update fields
        amrex::ParallelFor(
            tbx, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
                PerfectConductivity::ZeroNormalField(icomp, domain_lo, domain_hi, iv, n,
                                     fx, fx_nodal, fbndry_lo, fbndry_hi, boundary_type);
            },
            tby, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
                PerfectConductivity::ZeroNormalField(icomp, domain_lo, domain_hi, iv, n,
                                     fy, fy_nodal, fbndry_lo, fbndry_hi, boundary_type);
            },
            tbz, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
                amrex::ignore_unused(k);
#endif
#if (defined WARPX_DIM_1D_Z)
                amrex::ignore_unused(j,k);
#endif
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
                PerfectConductivity::ZeroNormalField(icomp, domain_lo, domain_hi, iv, n,
                                     fz, fz_nodal, fbndry_lo, fbndry_hi, boundary_type);
            }
        );

    }

}
