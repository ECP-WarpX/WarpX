#include "BoundaryConditions/WarpX_PEC.H"
#include "WarpX.H"
#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
using namespace amrex::literals;


bool
PEC::isAnyBoundaryPEC() {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if ( WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC) return true;
        if ( WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC) return true;
    }
    return false;
}

void
PEC::ApplyPECtoEfield (std::array<std::unique_ptr<amrex::MultiFab>, 3>& Efield, const int lev,
                       PatchType patch_type)
{
    auto& warpx = WarpX::GetInstance();
    amrex::Box domain_box = warpx.Geom(lev).Domain();
    amrex::IntVect ref_ratio = amrex::IntVect(1);
    if (patch_type == PatchType::coarse) {
        if (lev > 0 ) ref_ratio = WarpX::RefRatio(lev-1);
        domain_box.coarsen(ref_ratio);
    }
    amrex::IntVect domain_lo = domain_box.smallEnd();
    amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::IntVect shape_factor(AMREX_D_DECL(WarpX::nox, WarpX::noy, WarpX::noz));
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    shape_factor[1] = WarpX::noz;
#endif
    amrex::GpuArray<int, 3> fbndry_lo;
    amrex::GpuArray<int, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    amrex::IntVect Ex_stag = Efield[0]->ixType().toIntVect();
    amrex::IntVect Ey_stag = Efield[1]->ixType().toIntVect();
    amrex::IntVect Ez_stag = Efield[2]->ixType().toIntVect();
    int nComp_x = Efield[0]->nComp();
    int nComp_y = Efield[1]->nComp();
    int nComp_z = Efield[2]->nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Efield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Extract field data
        amrex::Array4<amrex::Real> const& Ex = Efield[0]->array(mfi);
        amrex::Array4<amrex::Real> const& Ey = Efield[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Ez = Efield[2]->array(mfi);

        // Extract tileboxes for which to loop
        amrex::Box const& tex = mfi.tilebox(Efield[0]->ixType().toIntVect(), shape_factor);
        amrex::Box const& tey = mfi.tilebox(Efield[1]->ixType().toIntVect(), shape_factor);
        amrex::Box const& tez = mfi.tilebox(Efield[2]->ixType().toIntVect(), shape_factor);

        // loop over cells and update fields
        amrex::ParallelFor(
            tex, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
                PEC::SetTangentialEfield(icomp, domain_lo, domain_hi, iv, n,
                                           Ex, Ex_stag, fbndry_lo, fbndry_hi);
            },
            tey, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
                PEC::SetTangentialEfield(icomp, domain_lo, domain_hi, iv, n,
                                           Ey, Ey_stag, fbndry_lo, fbndry_hi);
            },
            tez, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
                PEC::SetTangentialEfield(icomp, domain_lo, domain_hi, iv, n,
                                           Ez, Ez_stag, fbndry_lo, fbndry_hi);
            }
        );

    }
}


void
PEC::ApplyPECtoBfield (std::array<std::unique_ptr<amrex::MultiFab>, 3>& Bfield, const int lev,
                       PatchType patch_type)
{
    auto& warpx = WarpX::GetInstance();
    amrex::Box domain_box = warpx.Geom(lev).Domain();
    amrex::IntVect ref_ratio = amrex::IntVect(1);
    if (patch_type == PatchType::coarse) {
        if (lev > 0 ) ref_ratio = WarpX::RefRatio(lev-1);
        domain_box.coarsen(ref_ratio);
    }
    amrex::IntVect domain_lo = domain_box.smallEnd();
    amrex::IntVect domain_hi = domain_box.bigEnd();
    amrex::IntVect shape_factor(AMREX_D_DECL(WarpX::nox, WarpX::noy, WarpX::noz));
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    shape_factor[1] = WarpX::noz;
#endif
    amrex::GpuArray<int, 3> fbndry_lo;
    amrex::GpuArray<int, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    amrex::IntVect Bx_stag = Bfield[0]->ixType().toIntVect();
    amrex::IntVect By_stag = Bfield[1]->ixType().toIntVect();
    amrex::IntVect Bz_stag = Bfield[2]->ixType().toIntVect();
    const int nComp_x = Bfield[0]->nComp();
    const int nComp_y = Bfield[1]->nComp();
    const int nComp_z = Bfield[2]->nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Bfield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Extract field data
        amrex::Array4<amrex::Real> const& Bx = Bfield[0]->array(mfi);
        amrex::Array4<amrex::Real> const& By = Bfield[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Bz = Bfield[2]->array(mfi);

        // Extract tileboxes for which to loop
        amrex::Box const& tbx = mfi.tilebox(Bfield[0]->ixType().toIntVect(), shape_factor);
        amrex::Box const& tby = mfi.tilebox(Bfield[1]->ixType().toIntVect(), shape_factor);
        amrex::Box const& tbz = mfi.tilebox(Bfield[2]->ixType().toIntVect(), shape_factor);

        // loop over cells and update fields
        amrex::ParallelFor(
            tbx, nComp_x,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
                PEC::SetNormalBfield(icomp, domain_lo, domain_hi, iv, n,
                                     Bx, Bx_stag, fbndry_lo, fbndry_hi);
            },
#if (AMREX_SPACEDIM == 3)
            tby, nComp_y,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
                PEC::SetNormalBfield(icomp, domain_lo, domain_hi, iv, n,
                                     By, By_stag, fbndry_lo, fbndry_hi);
            },
#endif
            tbz, nComp_z,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
                PEC::SetNormalBfield(icomp, domain_lo, domain_hi, iv, n,
                                     Bz, Bz_stag, fbndry_lo, fbndry_hi);
            }
        );

    }

}
