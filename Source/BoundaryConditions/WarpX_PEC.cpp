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
#if (AMREX_SPACEDIM == 2)
    shape_factor[1] = WarpX::noz;
#endif
    amrex::GpuArray<int, 3> fbndry_lo;
    amrex::GpuArray<int, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    amrex::Print() << " domain box : " << domain_box << "\n";
    amrex::Print() << " domain_lo : " << domain_lo[0] << " " << domain_lo[1] << " " << domain_lo[2] << "\n";
    amrex::Print() << " domain_hi : " << domain_hi[0] << " " << domain_hi[1] << " " << domain_hi[2] << "\n";
    amrex::IntVect Ex_stag = Efield[0]->ixType().toIntVect();
    amrex::IntVect Ey_stag = Efield[1]->ixType().toIntVect();
    amrex::IntVect Ez_stag = Efield[2]->ixType().toIntVect();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Efield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Extract field data
        amrex::Array4<amrex::Real> const& Ex = Efield[0]->array(mfi);
        amrex::Array4<amrex::Real> const& Ey = Efield[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Ez = Efield[2]->array(mfi);

        // Extract tileboxes for which to loop
        amrex::Box const& tex = mfi.tilebox(Efield[0]->ixType().toIntVect());
        amrex::Box const& tey = mfi.tilebox(Efield[1]->ixType().toIntVect());
        amrex::Box const& tez = mfi.tilebox(Efield[2]->ixType().toIntVect());

        // loop over cells and update fields
        amrex::ParallelFor(tex, tey, tez,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
//                PEC::ZeroTangentialEfield(icomp, domain_lo, domain_hi, iv, Ex(i,j,k), Ex_stag);
                PEC::SetTangentialEfield(icomp, domain_lo, domain_hi, iv, shape_factor,
                                           Ex, Ex_stag, fbndry_lo, fbndry_hi);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
//                PEC::ZeroTangentialEfield(icomp, domain_lo, domain_hi, iv, Ey(i,j,k), Ey_stag);
                PEC::SetTangentialEfield(icomp, domain_lo, domain_hi, iv, shape_factor,
                                           Ey, Ey_stag, fbndry_lo, fbndry_hi);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
//                PEC::ZeroTangentialEfield(icomp, domain_lo, domain_hi, iv, Ez(i,j,k), Ez_stag);
                PEC::SetTangentialEfield(icomp, domain_lo, domain_hi, iv, shape_factor,
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
#if (AMREX_SPACEDIM == 2)
    shape_factor[1] = WarpX::noz;
#endif
    amrex::GpuArray<int, 3> fbndry_lo;
    amrex::GpuArray<int, 3> fbndry_hi;
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        fbndry_lo[idim] = WarpX::field_boundary_lo[idim];
        fbndry_hi[idim] = WarpX::field_boundary_hi[idim];
    }
    amrex::Print() << " domain box : " << domain_box << "\n";
    amrex::Print() << " domain_lo : " << domain_lo[0] << " " << domain_lo[1] << " " << domain_lo[2] << "\n";
    amrex::Print() << " domain_hi : " << domain_hi[0] << " " << domain_hi[1] << " " << domain_hi[2] << "\n";
    amrex::IntVect Bx_stag = Bfield[0]->ixType().toIntVect();
    amrex::IntVect By_stag = Bfield[1]->ixType().toIntVect();
    amrex::IntVect Bz_stag = Bfield[2]->ixType().toIntVect();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Bfield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Extract field data
        amrex::Array4<amrex::Real> const& Bx = Bfield[0]->array(mfi);
        amrex::Array4<amrex::Real> const& By = Bfield[1]->array(mfi);
        amrex::Array4<amrex::Real> const& Bz = Bfield[2]->array(mfi);

        // Extract tileboxes for which to loop
        amrex::Box const& tbx = mfi.tilebox(Bfield[0]->ixType().toIntVect());
        amrex::Box const& tby = mfi.tilebox(Bfield[1]->ixType().toIntVect());
        amrex::Box const& tbz = mfi.tilebox(Bfield[2]->ixType().toIntVect());

        // loop over cells and update fields
        amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 0;
                //PEC::ZeroNormalBfield(icomp, domain_lo, domain_hi, iv, Bx(i,j,k), Bx_stag);
                PEC::SetNormalBfield(icomp, domain_lo, domain_hi, iv, shape_factor,
                                     Bx, Bx_stag, fbndry_lo, fbndry_hi);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 1;
                //PEC::ZeroNormalBfield(icomp, domain_lo, domain_hi, iv, By(i,j,k), By_stag);
                PEC::SetNormalBfield(icomp, domain_lo, domain_hi, iv, shape_factor,
                                     By, By_stag, fbndry_lo, fbndry_hi);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                const int icomp = 2;
                //PEC::ZeroNormalBfield(icomp, domain_lo, domain_hi, iv, Bz(i,j,k), Bz_stag);
                PEC::SetNormalBfield(icomp, domain_lo, domain_hi, iv, shape_factor,
                                     Bz, Bz_stag, fbndry_lo, fbndry_hi);
            }
        );

    }

}
