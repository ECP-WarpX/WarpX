#include "WarpX.H"
#include "WarpX_PEC.H"
#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
using namespace amrex::literals;

void WarpX::ApplyEfieldBoundary(const int lev, PatchType patch_type)
{    
    if (PEC::isAnyBoundaryPEC()) {
        if (patch_type == PatchType::fine) {
            PEC::ApplyPECtoEfield( Efield_fp[lev], lev, patch_type);
        } else {
            PEC::ApplyPECtoEfield( Efield_cp[lev], lev, patch_type);
        }
    }
}

void WarpX::ApplyBfieldBoundary (const int lev, PatchType patch_type)
{
    if (PEC::isAnyBoundaryPEC()) {
        if (patch_type == PatchType::fine) {
            PEC::ApplyPECtoBfield( Bfield_fp[lev], lev, patch_type);
        } else {
            PEC::ApplyPECtoBfield( Bfield_cp[lev], lev, patch_type);
        }
    }
}







//    amrex::Print() << " CALL PEC is " << callPEC << "\n";
//    amrex::Box const& domain_box = Geom(0).Domain();
//    amrex::IntVect domain_lo = domain_box.smallEnd();
//    amrex::IntVect domain_hi = domain_box.bigEnd();
//    amrex::Print() << " domain box : " << domain_box << "\n";
//    amrex::Print() << " domain_lo : " << domain_lo[0] << " " << domain_lo[1] << " " << domain_lo[2] << "\n";
//    amrex::Print() << " domain_hi : " << domain_hi[0] << " " << domain_hi[1] << " " << domain_hi[2] << "\n";
//    amrex::IntVect Ex_stag = Efield[0]->ixType().toIntVect();
//    amrex::IntVect Ey_stag = Efield[1]->ixType().toIntVect();
//    amrex::IntVect Ez_stag = Efield[2]->ixType().toIntVect();
//#ifdef AMREX_USE_OMP
//#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
//#endif
//    for (amrex::MFIter mfi(*Efield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//        // Extract field data
//        amrex::Array4<amrex::Real> const& Ex = Efield[0]->array(mfi);
//        amrex::Array4<amrex::Real> const& Ey = Efield[1]->array(mfi);
//        amrex::Array4<amrex::Real> const& Ez = Efield[2]->array(mfi);
//
//        // Extract tileboxes for which to loop
//        amrex::Box const& tex = mfi.tilebox(Efield[0]->ixType().toIntVect());
//        amrex::Box const& tey = mfi.tilebox(Efield[1]->ixType().toIntVect());
//        amrex::Box const& tez = mfi.tilebox(Efield[2]->ixType().toIntVect());
//
//        // loop over cells and update fields
//        amrex::ParallelFor(tex, tey, tez,
//            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
//                amrex::IntVect iv(AMREX_D_DECL(i,j,k));                    
//                const int icomp = 0;
//                PEC::ZeroTangentialEfield(icomp, domain_lo, domain_hi, iv, Ex(i,j,k), Ex_stag);
//            },
//            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
//                amrex::IntVect iv(AMREX_D_DECL(i,j,k));                    
//                const int icomp = 1;
//                PEC::ZeroTangentialEfield(icomp, domain_lo, domain_hi, iv, Ey(i,j,k), Ey_stag);
//            },
//            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
//                amrex::IntVect iv(AMREX_D_DECL(i,j,k));                    
//                const int icomp = 2;
//                PEC::ZeroTangentialEfield(icomp, domain_lo, domain_hi, iv, Ez(i,j,k), Ez_stag);
//            }
//        );
//
//    }
//}

