#include "WarpX.H"

#include "WarpX_PEC.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <memory>

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

