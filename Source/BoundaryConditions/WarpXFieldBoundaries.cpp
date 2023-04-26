#include "WarpX.H"
#include "BoundaryConditions/PML.H"
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#include "Evolve/WarpXDtType.H"
#include "WarpX_PEC.H"

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <array>
#include <memory>
using namespace amrex::literals;
using namespace amrex;

void
WarpX::ApplyEfieldBoundary (const amrex::IntVect& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        ApplyEfieldBoundary(lev, ng);
    }
}

void
WarpX::ApplyEfieldBoundary (const int lev, const amrex::IntVect& ng)
{
    ApplyEfieldBoundary(lev, PatchType::fine, ng);
    if (lev > 0) ApplyEfieldBoundary(lev, PatchType::coarse, ng);
}

void
WarpX::ApplyEfieldBoundary (const int lev, const PatchType patch_type, const amrex::IntVect& ng)
{
    if (PEC::isAnyBoundaryPEC()) {
        if (patch_type == PatchType::fine) {
            PEC::ApplyPECtoEfield( { get_pointer_Efield_fp(lev, 0),
                                     get_pointer_Efield_fp(lev, 1),
                                     get_pointer_Efield_fp(lev, 2) }, lev, patch_type, ng);
            if (WarpX::isAnyBoundaryPML()) {
                // apply pec on split E-fields in PML region
                const bool split_pml_field = true;
                PEC::ApplyPECtoEfield( pml[lev]->GetE_fp(), lev, patch_type, ng, split_pml_field);
            }
        } else {
            PEC::ApplyPECtoEfield( { get_pointer_Efield_cp(lev, 0),
                                     get_pointer_Efield_cp(lev, 1),
                                     get_pointer_Efield_cp(lev, 2) }, lev, patch_type, ng);
            if (WarpX::isAnyBoundaryPML()) {
                // apply pec on split E-fields in PML region
                const bool split_pml_field = true;
                PEC::ApplyPECtoEfield( pml[lev]->GetE_cp(), lev, patch_type, ng, split_pml_field);
            }
        }
    }
}

void
WarpX::ApplyBfieldBoundary (const amrex::IntVect& ng, const DtType a_dt_type)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        ApplyBfieldBoundary(lev, ng, a_dt_type);
    }
}

void
WarpX::ApplyBfieldBoundary (const int lev, const amrex::IntVect& ng, const DtType a_dt_type)
{
    ApplyBfieldBoundary(lev, PatchType::fine, ng, a_dt_type);
    if (lev > 0) ApplyBfieldBoundary(lev, PatchType::coarse, ng, a_dt_type);
}

void WarpX::ApplyBfieldBoundary (const int lev, const PatchType patch_type, const amrex::IntVect& ng, const DtType a_dt_type)
{
    if (PEC::isAnyBoundaryPEC()) {
        if (patch_type == PatchType::fine) {
            PEC::ApplyPECtoBfield( { get_pointer_Bfield_fp(lev, 0),
                                     get_pointer_Bfield_fp(lev, 1),
                                     get_pointer_Bfield_fp(lev, 2) }, lev, patch_type, ng);
        } else {
            PEC::ApplyPECtoBfield( { get_pointer_Bfield_cp(lev, 0),
                                     get_pointer_Bfield_cp(lev, 1),
                                     get_pointer_Bfield_cp(lev, 2) }, lev, patch_type, ng);
        }
    }

    // Silver-Mueller boundaries are only applied on the first half-push of B
    // This is because the formula used for Silver-Mueller assumes that
    // E and B are staggered in time, which is only true after the first half-push
    if (lev == 0) {
        if (a_dt_type == DtType::FirstHalf) {
            bool applySilverMueller = false;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if ( (WarpX::field_boundary_lo[idim] == FieldBoundaryType::Absorbing_SilverMueller) ||
                   (WarpX::field_boundary_hi[idim] == FieldBoundaryType::Absorbing_SilverMueller) ) {
                    applySilverMueller = true;
                }
            }
            if(applySilverMueller) m_fdtd_solver_fp[0]->ApplySilverMuellerBoundary(
                                         Efield_fp[lev], Bfield_fp[lev],
                                         Geom(lev).Domain(), dt[lev],
                                         WarpX::field_boundary_lo,
                                         WarpX::field_boundary_hi);
        }
    }
}

void WarpX::ApplyRhofieldBoundary (const int lev, MultiFab* rho,
                                   PatchType patch_type)
{
    if (PEC::isAnyBoundaryPEC()) PEC::ApplyPECtoRhofield(rho, lev, patch_type);
}

void WarpX::ApplyJfieldBoundary (const int lev, amrex::MultiFab* Jx,
                                 amrex::MultiFab* Jy, amrex::MultiFab* Jz,
                                 PatchType patch_type)
{
    if (PEC::isAnyBoundaryPEC()) PEC::ApplyPECtoJfield(Jx, Jy, Jz, lev, patch_type);
}
