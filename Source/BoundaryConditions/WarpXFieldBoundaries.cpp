#include "WarpX.H"
#include "BoundaryConditions/PML.H"
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#include "Evolve/WarpXDtType.H"
#include "WarpX_PerfectConductivity.H"

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

void WarpX::ApplyEfieldBoundary(const int lev, PatchType patch_type)
{
    std::array<amrex::MultiFab*, 3> Evector;
    std::array<amrex::MultiFab*, 3> EPMLvector;
    if (patch_type == PatchType::fine) {
        Evector = { get_pointer_Efield_fp(lev, 0), get_pointer_Efield_fp(lev, 1), get_pointer_Efield_fp(lev, 2) };
        if (WarpX::isAnyBoundaryPML() && pml[lev]) {
            EPMLvector = pml[lev]->GetE_fp();
        }
    } else {
        Evector = { get_pointer_Efield_cp(lev, 0), get_pointer_Efield_cp(lev, 1), get_pointer_Efield_cp(lev, 2) };
        if (WarpX::isAnyBoundaryPML() && pml[lev]) {
            EPMLvector = pml[lev]->GetE_cp();
        }
    }

    if (PerfectConductivity::isAnyBoundaryPEC()) {
        PerfectConductivity::ApplyTangentialZeroToField(Evector, lev, FieldBoundaryType::PEC, patch_type);
        if (WarpX::isAnyBoundaryPML() && pml[lev]) {
            // apply pec on split E-fields in PML region
            const bool split_pml_field = true;
            PerfectConductivity::ApplyTangentialZeroToField( EPMLvector, lev, FieldBoundaryType::PEC, patch_type, split_pml_field);
        }
    }

    if (PerfectConductivity::isAnyBoundaryPMC()) {
        PerfectConductivity::ApplyNormalZeroToField(Evector, lev, FieldBoundaryType::PMC, patch_type);
        if (WarpX::isAnyBoundaryPML() && pml[lev]) {
            // apply pec on split E-fields in PML region
            const bool split_pml_field = true;
            PerfectConductivity::ApplyNormalZeroToField( EPMLvector, lev, FieldBoundaryType::PMC, patch_type, split_pml_field);
        }
    }

}

void WarpX::ApplyBfieldBoundary (const int lev, PatchType patch_type, DtType a_dt_type)
{
    std::array<amrex::MultiFab*, 3> Bvector;
    if (patch_type == PatchType::fine) {
        Bvector = { get_pointer_Bfield_fp(lev, 0), get_pointer_Bfield_fp(lev, 1), get_pointer_Bfield_fp(lev, 2) };
    } else {
        Bvector = { get_pointer_Bfield_cp(lev, 0), get_pointer_Bfield_cp(lev, 1), get_pointer_Bfield_cp(lev, 2) };
    }
    if (PerfectConductivity::isAnyBoundaryPEC()) {
        PerfectConductivity::ApplyNormalZeroToField(Bvector, lev, FieldBoundaryType::PEC, patch_type);
    }
    if (PerfectConductivity::isAnyBoundaryPMC()) {
        PerfectConductivity::ApplyTangentialZeroToField(Bvector, lev, FieldBoundaryType::PMC, patch_type);
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
