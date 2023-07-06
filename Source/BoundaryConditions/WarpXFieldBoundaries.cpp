#include "WarpX.H"
#include "BoundaryConditions/PML.H"
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
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

void WarpX::ApplyEfieldBoundary(const int lev, PatchType patch_type)
{
    if (PEC::isAnyBoundaryPEC()) {
        if (patch_type == PatchType::fine) {
            PEC::ApplyPECtoEfield( { get_pointer_Efield_fp(lev, 0),
                                     get_pointer_Efield_fp(lev, 1),
                                     get_pointer_Efield_fp(lev, 2) }, lev, patch_type);
            if (WarpX::isAnyBoundaryPML()) {
                // apply pec on split E-fields in PML region
                const bool split_pml_field = true;
                PEC::ApplyPECtoEfield( pml[lev]->GetE_fp(), lev, patch_type, split_pml_field);
            }
        } else {
            PEC::ApplyPECtoEfield( { get_pointer_Efield_cp(lev, 0),
                                     get_pointer_Efield_cp(lev, 1),
                                     get_pointer_Efield_cp(lev, 2) }, lev, patch_type);
            if (WarpX::isAnyBoundaryPML()) {
                // apply pec on split E-fields in PML region
                const bool split_pml_field = true;
                PEC::ApplyPECtoEfield( pml[lev]->GetE_cp(), lev, patch_type, split_pml_field);
            }
        }
    }

#ifdef WARPX_DIM_RZ
    if (patch_type == PatchType::fine) {
        ApplyFieldBoundaryOnAxis(get_pointer_Efield_fp(lev, 0),
                                 get_pointer_Efield_fp(lev, 1),
                                 get_pointer_Efield_fp(lev, 2), lev);
    } else {
        ApplyFieldBoundaryOnAxis(get_pointer_Efield_cp(lev, 0),
                                 get_pointer_Efield_cp(lev, 1),
                                 get_pointer_Efield_cp(lev, 2), lev);
    }
#endif
}

void WarpX::ApplyBfieldBoundary (const int lev, PatchType patch_type, DtType a_dt_type)
{
    if (PEC::isAnyBoundaryPEC()) {
        if (patch_type == PatchType::fine) {
            PEC::ApplyPECtoBfield( { get_pointer_Bfield_fp(lev, 0),
                                     get_pointer_Bfield_fp(lev, 1),
                                     get_pointer_Bfield_fp(lev, 2) }, lev, patch_type);
        } else {
            PEC::ApplyPECtoBfield( { get_pointer_Bfield_cp(lev, 0),
                                     get_pointer_Bfield_cp(lev, 1),
                                     get_pointer_Bfield_cp(lev, 2) }, lev, patch_type);
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

#ifdef WARPX_DIM_RZ
    if (patch_type == PatchType::fine) {
        ApplyFieldBoundaryOnAxis(get_pointer_Bfield_fp(lev, 0),
                                 get_pointer_Bfield_fp(lev, 1),
                                 get_pointer_Bfield_fp(lev, 2), lev);
    } else {
        ApplyFieldBoundaryOnAxis(get_pointer_Bfield_cp(lev, 0),
                                 get_pointer_Bfield_cp(lev, 1),
                                 get_pointer_Bfield_cp(lev, 2), lev);
    }
#endif
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

#ifdef WARPX_DIM_RZ
// Applies the boundary conditions that are specific to the axis when in RZ.
void
WarpX::ApplyFieldBoundaryOnAxis (amrex::MultiFab* Er, amrex::MultiFab* Et, amrex::MultiFab* Ez, int lev)
{
    const amrex::IntVect ngE = get_ng_fieldgather();
    /* const amrex::IntVect ngE = Er->nGrowVect(); */

    constexpr int NODE = amrex::IndexType::NODE;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(*Er, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

        amrex::Box const & tilebox = mfi.tilebox();

        // Lower corner of tile box physical domain
        // Note that this is done before the tilebox.grow so that
        // these do not include the guard cells.
        const std::array<amrex::Real, 3>& xyzmin = LowerCorner(tilebox, lev, 0._rt);
        const amrex::Real rmin = xyzmin[0];

        // Skip blocks that don't touch the axis
        if (rmin > 0._rt) continue;

        amrex::Array4<amrex::Real> const& Er_arr = Er->array(mfi);
        amrex::Array4<amrex::Real> const& Et_arr = Et->array(mfi);
        amrex::Array4<amrex::Real> const& Ez_arr = Ez->array(mfi);

        amrex::Box tbr = amrex::convert( tilebox, Er->ixType().toIntVect() );
        amrex::Box tbt = amrex::convert( tilebox, Et->ixType().toIntVect() );
        amrex::Box tbz = amrex::convert( tilebox, Ez->ixType().toIntVect() );

        // For ishift, 1 means cell centered, 0 means node centered
        int const ishift_r = (tbr.type(0) != NODE);
        int const ishift_t = (tbt.type(0) != NODE);
        int const ishift_z = (tbz.type(0) != NODE);

        // Set tileboxes to only include the axis guard cells
        // (including the corners in z).
        tbr.setRange(0, -ngE[0], ngE[0]);
        tbt.setRange(0, -ngE[0], ngE[0]);
        tbz.setRange(0, -ngE[0], ngE[0]);
        tbr.grow(1, ngE[1]);
        tbt.grow(1, ngE[1]);
        tbz.grow(1, ngE[1]);

        const int nmodes = n_rz_azimuthal_modes;

        amrex::ParallelFor(tbr, tbt, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
        {
            Er_arr(i,j,0,0) = -Er_arr(-i-ishift_r,j,0,0);

            for (int imode=1 ; imode < nmodes ; imode++) {
                Er_arr(i,j,0,2*imode-1) = std::pow(-1, imode+1)*Er_arr(-i-ishift_r,j,0,2*imode-1);
                Er_arr(i,j,0,2*imode) = std::pow(-1, imode+1)*Er_arr(-i-ishift_r,j,0,2*imode);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
        {
            Et_arr(i,j,0,0) = -Et_arr(-i-ishift_t,j,0,0);

            for (int imode=1 ; imode < nmodes ; imode++) {
                Et_arr(i,j,0,2*imode-1) = std::pow(-1, imode+1)*Et_arr(-i-ishift_t,j,0,2*imode-1);
                Et_arr(i,j,0,2*imode) = std::pow(-1, imode+1)*Et_arr(-i-ishift_t,j,0,2*imode);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/)
        {
            Ez_arr(i,j,0,0) = Ez_arr(-i-ishift_z,j,0,0);

            for (int imode=1 ; imode < nmodes ; imode++) {
                Ez_arr(i,j,0,2*imode-1) = -std::pow(-1, imode+1)*Ez_arr(-i-ishift_z,j,0,2*imode-1);
                Ez_arr(i,j,0,2*imode) = -std::pow(-1, imode+1)*Ez_arr(-i-ishift_z,j,0,2*imode);
            }

        });
    }
}
#endif

void WarpX::ApplyElectronPressureBoundary (const int lev, PatchType patch_type)
{
    if (PEC::isAnyBoundaryPEC()) {
        if (patch_type == PatchType::fine) {
            PEC::ApplyPECtoElectronPressure(
                m_hybrid_pic_model->get_pointer_electron_pressure_fp(lev), lev, patch_type
            );
        } else {
            amrex::Abort(Utils::TextMsg::Err(
            "ApplyElectronPressureBoundary: Only one level implemented for hybrid solver."));
        }
    }
}
