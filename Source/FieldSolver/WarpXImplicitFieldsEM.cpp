/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "Evolve/WarpXDtType.H"
#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#if defined(WARPX_USE_PSATD)
#   include "FieldSolver/SpectralSolver/SpectralFieldData.H"
#   ifdef WARPX_DIM_RZ
#       include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#       include "BoundaryConditions/PML_RZ.H"
#   else
#       include "FieldSolver/SpectralSolver/SpectralSolver.H"
#   endif
#endif
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpXPushFieldsEM_K.H"
#include "WarpX_FDTD.H"

#include <AMReX.H>
#ifdef AMREX_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_MFIter.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <cmath>
#include <memory>

using namespace amrex;

void
WarpX::ComputeRHSB (amrex::Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        ComputeRHSB(lev, a_dt);
    }
}

void
WarpX::ComputeRHSB (int lev, amrex::Real a_dt)
{
    WARPX_PROFILE("WarpX::ComputeRHSB()");
    ComputeRHSB(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        ComputeRHSB(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::ComputeRHSB (int lev, PatchType patch_type, amrex::Real a_dt)
{

    // set RHS to zero value
    Bfield_rhs[lev][0]->setVal(0.0);
    Bfield_rhs[lev][1]->setVal(0.0);
    Bfield_rhs[lev][2]->setVal(0.0);

    // Compute Bfield_rhs in regular cells by calling EvolveB
    if (patch_type == PatchType::fine) {
        m_fdtd_solver_fp[lev]->EvolveB(Bfield_rhs[lev], Efield_fp[lev], G_fp[lev],
                                       m_face_areas[lev], m_area_mod[lev], ECTRhofield[lev], Venl[lev],
                                       m_flag_info_face[lev], m_borrowing[lev], lev, a_dt);
    } else {
        m_fdtd_solver_cp[lev]->EvolveB(Bfield_rhs[lev], Efield_cp[lev], G_cp[lev],
                                       m_face_areas[lev], m_area_mod[lev], ECTRhofield[lev], Venl[lev],
                                       m_flag_info_face[lev], m_borrowing[lev], lev, a_dt);
    }

    // Compute Bfield_rhs in PML cells by calling EvolveBPML
    if (do_pml && pml[lev]->ok()) {
        if (patch_type == PatchType::fine) {
            //m_fdtd_solver_fp[lev]->EvolveBPML(
            //    pml[lev]->GetRHSB_fp(), pml[lev]->GetE_fp(), a_dt, WarpX::do_dive_cleaning);
        } else {
            //m_fdtd_solver_cp[lev]->EvolveBPML(
            //    pml[lev]->GetRHSB_cp(), pml[lev]->GetE_cp(), a_dt, WarpX::do_dive_cleaning);
        }
    }

}


void
WarpX::ComputeRHSE (amrex::Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        ComputeRHSE(lev, a_dt);
    }
}

void
WarpX::ComputeRHSE (int lev, amrex::Real a_dt)
{
    WARPX_PROFILE("WarpX::ComputeRHSE()");
    ComputeRHSE(lev, PatchType::fine, a_dt);
    if (lev > 0)
    {
        ComputeRHSE(lev, PatchType::coarse, a_dt);
    }
}

void
WarpX::ComputeRHSE (int lev, PatchType patch_type, amrex::Real a_dt)
{
    // set RHS to zero value
    Efield_rhs[lev][0]->setVal(0.0);
    Efield_rhs[lev][1]->setVal(0.0);
    Efield_rhs[lev][2]->setVal(0.0);

    // Compute Efield_rhs in regular cells by calling EvolveE
    if (patch_type == PatchType::fine) {
        m_fdtd_solver_fp[lev]->EvolveE( Efield_rhs[lev], Bfield_fp[lev],
                                        current_fp[lev], m_edge_lengths[lev],
                                        m_face_areas[lev], ECTRhofield[lev],
                                        F_fp[lev], lev, a_dt );
    } else {
        m_fdtd_solver_cp[lev]->EvolveE( Efield_rhs[lev], Bfield_cp[lev],
                                        current_cp[lev], m_edge_lengths[lev],
                                        m_face_areas[lev], ECTRhofield[lev],
                                        F_cp[lev], lev, a_dt );
    }

    // Compute Efield_rhs in PML cells by calling EvolveEPML
    if (do_pml && pml[lev]->ok()) {
        if (patch_type == PatchType::fine) {
            //m_fdtd_solver_fp[lev]->EvolveEPML(
            //    pml[lev]->GetRHSE_fp(), pml[lev]->GetB_fp(),
            //    pml[lev]->Getj_fp(), pml[lev]->Get_edge_lengths(),
            //    pml[lev]->GetF_fp(),
            //    pml[lev]->GetMultiSigmaBox_fp(),
            //    a_dt, pml_has_particles );
        } else {
            //m_fdtd_solver_cp[lev]->EvolveEPML(
            //    pml[lev]->GetRHSE_cp(), pml[lev]->GetB_cp(),
            //    pml[lev]->Getj_cp(), pml[lev]->Get_edge_lengths(),
            //    pml[lev]->GetF_cp(),
            //    pml[lev]->GetMultiSigmaBox_cp(),
            //    a_dt, pml_has_particles );
        }
    }

}
