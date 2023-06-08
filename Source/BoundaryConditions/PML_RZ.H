/* Copyright 2021 David Grote
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PML_RZ_H_
#define WARPX_PML_RZ_H_

#include "PML_RZ_fwd.H"

#ifdef WARPX_USE_PSATD
#   include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#endif

#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <string>

enum struct PatchType : int;

class PML_RZ
{
public:
    PML_RZ (const int lev, const amrex::BoxArray& grid_ba, const amrex::DistributionMapping& grid_dm,
            const amrex::Geometry* geom, const int ncell, const int do_pml_in_domain);

    void ApplyDamping(amrex::MultiFab* Et_fp, amrex::MultiFab* Ez_fp,
                      amrex::MultiFab* Bt_fp, amrex::MultiFab* Bz_fp,
                      amrex::Real dt);

    std::array<amrex::MultiFab*,2> GetE_fp ();
    std::array<amrex::MultiFab*,2> GetB_fp ();

#ifdef WARPX_USE_PSATD
    void PushPSATD (const int lev);
#endif

    void FillBoundaryE ();
    void FillBoundaryB ();
    void FillBoundaryE (PatchType patch_type);
    void FillBoundaryB (PatchType patch_type);

    void CheckPoint (const std::string& dir) const;
    void Restart (const std::string& dir);

    ~PML_RZ () = default;

private:

    const int m_ncell;
    const int m_do_pml_in_domain;
    const amrex::Geometry* m_geom;

    // Only contains Er and Et, and Br and Bt
    std::array<std::unique_ptr<amrex::MultiFab>,2> pml_E_fp;
    std::array<std::unique_ptr<amrex::MultiFab>,2> pml_B_fp;

#ifdef WARPX_USE_PSATD
    void PushPMLPSATDSinglePatchRZ ( const int lev,
                SpectralSolverRZ& solver,
                std::array<std::unique_ptr<amrex::MultiFab>,2>& pml_E,
                std::array<std::unique_ptr<amrex::MultiFab>,2>& pml_B);
#endif

};

#endif
