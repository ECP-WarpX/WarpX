/* Copyright 2019 Andrew Myers, David Grote, Maxence Thevenet
 * Remi Lehe, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * This file is a legacy file and will be removed soon.
 * Please do NOT add new bindings here! Please see the other files
 * in this directory for the new pybind11-based bindings!
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_WRAPPERS_H_
#define WARPX_WRAPPERS_H_

#include "Particles/WarpXParticleContainer.H"
#include "Evolve/WarpXDtType.H"
#include <AMReX_Config.H>
#include <AMReX_REAL.H>

#ifdef AMREX_USE_MPI
#   include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

    int warpx_Real_size();
    int warpx_ParticleReal_size();

    int warpx_nSpecies();

    bool warpx_use_fdtd_nci_corr();

    int warpx_galerkin_interpolation();

    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm);

    typedef void(*WARPX_CALLBACK_PY_FUNC_0)();

    void warpx_ConvertLabParamsToBoost();

    void warpx_ReadBCParams();

    void warpx_CheckGriddingForRZSpectral();

    amrex::Real warpx_getCellSize(int dir, int lev);

    amrex::Real warpx_sumParticleCharge(const char* char_species_name, bool local);

  void warpx_ComputeDt ();
  void warpx_MoveWindow (int step, bool move_j);

  void warpx_EvolveE (amrex::Real dt);
  void warpx_EvolveB (amrex::Real dt, DtType a_dt_type);
  void warpx_FillBoundaryE ();
  void warpx_FillBoundaryB ();
  void warpx_SyncCurrent (
      const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
      const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp,
      const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_buffer);
  void warpx_UpdateAuxilaryData ();
  void warpx_PushParticlesandDepose (amrex::Real cur_time);

  void warpx_setistep (int lev, int ii);
  void warpx_sett_new (int lev, amrex::Real time);
  amrex::Real warpx_getdt (int lev);

  int warpx_maxStep ();
  amrex::Real warpx_stopTime ();

  int warpx_finestLevel ();

  void mypc_Redistribute ();

#ifdef __cplusplus
}
#endif

#endif
