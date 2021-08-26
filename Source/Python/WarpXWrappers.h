/* Copyright 2019 Andrew Myers, David Grote, Maxence Thevenet
 * Remi Lehe, Weiqun Zhang
 *
 * This file is part of WarpX.
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

    int warpx_nComps();

    int warpx_nCompsSpecies(const char* char_species_name);

    int warpx_SpaceDim();

    void amrex_init (int argc, char* argv[]);

#ifdef AMREX_USE_MPI
    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm);
#endif

    void amrex_finalize (int finalize_mpi);

    void warpx_init ();

    void warpx_finalize ();

    typedef void(*WARPX_CALLBACK_PY_FUNC_0)();

    void warpx_set_callback_py_afterinit (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_beforeEsolve (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterEsolve (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_beforedeposition (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterdeposition (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_particlescraper (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_particleloader (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_beforestep (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterstep (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterrestart (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_particleinjection (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_appliedfields (WARPX_CALLBACK_PY_FUNC_0);

    void warpx_evolve (int numsteps);  // -1 means the inputs parameter will be used.

    void warpx_addNParticles(const char* char_species_name,
                             int lenx,
                             amrex::ParticleReal const * x,
                             amrex::ParticleReal const * y,
                             amrex::ParticleReal const * z,
                             amrex::ParticleReal const * vx,
                             amrex::ParticleReal const * vy,
                             amrex::ParticleReal const * vz,
                             int nattr,
                             amrex::ParticleReal const * attr,
                             int uniqueparticles);

    void warpx_ConvertLabParamsToBoost();

    void warpx_ReadBCParams();

    void warpx_CheckGriddingForRZSpectral();

    amrex::Real warpx_getProbLo(int dir);

    amrex::Real warpx_getProbHi(int dir);

    long warpx_getNumParticles(const char* char_species_name);

    amrex::ParticleReal** warpx_getParticleStructs(
        const char* char_species_name, int lev, int* num_tiles,
        int** particles_per_tile);

    amrex::ParticleReal** warpx_getParticleArrays(
        const char* char_species_name, const char* char_comp_name, int lev,
        int* num_tiles, int** particles_per_tile);

    int warpx_getParticleCompIndex(
        const char* char_species_name, const char* char_comp_name);

    void warpx_addRealComp(
        const char* char_species_name, const char* char_comp_name, bool comm);

  void warpx_ComputeDt ();
  void warpx_MoveWindow (int step, bool move_j);

  void warpx_EvolveE (amrex::Real dt);
  void warpx_EvolveB (amrex::Real dt, DtType a_dt_type);
  void warpx_FillBoundaryE ();
  void warpx_FillBoundaryB ();
  void warpx_SyncCurrent ();
  void warpx_UpdateAuxilaryData ();
  void warpx_PushParticlesandDepose (amrex::Real cur_time);

  int warpx_getistep (int lev);
  void warpx_setistep (int lev, int ii);
  amrex::Real warpx_gett_new (int lev);
  void warpx_sett_new (int lev, amrex::Real time);
  amrex::Real warpx_getdt (int lev);

  int warpx_maxStep ();
  amrex::Real warpx_stopTime ();

  int warpx_finestLevel ();

  int warpx_getMyProc ();
  int warpx_getNProcs ();


  void mypc_Redistribute ();

#ifdef __cplusplus
}
#endif

#endif
