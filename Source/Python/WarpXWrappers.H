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

    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm);

    void amrex_finalize (int finalize_mpi);

    void warpx_init ();

    void warpx_finalize ();

    typedef void(*WARPX_CALLBACK_PY_FUNC_0)();

    void warpx_set_callback_py (const char* char_callback_name,
                                WARPX_CALLBACK_PY_FUNC_0 callback);
    void warpx_clear_callback_py (const char* char_callback_name);

    void warpx_evolve (int numsteps);  // -1 means the inputs parameter will be used.

    void warpx_addNParticles(const char* char_species_name,
                             int lenx,
                             amrex::ParticleReal const * x,
                             amrex::ParticleReal const * y,
                             amrex::ParticleReal const * z,
                             amrex::ParticleReal const * vx,
                             amrex::ParticleReal const * vy,
                             amrex::ParticleReal const * vz,
                             const int nattr_real,
                             amrex::ParticleReal const * attr_real,
                             const int nattr_int,
                             int const * attr_int,
                             int uniqueparticles);

    void warpx_ConvertLabParamsToBoost();

    void warpx_ReadBCParams();

    void warpx_CheckGriddingForRZSpectral();

    amrex::Real warpx_getProbLo(int dir);

    amrex::Real warpx_getProbHi(int dir);

    amrex::Real warpx_getCellSize(int dir, int lev);

    long warpx_getNumParticles(const char* char_species_name, const bool local);

    void warpx_convert_id_to_long (amrex::Long* ids, const WarpXParticleContainer::ParticleType* pstructs, int size);

    void warpx_convert_cpu_to_int (int* cpus, const WarpXParticleContainer::ParticleType* pstructs, int size);

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

    amrex::Real warpx_sumParticleCharge(const char* char_species_name, const bool local);

    int warpx_getParticleBoundaryBufferSize(const char* species_name, int boundary);

    int** warpx_getParticleBoundaryBufferScrapedSteps(
        const char* species_name, int boundary, int lev,
        int* num_tiles, int** particles_per_tile);

    amrex::ParticleReal** warpx_getParticleBoundaryBuffer(
        const char* species_name, int boundary, int lev,
        int* num_tiles, int** particles_per_tile, const char* comp_name);

    amrex::ParticleReal** warpx_getParticleBoundaryBufferStructs(
            const char* species_name, int boundary, int lev,
            int* num_tiles, int** particles_per_tile);

    void warpx_clearParticleBoundaryBuffer ();

    /**
     * \brief This function is used to deposit a given species' charge density
     * in the rho_fp multifab which can then be accessed from python via
     * pywarpx.fields.RhoFPWrapper()
     *
     * @param[in] species_name specifying the name of the species to deposit
     * @param[in] lev mesh refinement level
     */
    void warpx_depositChargeDensity (const char* species_name, int lev);

  void warpx_ComputeDt ();
  void warpx_MoveWindow (int step, bool move_j);

  void warpx_EvolveE (amrex::Real dt);
  void warpx_EvolveB (amrex::Real dt, DtType a_dt_type);
  void warpx_FillBoundaryE ();
  void warpx_FillBoundaryB ();
  void warpx_SyncRho ();
  void warpx_SyncCurrent (
      const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
      const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp);
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

  amrex::Real** warpx_getEfield (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getEfieldCP (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getEfieldFP (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);

  amrex::Real** warpx_getBfield (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getBfieldCP (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getBfieldFP (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);

  amrex::Real** warpx_getCurrentDensity (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getCurrentDensityCP (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getCurrentDensityFP (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);

  int* warpx_getEfieldLoVects (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getEfieldCPLoVects (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getEfieldFPLoVects (int lev, int direction, int *return_size, int **ngrowvect);

  int* warpx_getBfieldLoVects (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getBfieldCPLoVects (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getBfieldFPLoVects (int lev, int direction, int *return_size, int **ngrowvect);

  int* warpx_getCurrentDensityLoVects (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getCurrentDensityCPLoVects (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getCurrentDensityFPLoVects (int lev, int direction, int *return_size, int **ngrowvect);

  amrex::Real** warpx_getEdgeLengths (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  int* warpx_getEdgeLengthsLoVects (int lev, int direction, int *return_size, int **ngrowvect);

  amrex::Real** warpx_getFaceAreas (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  int* warpx_getFaceAreasLoVects (int lev, int direction, int *return_size, int **ngrowvect);

  int* warpx_getEx_nodal_flag ();
  int* warpx_getEy_nodal_flag ();
  int* warpx_getEz_nodal_flag ();
  int* warpx_getBx_nodal_flag ();
  int* warpx_getBy_nodal_flag ();
  int* warpx_getBz_nodal_flag ();
  int* warpx_getJx_nodal_flag ();
  int* warpx_getJy_nodal_flag ();
  int* warpx_getJz_nodal_flag ();
  int* warpx_getRho_nodal_flag ();
  int* warpx_getPhi_nodal_flag ();
  int* warpx_getF_nodal_flag ();
  int* warpx_getG_nodal_flag ();
  int* warpx_get_edge_lengths_x_nodal_flag ();
  int* warpx_get_edge_lengths_y_nodal_flag ();
  int* warpx_get_edge_lengths_z_nodal_flag ();
  int* warpx_get_face_areas_x_nodal_flag ();
  int* warpx_get_face_areas_y_nodal_flag ();
  int* warpx_get_face_areas_z_nodal_flag ();

  amrex::Real** warpx_getChargeDensityCP (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getChargeDensityFP (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  int* warpx_getChargeDensityCPLoVects (int lev, int *return_size, int **ngrowvect);
  int* warpx_getChargeDensityFPLoVects (int lev, int *return_size, int **ngrowvect);

  amrex::Real** warpx_getPhiFP (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);

  int* warpx_getPhiFPLoVects (int lev, int *return_size, int **ngrowvect);

  amrex::Real** warpx_getFfieldCP (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getFfieldFP (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  int* warpx_getFfieldCPLoVects (int lev, int *return_size, int **ngrowvect);
  int* warpx_getFfieldFPLoVects (int lev, int *return_size, int **ngrowvect);
  amrex::Real** warpx_getGfieldCP (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getGfieldFP (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  int* warpx_getGfieldCPLoVects (int lev, int *return_size, int **ngrowvect);
  int* warpx_getGfieldFPLoVects (int lev, int *return_size, int **ngrowvect);

  amrex::Real** warpx_getEfieldCP_PML (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getEfieldFP_PML (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getBfieldCP_PML (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getBfieldFP_PML (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getCurrentDensityCP_PML (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getCurrentDensityFP_PML (int lev, int direction, int *return_size, int *ncomps, int **ngrowvect, int **shapes);

  int* warpx_getEfieldCPLoVects_PML (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getEfieldFPLoVects_PML (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getBfieldCPLoVects_PML (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getBfieldFPLoVects_PML (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getCurrentDensityCPLoVects_PML (int lev, int direction, int *return_size, int **ngrowvect);
  int* warpx_getCurrentDensityFPLoVects_PML (int lev, int direction, int *return_size, int **ngrowvect);

  amrex::Real** warpx_getFfieldCP_PML (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getFfieldFP_PML (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  int* warpx_getFfieldCPLoVects_PML (int lev, int *return_size, int **ngrowvect);
  int* warpx_getFfieldFPLoVects_PML (int lev, int *return_size, int **ngrowvect);
  amrex::Real** warpx_getGfieldCP_PML (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  amrex::Real** warpx_getGfieldFP_PML (int lev, int *return_size, int *ncomps, int **ngrowvect, int **shapes);
  int* warpx_getGfieldCPLoVects_PML (int lev, int *return_size, int **ngrowvect);
  int* warpx_getGfieldFPLoVects_PML (int lev, int *return_size, int **ngrowvect);

  int* warpx_getF_pml_nodal_flag ();
  int* warpx_getG_pml_nodal_flag ();

#ifdef __cplusplus
}
#endif

#endif
