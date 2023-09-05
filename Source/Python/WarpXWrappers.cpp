/* Copyright 2019 Andrew Myers, Axel Huebl, David Grote
 * Luca Fedeli, Maxence Thevenet, Remi Lehe
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BoundaryConditions/PML.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Initialization/WarpXAMReXInit.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"
#include "WarpXWrappers.H"
#include "WarpX_py.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#include <AMReX_StructOfArrays.H>

#include <array>
#include <cstdlib>

    int warpx_Real_size()
    {
        return (int)sizeof(amrex::Real);
    }

    int warpx_ParticleReal_size()
    {
        return (int)sizeof(amrex::ParticleReal);
    }

    int warpx_nSpecies()
    {
        const auto & mypc = WarpX::GetInstance().GetPartContainer();
        return mypc.nSpecies();
    }

    bool warpx_use_fdtd_nci_corr()
    {
        return WarpX::use_fdtd_nci_corr;
    }

    int warpx_galerkin_interpolation()
    {
        return WarpX::galerkin_interpolation;
    }

    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm /* mpicomm */)
    {
      warpx::initialization::amrex_init(argc, argv, true);
    }

    void warpx_ConvertLabParamsToBoost()
    {
      ConvertLabParamsToBoost();
    }

    void warpx_ReadBCParams()
    {
      ReadBCParams();
    }

    void warpx_CheckGriddingForRZSpectral()
    {
        CheckGriddingForRZSpectral();
    }

    amrex::Real warpx_getCellSize(int dir, int lev) {
        const std::array<amrex::Real,3>& dx = WarpX::CellSize(lev);
        return dx[dir];
    }

    amrex::Real warpx_sumParticleCharge(const char* char_species_name, const bool local)
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);
        return myspc.sumParticleCharge(local);
    }

    void warpx_ComputeDt () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.ComputeDt();
    }
    void warpx_MoveWindow (int step,bool move_j) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.MoveWindow(step, move_j);
    }

    void warpx_EvolveE (amrex::Real dt) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveE(dt);
    }
    void warpx_EvolveB (amrex::Real dt, DtType a_dt_type) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveB(dt, a_dt_type);
    }
    void warpx_FillBoundaryE () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryE(warpx.getngEB());
    }
    void warpx_FillBoundaryB () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryB(warpx.getngEB());
    }
    void warpx_SyncCurrent (
        const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
        const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp,
        const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_buffer) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.SyncCurrent(J_fp, J_cp, J_buffer);
    }
    void warpx_UpdateAuxilaryData () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.UpdateAuxilaryData();
    }
    void warpx_PushParticlesandDepose (amrex::Real cur_time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.PushParticlesandDepose(cur_time);
    }

    void warpx_setistep (int lev, int ii) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.setistep(lev, ii);
    }
    void warpx_sett_new (int lev, amrex::Real time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.sett_new(lev, time);
    }
    amrex::Real warpx_getdt (int lev) {
        const WarpX& warpx = WarpX::GetInstance();
        return warpx.getdt(lev);
    }

    int warpx_maxStep () {
        const WarpX& warpx = WarpX::GetInstance();
        return warpx.maxStep();
    }
    amrex::Real warpx_stopTime () {
        const WarpX& warpx = WarpX::GetInstance();
        return warpx.stopTime();
    }

    int warpx_finestLevel () {
        const WarpX& warpx = WarpX::GetInstance();
        return warpx.finestLevel();
    }

    void mypc_Redistribute () {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        mypc.Redistribute();
    }
