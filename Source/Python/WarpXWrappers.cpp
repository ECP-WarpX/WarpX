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

    int warpx_nComps()
    {
        return PIdx::nattribs;
    }

    int warpx_SpaceDim()
    {
        return AMREX_SPACEDIM;
    }

    void amrex_init (int argc, char* argv[])
    {
        warpx::initialization::amrex_init(argc, argv);
    }

    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm /* mpicomm */)
    {
      warpx::initialization::amrex_init(argc, argv, true);
    }

    void amrex_finalize (int /*finalize_mpi*/)
    {
        amrex::Finalize();
    }

    void warpx_init ()
    {
        WarpX& warpx = WarpX::GetInstance();
        warpx.InitData();
        ExecutePythonCallback("afterinit");
        ExecutePythonCallback("particleloader");
    }

    void warpx_finalize ()
    {
        WarpX::ResetInstance();
    }

    void warpx_set_callback_py (
        const char* char_callback_name, WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        const std::string callback_name(char_callback_name);
        warpx_callback_py_map[callback_name] = callback;
    }

    void warpx_clear_callback_py (const char* char_callback_name)
    {
        const std::string callback_name(char_callback_name);
        warpx_callback_py_map.erase(callback_name);
    }

    void warpx_evolve (int numsteps)
    {
        WarpX& warpx = WarpX::GetInstance();
        warpx.Evolve(numsteps);
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

    amrex::Real warpx_getProbLo(int dir)
    {
      WarpX& warpx = WarpX::GetInstance();
      const amrex::Geometry& geom = warpx.Geom(0);
      return geom.ProbLo(dir);
    }

    amrex::Real warpx_getProbHi(int dir)
    {
      WarpX& warpx = WarpX::GetInstance();
      const amrex::Geometry& geom = warpx.Geom(0);
      return geom.ProbHi(dir);
    }

    amrex::Real warpx_getCellSize(int dir, int lev) {
        const std::array<amrex::Real,3>& dx = WarpX::CellSize(lev);
        return dx[dir];
    }

    amrex::ParticleReal** warpx_getParticleStructs(
            const char* char_species_name, int lev,
            int* num_tiles, int** particles_per_tile) {
        const auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);

        *num_tiles = myspc.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(typename WarpXParticleContainer::ParticleType*)));
        int i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {
            auto& aos = pti.GetArrayOfStructs();
            data[i] = (amrex::ParticleReal*) aos.data();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    amrex::ParticleReal** warpx_getParticleArrays (
            const char* char_species_name, const char* char_comp_name,
            int lev, int* num_tiles, int** particles_per_tile ) {

        const auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);

        const int comp = warpx_getParticleCompIndex(char_species_name, char_comp_name);

        *num_tiles = myspc.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(amrex::ParticleReal*)));
        int i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (amrex::ParticleReal*) soa.GetRealData(comp).dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    void warpx_convert_id_to_long (amrex::Long* ids, const WarpXParticleContainer::ParticleType* pstructs, int size)
    {
        amrex::Long* d_ptr = nullptr;
#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<amrex::Long> d_ids(size);
        d_ptr = d_ids.data();
#else
        d_ptr = ids;
#endif
        amrex::ParallelFor(size, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            d_ptr[i] = pstructs[i].id();
        });
#ifdef AMREX_USE_GPU
        amrex::Gpu::dtoh_memcpy(ids, d_ptr, size*sizeof(amrex::Long));
#endif
    }

    void warpx_convert_cpu_to_int (int* cpus, const WarpXParticleContainer::ParticleType* pstructs, int size)
    {
        int* d_ptr = nullptr;
#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<int> d_cpus(size);
        d_ptr = d_cpus.data();
#else
        d_ptr = cpus;
#endif
        amrex::ParallelFor(size, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            d_ptr[i] = pstructs[i].cpu();
        });
#ifdef AMREX_USE_GPU
        amrex::Gpu::dtoh_memcpy(cpus, d_ptr, size*sizeof(int));
#endif
    }

    int warpx_getParticleCompIndex (
         const char* char_species_name, const char* char_comp_name )
    {
        const auto & mypc = WarpX::GetInstance().GetPartContainer();

        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);

        const std::string comp_name(char_comp_name);
        auto particle_comps = myspc.getParticleComps();

        return particle_comps.at(comp_name);
    }

    amrex::Real warpx_sumParticleCharge(const char* char_species_name, const bool local)
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);
        return myspc.sumParticleCharge(local);
    }

    int** warpx_getParticleBoundaryBufferScrapedSteps(const char* species_name, int boundary, int lev,
                     int* num_tiles, int** particles_per_tile)
    {
        const std::string name(species_name);
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        auto& particle_buffer = particle_buffers.getParticleBuffer(species_name, boundary);

        const int comp = particle_buffer.NumIntComps() - 1;

        *num_tiles = particle_buffer.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<int**>(malloc(*num_tiles*sizeof(int*)));
        int i = 0;
        for (amrex::ParIter<0,0,PIdx::nattribs, 0, amrex::PinnedArenaAllocator> pti(particle_buffer, lev); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (int*) soa.GetIntData(comp).dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }

        return data;
    }

    amrex::ParticleReal** warpx_getParticleBoundaryBuffer(const char* species_name, int boundary, int lev,
                     int* num_tiles, int** particles_per_tile, const char* comp_name)
    {
        const std::string name(species_name);
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        auto& particle_buffer = particle_buffers.getParticleBuffer(species_name, boundary);

        const int comp = warpx_getParticleCompIndex(species_name, comp_name);

        *num_tiles = particle_buffer.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(amrex::ParticleReal*)));
        int i = 0;
        for (amrex::ParIter<0,0,PIdx::nattribs, 0, amrex::PinnedArenaAllocator> pti(particle_buffer, lev); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (amrex::ParticleReal*) soa.GetRealData(comp).dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }

        return data;
    }

    amrex::ParticleReal** warpx_getParticleBoundaryBufferStructs(const char* species_name, int boundary, int lev,
                     int* num_tiles, int** particles_per_tile)
    {
        const std::string name(species_name);
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        auto& particle_buffer = particle_buffers.getParticleBuffer(species_name, boundary);

        *num_tiles = particle_buffer.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(typename WarpXParticleContainer::ParticleType*)));
        int i = 0;
        for (amrex::ParIter<0,0,PIdx::nattribs, 0, amrex::PinnedArenaAllocator> pti(particle_buffer, lev); pti.isValid(); ++pti, ++i) {
            auto& aos = pti.GetArrayOfStructs();
            data[i] = (amrex::ParticleReal*) aos.data();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    void warpx_clearParticleBoundaryBuffer () {
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        particle_buffers.clearParticles();
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

    void warpx_setPotentialEB (const char * char_potential) {
        WarpX& warpx = WarpX::GetInstance();
        const std::string potential(char_potential);
        warpx.m_poisson_boundary_handler.setPotentialEB(potential);
    }

    void mypc_Redistribute () {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        mypc.Redistribute();
    }
