/* Copyright 2024 Marco Acciarri (Helion Energy Inc.)
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "QdsmcParticleContainer.H"

#include "Particles/Deposition/ChargeDeposition.H"
#include "Particles/Deposition/CurrentDeposition.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Pusher/UpdatePosition.H"
#include "Particles/ParticleBoundaries_K.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_Dim3.H>
#include <AMReX_Extension.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAllocators.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleContainerBase.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParticleTransformation.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_Utility.H>


#include <algorithm>
#include <cmath>

using namespace amrex;

QdsmcParticleContainer::QdsmcParticleContainer (AmrCore* amr_core)
    : ParticleContainerPureSoA<QdsmcPIdx::nattribs, 0>(amr_core->GetParGDB())
{
    SetParticleSize();
}

void
QdsmcParticleContainer::AddNParticles (int lev, long n,
                                       amrex::Vector<amrex::ParticleReal> const & x,
                                       amrex::Vector<amrex::ParticleReal> const & y,
                                       amrex::Vector<amrex::ParticleReal> const & z,
                                       amrex::Vector<amrex::ParticleReal> const & vx,
                                       amrex::Vector<amrex::ParticleReal> const & vy,
                                       amrex::Vector<amrex::ParticleReal> const & vz,
                                       amrex::Vector<amrex::ParticleReal> const & entropy,
                                       amrex::Vector<amrex::ParticleReal> const & np_real)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0, "QdsmcParticleContainer::AddNParticles: only lev=0 is supported yet.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(x.size() == n,"x.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(y.size() == n,"y.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(z.size() == n,"z.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(vx.size() == n,"vx.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(vy.size() == n,"vy.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(vz.size() == n,"vz.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(entropy.size() == n,"entropy.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(np_real.size() == n,"np_real.size() != # of qdsmc particles to add");

    if (n <= 0){
        Redistribute();
        return;
    }

    // have to resize here, not in the constructor because grids have not
    // been built when constructor was called.
    reserveData();
    resizeData();

    /*
     * Creates a temporary tile to obtain data from simulation. This data
     * is then coppied to the permament tile which is stored on the particle
     * (particle_tile).
     */
    using PinnedTile = typename ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType;

    PinnedTile pinned_tile;
    pinned_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());

    for (int i = 0; i < n; i++)
    {
        auto & idcpu_data = pinned_tile.GetStructOfArrays().GetIdCPUData();
        idcpu_data.push_back(amrex::SetParticleIDandCPU(ParticleType::NextID(), ParallelDescriptor::MyProc()));
    }

    // write Real attributes (SoA) to particle initialized zero
    DefineAndReturnParticleTile(0, 0, 0);

    // for RZ write theta value
#ifdef WARPX_DIM_RZ
    pinned_tile.push_back_real(QdsmcPIdx::theta, np, 0.0);
#endif
#if !defined (WARPX_DIM_1D_Z)
    pinned_tile.push_back_real(QdsmcPIdx::x, x);
#endif
#if defined (WARPX_DIM_3D)
    pinned_tile.push_back_real(QdsmcPIdx::y, y);
#endif
    pinned_tile.push_back_real(QdsmcPIdx::z, z);
    pinned_tile.push_back_real(QdsmcPIdx::vx, n, 0.0);
    pinned_tile.push_back_real(QdsmcPIdx::vy, n, 0.0);
    pinned_tile.push_back_real(QdsmcPIdx::vz, n, 0.0);
    pinned_tile.push_back_real(QdsmcPIdx::entropy, n, 0.0);
    pinned_tile.push_back_real(QdsmcPIdx::np_real, n, 0.0);

    // Move particles to their appropriate tiles
    Redistribute();

    // Remove particles that are inside the embedded boundaries
/*
#ifdef AMREX_USE_EB
    if (EB::enabled()) {
        auto & warpx = WarpX::GetInstance();
        scrapeParticlesAtEB(
            *this,
            warpx.m_fields.get_mr_levels(FieldType::distance_to_eb, warpx.finestLevel()),
            ParticleBoundaryProcess::Absorb());
        deleteInvalidParticles();
    }
#endif
*/
}
