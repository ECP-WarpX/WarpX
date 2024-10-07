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
    : ParticleContainerPureSoA<QdsmcPart::nattribs, 0>(amr_core->GetParGDB())
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
                                       //amrex::Long id)

{
    using namespace amrex::literals;
    using warpx::fields::FieldType;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0, "QdsmcParticleContainer::AddNParticles: only lev=0 is supported yet.");

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(x.size() == n,"x.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(y.size() == n,"y.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(z.size() == n,"z.size() != # of qdsmc particles to add");

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(vx.size() == n,"vx.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(vy.size() == n,"vy.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(vz.size() == n,"vz.size() != # of qdsmc particles to add");

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(entropy.size() == n,"entropy.size() != # of qdsmc particles to add");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(np_real.size() == n,"np_real.size() != # of qdsmc particles to add");

    //  Add to grid 0 and tile 0
    // Redistribute() will move them to proper places.
    auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);

    using PinnedTile = typename ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType;
    PinnedTile pinned_tile;
    pinned_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());

    long ibegin = 0;
    long iend = n;
    const std::size_t np = n;

#ifdef WARPX_DIM_RZ
   amrex::Vector<amrex::ParticleReal> r(np);
    amrex::Vector<amrex::ParticleReal> theta(np);
#endif

    for (auto i = 0; i < np; ++i)
    {
        auto & idcpu_data = pinned_tile.GetStructOfArrays().GetIdCPUData();
        //amrex::Long current_id = id;  // copy input
        //if (id == -1) {
        //    current_id = ParticleType::NextID();
        //}
        current_id = ParticleType::NextID();
        idcpu_data.push_back(amrex::SetParticleIDandCPU(current_id, ParallelDescriptor::MyProc()));

#ifdef WARPX_DIM_RZ
        r[i-ibegin] = std::sqrt(x[i]*x[i] + y[i]*y[i]);
        theta[i-ibegin] = std::atan2(y[i], x[i]);
#endif
    }

    if (np > 0)
    {
#if defined(WARPX_DIM_3D)
        pinned_tile.push_back_real(QdsmcPart::x, x.data() + ibegin, x.data() + iend);
        pinned_tile.push_back_real(QdsmcPart::y, y.data() + ibegin, y.data() + iend);
        pinned_tile.push_back_real(QdsmcPart::z, z.data() + ibegin, z.data() + iend);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        amrex::ignore_unused(y);
#ifdef WARPX_DIM_RZ
        pinned_tile.push_back_real(QdsmcPart::x, r.data(), r.data() + np);
#else
        pinned_tile.push_back_real(QdsmcPart::x, x.data() + ibegin, x.data() + iend);
#endif
        pinned_tile.push_back_real(QdsmcPart::z, z.data() + ibegin, z.data() + iend);
#else //AMREX_SPACEDIM == 1
        amrex::ignore_unused(x,y);
        pinned_tile.push_back_real(QdsmcPart::z, z.data() + ibegin, z.data() + iend);
#endif

        pinned_tile.push_back_real(QdsmcPart::entropy, entropy.data() + ibegin, entropy.data() + iend);
        pinned_tile.push_back_real(QdsmcPart::np_real, np_real.data() + ibegin, np_real.data() + iend);
        pinned_tile.push_back_real(QdsmcPart::vx, vx.data() + ibegin, vx.data() + iend);
        pinned_tile.push_back_real(QdsmcPart::vy, vy.data() + ibegin, vy.data() + iend);
        pinned_tile.push_back_real(QdsmcPart::vz, vz.data() + ibegin, vz.data() + iend);

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ){
            DefineAndReturnParticleTile(0, 0, 0);
        }

        pinned_tile.resize(np);

        auto old_np = particle_tile.numParticles();
        auto new_np = old_np + pinned_tile.numParticles();
        particle_tile.resize(new_np);
        amrex::copyParticles(
            particle_tile, pinned_tile, 0, old_np, pinned_tile.numParticles()
        );
    }

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
