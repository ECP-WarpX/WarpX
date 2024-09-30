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
QdsmcParticleContainer::AddNParticles (int /*lev*/, long n,
                                       amrex::Vector<amrex::ParticleReal> const & x,
                                       amrex::Vector<amrex::ParticleReal> const & y,
                                       amrex::Vector<amrex::ParticleReal> const & z,
                                       amrex::Vector<amrex::ParticleReal> const & ux,
                                       amrex::Vector<amrex::ParticleReal> const & uy,
                                       amrex::Vector<amrex::ParticleReal> const & uz,
                                       amrex::Vector<amrex::ParticleReal> const & KeNe,
                                       amrex::Long id)
//                                       const int nattr_real,
//                                       amrex::Vector<amrex::Vector<amrex::ParticleReal>> const & attr_real,
//                                       const int nattr_int,
//                                       amrex::Vector<amrex::Vector<int>> const & attr_int,
//                                       int uniqueparticles, amrex::Long id)
{
    using namespace amrex::literals;
    using warpx::fields::FieldType;

    //WARPX_ALWAYS_ASSERT_WITH_MESSAGE((PIdx::nattribs + nattr_real - 1) <= NumRealComps(),
    //                                 "Too many real attributes specified");
    
    /*
    long ibegin = 0;
    long iend = n;
    
    if (!uniqueparticles) {
        const int myproc = amrex::ParallelDescriptor::MyProc();
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const auto navg = n/nprocs;
        const auto nleft = n - navg * nprocs;
        if (myproc < nleft) {
            ibegin = myproc*(navg+1);
            iend = ibegin + navg+1;
        } else {
            ibegin = myproc*navg + nleft;
            iend = ibegin + navg;
        }
    }*/

    //  Add to grid 0 and tile 0
    // Redistribute() will move them to proper places.
    auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);

    using PinnedTile = typename ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType;
    PinnedTile pinned_tile;
    pinned_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());

    const std::size_t np = n;

#ifdef WARPX_DIM_RZ
    amrex::Vector<amrex::ParticleReal> r(np);
    amrex::Vector<amrex::ParticleReal> theta(np);
#endif

    for (auto i = ibegin; i < iend; ++i)
    {
        auto & idcpu_data = pinned_tile.GetStructOfArrays().GetIdCPUData();

        amrex::Long current_id = id;  // copy input
        if (id == -1) {
            current_id = ParticleType::NextID();
        }
        idcpu_data.push_back(amrex::SetParticleIDandCPU(current_id, ParallelDescriptor::MyProc()));

#ifdef WARPX_DIM_RZ
        r[i-ibegin] = std::sqrt(x[i]*x[i] + y[i]*y[i]);
        theta[i-ibegin] = std::atan2(y[i], x[i]);
#endif
    }

    if (np > 0)
    {
#if defined(WARPX_DIM_3D)
        pinned_tile.push_back_real(PIdx::x, x.data() + ibegin, x.data() + iend);
        pinned_tile.push_back_real(PIdx::y, y.data() + ibegin, y.data() + iend);
        pinned_tile.push_back_real(PIdx::z, z.data() + ibegin, z.data() + iend);
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        amrex::ignore_unused(y);
#ifdef WARPX_DIM_RZ
        pinned_tile.push_back_real(PIdx::x, r.data(), r.data() + np);
#else
        pinned_tile.push_back_real(PIdx::x, x.data() + ibegin, x.data() + iend);
#endif
        pinned_tile.push_back_real(PIdx::z, z.data() + ibegin, z.data() + iend);
#else //AMREX_SPACEDIM == 1
        amrex::ignore_unused(x,y);
        pinned_tile.push_back_real(PIdx::z, z.data() + ibegin, z.data() + iend);
#endif

        pinned_tile.push_back_real(PIdx::w, attr_real[0].data() + ibegin, attr_real[0].data() + iend);
        pinned_tile.push_back_real(PIdx::ux, ux.data() + ibegin, ux.data() + iend);
        pinned_tile.push_back_real(PIdx::uy, uy.data() + ibegin, uy.data() + iend);
        pinned_tile.push_back_real(PIdx::uz, uz.data() + ibegin, uz.data() + iend);

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ){
            DefineAndReturnParticleTile(0, 0, 0);
        }

        for (int comp = PIdx::uz+1; comp < PIdx::nattribs; ++comp)
        {
#ifdef WARPX_DIM_RZ
            if (comp == PIdx::theta) {
                pinned_tile.push_back_real(comp, theta.data(), theta.data() + np);
            }
            else {
                pinned_tile.push_back_real(comp, np, 0.0_prt);
            }
#else
            pinned_tile.push_back_real(comp, np, 0.0_prt);
#endif
        }

        // Initialize nattr_real - 1 runtime real attributes from data in the attr_real array
        for (int j = PIdx::nattribs; j < PIdx::nattribs + nattr_real - 1; ++j)
        {
            // get the next attribute from attr_real array
            pinned_tile.push_back_real(
                j, attr_real[j - PIdx::nattribs + 1].data() + ibegin, attr_real[j - PIdx::nattribs + 1].data() + iend
            );
        }

        // Initialize nattr_int runtime integer attributes from data in the attr_int array
        for (int j = 0; j < nattr_int; ++j)
        {
            // get the next attribute from attr_int array
            pinned_tile.push_back_int(j, attr_int[j].data() + ibegin, attr_int[j].data() + iend);
        }

        pinned_tile.resize(np);
        // Default initialize the other real and integer runtime attributes
        DefaultInitializeRuntimeAttributes(pinned_tile, nattr_real - 1, nattr_int);

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
}
