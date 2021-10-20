#include "FieldProbeParticleContainer.H"

#include "Particles/Deposition/ChargeDeposition.H"
#include "Particles/Deposition/CurrentDeposition.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Pusher/UpdatePosition.H"
#include "Particles/ParticleBoundaries_K.H"
#include "Utils/CoarsenMR.H"
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
#include <AMReX_TinyProfiler.H>
#include <AMReX_Utility.H>


#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#include <algorithm>
#include <cmath>

using namespace amrex;

FieldProbeParticleContainer::FieldProbeParticleContainer (AmrCore* amr_core)
	:ParticleContainer<0,0,static_cast<int>(ParticleVal::nattribs)>(amr_core->GetParGDB())
{
	SetParticleSize();
}

void
FieldProbeParticleContainer::AddNParticles (int /*lev*/,
		                    			    int np, const ParticleReal* x, const ParticleReal* y, const ParticleReal* z)
{
	auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);

	using PinnedTile = ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt,
	      				amrex::PinnedArenaAllocator>;
	PinnedTile pinned_tile;
	pinned_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());


	for (int i = 0; i < np; i++)
	{
		ParticleType p;
		p.id() = ParticleType::NextID();
		p.cpu() = ParallelDescriptor::MyProc();
#if (AMREX_SPACEDIM == 3)
		p.pos(0) = x[i];
		p.pos(1) = y[i];
		p.pos(2) = z[i];
#elif (AMREX_SPACEDIM == 2)
		amrex::ignore_unused(y);	
		p.pos(0) = x[i];
		p.pos(1) = 0;
		p.pos(2) = z[i];
#endif
		pinned_tile.push_back(p);

		//Particle Real attributes (SoA) initialized zero
/*		std::array<amrex::Real, 7> real_attribs;
		real_attribs[static_cast<int>(ParticleVal::Ex)] = 0;	//Ex
		real_attribs[static_cast<int>(ParticleVal::Ey)] = 0;	//Ey
		real_attribs[static_cast<int>(ParticleVal::Ez)] = 0;	//Ez
		real_attribs[static_cast<int>(ParticleVal::Bx)] = 0;	//Bx
		real_attribs[static_cast<int>(ParticleVal::By)] = 0;	//By
		real_attribs[static_cast<int>(ParticleVal::Bz)] = 0;	//Bz
		real_attribs[static_cast<int>(ParticleVal::S)] = 0;	//S
*/
        DefineAndReturnParticleTile(0, 0, 0);

		pinned_tile.push_back_real(static_cast<int>(ParticleVal::Ex), np, 0.0);
		pinned_tile.push_back_real(static_cast<int>(ParticleVal::Ey), np, 0.0);
		pinned_tile.push_back_real(static_cast<int>(ParticleVal::Ez), np, 0.0);
		pinned_tile.push_back_real(static_cast<int>(ParticleVal::Bx), np, 0.0);
		pinned_tile.push_back_real(static_cast<int>(ParticleVal::By), np, 0.0);
		pinned_tile.push_back_real(static_cast<int>(ParticleVal::Bz), np, 0.0);
		pinned_tile.push_back_real(static_cast<int>(ParticleVal::S), np, 0.0);
	}
	auto old_np = particle_tile.numParticles();
        auto new_np = old_np + pinned_tile.numParticles();
       	particle_tile.resize(new_np);
        amrex::copyParticles(
		particle_tile, pinned_tile, 0, old_np, pinned_tile.numParticles());
	Redistribute();

}
