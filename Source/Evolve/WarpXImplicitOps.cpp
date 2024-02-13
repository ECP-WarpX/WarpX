/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "Diagnostics/MultiDiagnostics.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "Evolve/WarpXDtType.H"
#include "Evolve/WarpXPushType.H"
#ifdef WARPX_USE_PSATD
#   ifdef WARPX_DIM_RZ
#       include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#   else
#       include "FieldSolver/SpectralSolver/SpectralSolver.H"
#   endif
#endif
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Python/callbacks.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <vector>

void
WarpX::PreRHSOpFromNonlinearIter ( const amrex::Real  a_cur_time, 
                                   const amrex::Real  a_full_dt, 
                                   const int          a_nl_iter )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_full_dt,a_nl_iter);

    // Advance the particle positions by 1/2 dt,
    // particle velocities by dt, then take average of old and new v,
    // deposit currents, giving J at n+1/2
    // This uses Efield_fp and Bfield_fp, the field at n+1/2 from the previous iteration.
    bool skip_current = false;
    PushType push_type = PushType::Implicit;
    PushParticlesandDeposit(a_cur_time, skip_current, push_type);

    SyncCurrentAndRho();

}

void
WarpX::SaveParticlesAtImplicitStepStart (WarpXParticleContainer& pc, int lev)
{

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {

    auto particle_comps = pc.getParticleComps();

    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {

        const auto getPosition = GetParticlePosition(pti);

        auto& attribs = pti.GetAttribs();
        amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

#if (AMREX_SPACEDIM >= 2)
        amrex::ParticleReal* x_n = pti.GetAttribs(particle_comps["x_n"]).dataPtr();
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
        amrex::ParticleReal* y_n = pti.GetAttribs(particle_comps["y_n"]).dataPtr();
#endif
        amrex::ParticleReal* z_n = pti.GetAttribs(particle_comps["z_n"]).dataPtr();
        amrex::ParticleReal* ux_n = pti.GetAttribs(particle_comps["ux_n"]).dataPtr();
        amrex::ParticleReal* uy_n = pti.GetAttribs(particle_comps["uy_n"]).dataPtr();
        amrex::ParticleReal* uz_n = pti.GetAttribs(particle_comps["uz_n"]).dataPtr();

        const long np = pti.numParticles();

        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            amrex::ParticleReal xp, yp, zp;
            getPosition(ip, xp, yp, zp);

#if (AMREX_SPACEDIM >= 2)
            x_n[ip] = xp;
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
            y_n[ip] = yp;
#endif
            z_n[ip] = zp;

            ux_n[ip] = ux[ip];
            uy_n[ip] = uy[ip];
            uz_n[ip] = uz[ip];

        });

    }
    }
}

void
WarpX::FinishImplicitParticleUpdate (WarpXParticleContainer& pc, int lev)
{
    using namespace amrex::literals;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {

    auto particle_comps = pc.getParticleComps();

    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {

        const auto getPosition = GetParticlePosition(pti);
        const auto setPosition = SetParticlePosition(pti);

        auto& attribs = pti.GetAttribs();
        amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

#if (AMREX_SPACEDIM >= 2)
        amrex::ParticleReal* x_n = pti.GetAttribs(particle_comps["x_n"]).dataPtr();
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
        amrex::ParticleReal* y_n = pti.GetAttribs(particle_comps["y_n"]).dataPtr();
#endif
        amrex::ParticleReal* z_n = pti.GetAttribs(particle_comps["z_n"]).dataPtr();
        amrex::ParticleReal* ux_n = pti.GetAttribs(particle_comps["ux_n"]).dataPtr();
        amrex::ParticleReal* uy_n = pti.GetAttribs(particle_comps["uy_n"]).dataPtr();
        amrex::ParticleReal* uz_n = pti.GetAttribs(particle_comps["uz_n"]).dataPtr();

        const long np = pti.numParticles();

        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            amrex::ParticleReal xp, yp, zp;
            getPosition(ip, xp, yp, zp);

#if (AMREX_SPACEDIM >= 2)
            xp = 2._rt*xp - x_n[ip];
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
            yp = 2._rt*yp - y_n[ip];
#endif
            zp = 2._rt*zp - z_n[ip];

            ux[ip] = 2._rt*ux[ip] - ux_n[ip];
            uy[ip] = 2._rt*uy[ip] - uy_n[ip];
            uz[ip] = 2._rt*uz[ip] - uz_n[ip];

            setPosition(ip, xp, yp, zp);
        });

    }
    }
}

void
WarpX::FinishImplicitFieldUpdate(amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& Field_fp,
                                 amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& Field_n)
{
    using namespace amrex::literals;

    for (int lev = 0; lev <= finest_level; ++lev) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
       for ( amrex::MFIter mfi(*Field_fp[lev][0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {

            amrex::Array4<amrex::Real> const& Fx = Field_fp[lev][0]->array(mfi);
            amrex::Array4<amrex::Real> const& Fy = Field_fp[lev][1]->array(mfi);
            amrex::Array4<amrex::Real> const& Fz = Field_fp[lev][2]->array(mfi);

            amrex::Array4<amrex::Real> const& Fx_n = Field_n[lev][0]->array(mfi);
            amrex::Array4<amrex::Real> const& Fy_n = Field_n[lev][1]->array(mfi);
            amrex::Array4<amrex::Real> const& Fz_n = Field_n[lev][2]->array(mfi);

            amrex::Box tbx = mfi.tilebox(Field_fp[lev][0]->ixType().toIntVect());
            amrex::Box tby = mfi.tilebox(Field_fp[lev][1]->ixType().toIntVect());
            amrex::Box tbz = mfi.tilebox(Field_fp[lev][2]->ixType().toIntVect());

            amrex::ParallelFor(
            tbx, ncomps, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Fx(i,j,k,n) = 2._rt*Fx(i,j,k,n) - Fx_n(i,j,k,n);
            },
            tby, ncomps, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Fy(i,j,k,n) = 2._rt*Fy(i,j,k,n) - Fy_n(i,j,k,n);
            },
            tbz, ncomps, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Fz(i,j,k,n) = 2._rt*Fz(i,j,k,n) - Fz_n(i,j,k,n);
            });
        }
    }
}
