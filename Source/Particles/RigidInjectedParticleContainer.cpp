/* Copyright 2019-2020 Andrew Myers, David Grote, Glenn Richardson
 * Ligia Diana Amorim, Luca Fedeli, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifdef _OPENMP
#   include <omp.h>
#endif

#include "RigidInjectedParticleContainer.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Pusher/UpdateMomentumBoris.H"
#include "Pusher/UpdateMomentumVay.H"
#include "Pusher/UpdateMomentumBorisWithRadiationReaction.H"
#include "Pusher/UpdateMomentumHigueraCary.H"
#include "Pusher/GetAndSetPosition.H"
#include "Gather/ScaleFields.H"
#include "Gather/FieldGather.H"

#include <limits>
#include <sstream>
#include <algorithm>


using namespace amrex;

RigidInjectedParticleContainer::RigidInjectedParticleContainer (AmrCore* amr_core, int ispecies,
                                                                const std::string& name)
    : PhysicalParticleContainer(amr_core, ispecies, name)
{

    ParmParse pp(species_name);

    pp.get("zinject_plane", zinject_plane);
    pp.query("projected", projected);
    pp.query("focused", focused);
    pp.query("rigid_advance", rigid_advance);

}

void RigidInjectedParticleContainer::InitData()
{
    done_injecting.resize(finestLevel()+1, 0);
    zinject_plane_levels.resize(finestLevel()+1, zinject_plane/WarpX::gamma_boost);

    AddParticles(0); // Note - add on level 0

    // Particles added by AddParticles should already be in the boosted frame
    RemapParticles();

    Redistribute();  // We then redistribute
}

void
RigidInjectedParticleContainer::RemapParticles()
{

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(projected, "ERROR: projected = false is not supported with this particle loading");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!focused, "ERROR: focused = true is not supported with this particle loading");

    // For rigid_advance == false, nothing needs to be done

    if (rigid_advance) {

        // The particle z positions are adjusted to account for the difference between
        // advancing with vzbar and wih vz[i] before injection

        // For now, start with the assumption that this will only happen
        // at the start of the simulation.
        const Real t_lab = 0.;

        const Real uz_boost = WarpX::gamma_boost*WarpX::beta_boost*PhysConst::c;
        const Real csqi = 1./(PhysConst::c*PhysConst::c);

        vzbeam_ave_boosted = meanParticleVelocity(false)[2];

        for (int lev = 0; lev <= finestLevel(); lev++) {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                // Get the average beam velocity in the boosted frame.
                // Note that the particles are already in the boosted frame.
                // This value is saved to advance the particles not injected yet

                for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
                {
                    const auto& attribs = pti.GetAttribs();
                    const auto uxp = attribs[PIdx::ux].dataPtr();
                    const auto uyp = attribs[PIdx::uy].dataPtr();
                    const auto uzp = attribs[PIdx::uz].dataPtr();

                    const auto GetPosition = GetParticlePosition(pti);
                          auto SetPosition = SetParticlePosition(pti);

                    // Loop over particles
                    const long np = pti.numParticles();
                    const Real lvzbeam_ave_boosted = vzbeam_ave_boosted;
                    const Real gamma_boost = WarpX::gamma_boost;
                    amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long i)
                    {
                        ParticleReal xp, yp, zp;
                        GetPosition(i, xp, yp, zp);

                        const Real gammapr = std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])*csqi);
                        const Real vzpr = uzp[i]/gammapr;

                        // Back out the value of z_lab
                        const Real z_lab = (zp + uz_boost*t_lab + gamma_boost*t_lab*vzpr)/(gamma_boost + uz_boost*vzpr*csqi);

                        // Time of the particle in the boosted frame given its position in the lab frame at t=0.
                        const Real tpr = gamma_boost*t_lab - uz_boost*z_lab*csqi;

                        // Adjust the position, taking away its motion from its own velocity and adding
                        // the motion from the average velocity
                        zp += tpr*vzpr - tpr*lvzbeam_ave_boosted;

                        SetPosition(i, xp, yp, zp);
                    });
                }
            }
        }
    }
}

void
RigidInjectedParticleContainer::BoostandRemapParticles()
{

    // Boost the particles into the boosted frame and map the particles
    // to the t=0 in the boosted frame. If using rigid_advance, the z position
    // is adjusted using vzbar, otherwise using vz[i]

    if (rigid_advance) {
        // Get the average beam velocity in the boosted frame
        // This value is saved to advance the particles not injected yet
        const Real vzbeam_ave_lab = meanParticleVelocity(false)[2];
        vzbeam_ave_boosted = (vzbeam_ave_lab - WarpX::beta_boost*PhysConst::c)/(1. - vzbeam_ave_lab*WarpX::beta_boost/PhysConst::c);
    }

    const Real csqi = 1./(PhysConst::c*PhysConst::c);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (WarpXParIter pti(*this, 0); pti.isValid(); ++pti)
        {

            auto& attribs = pti.GetAttribs();
            const auto uxp = attribs[PIdx::ux].dataPtr();
            const auto uyp = attribs[PIdx::uy].dataPtr();
            const auto uzp = attribs[PIdx::uz].dataPtr();

            const auto GetPosition = GetParticlePosition(pti);
                  auto SetPosition = SetParticlePosition(pti);

            // Loop over particles
            const long np = pti.numParticles();
            const Real lvzbeam_ave_boosted = vzbeam_ave_boosted;
            const Real gamma_boost = WarpX::gamma_boost;
            const Real beta_boost = WarpX::beta_boost;
            const bool lprojected = projected;
            const bool lfocused = focused;
            const bool lrigid_advance = rigid_advance;
            const Real lzinject_plane = zinject_plane;
            amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long i)
            {
                ParticleReal xp, yp, zp;
                GetPosition(i, xp, yp, zp);

                const Real gamma_lab = std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])*csqi);

                const Real vx_lab = uxp[i]/gamma_lab;
                const Real vy_lab = uyp[i]/gamma_lab;
                const Real vz_lab = uzp[i]/gamma_lab;

                // t0_lab is the time in the lab frame that the particles reaches z=0
                // The location and time (z=0, t=0) is a synchronization point between the
                // lab and boosted frames.
                const Real t0_lab = -zp/vz_lab;

                if (!lprojected) {
                    xp += t0_lab*vx_lab;
                    yp += t0_lab*vy_lab;
                }
                if (lfocused) {
                    // Correct for focusing effect from shift from z=0 to zinject
                    const Real tfocus = -lzinject_plane*gamma_boost/vz_lab;
                    xp -= tfocus*vx_lab;
                    yp -= tfocus*vy_lab;
                }

                // Time of the particle in the boosted frame given its position in the lab frame at t=0.
                const Real tpr = -gamma_boost*beta_boost*zp/PhysConst::c;

                // Position of the particle in the boosted frame given its position in the lab frame at t=0.
                const Real zpr = gamma_boost*zp;

                // Momentum of the particle in the boosted frame (assuming that it is fixed).
                uzp[i] = gamma_boost*(uzp[i] - beta_boost*PhysConst::c*gamma_lab);

                // Put the particle at the location in the boosted frame at boost frame t=0,
                if (lrigid_advance) {
                    // with the particle moving at the average velocity
                    zp = zpr - lvzbeam_ave_boosted*tpr;
                }
                else {
                    // with the particle moving with its own velocity
                    const Real gammapr = std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])*csqi);
                    const Real vzpr = uzp[i]/gammapr;
                    zp = zpr - vzpr*tpr;
                }

                SetPosition(i, xp, yp, zp);
            });
        }
    }
}

void
RigidInjectedParticleContainer::PushPX (WarpXParIter& pti,
                                        amrex::FArrayBox const * exfab,
                                        amrex::FArrayBox const * eyfab,
                                        amrex::FArrayBox const * ezfab,
                                        amrex::FArrayBox const * bxfab,
                                        amrex::FArrayBox const * byfab,
                                        amrex::FArrayBox const * bzfab,
                                        const int ngE, const int e_is_nodal,
                                        const long offset,
                                        const long np_to_push,
                                        int lev, int gather_lev,
                                        amrex::Real dt, ScaleFields /*scaleFields*/,
                                        DtType a_dt_type)
{
    auto& attribs = pti.GetAttribs();
    auto& uxp = attribs[PIdx::ux];
    auto& uyp = attribs[PIdx::uy];
    auto& uzp = attribs[PIdx::uz];

    // Save the position and momenta, making copies
    Gpu::DeviceVector<ParticleReal> xp_save, yp_save, zp_save;
    RealVector uxp_save, uyp_save, uzp_save;

    const auto GetPosition = GetParticlePosition(pti);
          auto SetPosition = SetParticlePosition(pti);

    ParticleReal* const AMREX_RESTRICT ux = uxp.dataPtr();
    ParticleReal* const AMREX_RESTRICT uy = uyp.dataPtr();
    ParticleReal* const AMREX_RESTRICT uz = uzp.dataPtr();

    if (!done_injecting_lev)
    {
        // If the old values are not already saved, create copies here.
        const auto np = pti.numParticles();

        xp_save.resize(np);
        yp_save.resize(np);
        zp_save.resize(np);

        uxp_save.resize(np);
        uyp_save.resize(np);
        uzp_save.resize(np);

        amrex::Real* const AMREX_RESTRICT xp_save_ptr = xp_save.dataPtr();
        amrex::Real* const AMREX_RESTRICT yp_save_ptr = yp_save.dataPtr();
        amrex::Real* const AMREX_RESTRICT zp_save_ptr = zp_save.dataPtr();

        amrex::Real* const AMREX_RESTRICT uxp_save_ptr = uxp_save.dataPtr();
        amrex::Real* const AMREX_RESTRICT uyp_save_ptr = uyp_save.dataPtr();
        amrex::Real* const AMREX_RESTRICT uzp_save_ptr = uzp_save.dataPtr();

        amrex::ParallelFor( np,
                            [=] AMREX_GPU_DEVICE (long i) {
                                ParticleReal xp, yp, zp;
                                GetPosition(i, xp, yp, zp);
                                xp_save_ptr[i] = xp;
                                yp_save_ptr[i] = yp;
                                zp_save_ptr[i] = zp;
                                uxp_save_ptr[i] = ux[i];
                                uyp_save_ptr[i] = uy[i];
                                uzp_save_ptr[i] = uz[i];
                            });
    }

    const bool do_scale = not done_injecting_lev;
    const Real v_boost = WarpX::beta_boost*PhysConst::c;
    PhysicalParticleContainer::PushPX(pti, exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                                      ngE, e_is_nodal, offset, np_to_push, lev, gather_lev, dt,
                                      ScaleFields(do_scale, dt, zinject_plane_lev_previous,
                                                  vzbeam_ave_boosted, v_boost),
                                      a_dt_type);

    if (!done_injecting_lev) {

        ParticleReal* AMREX_RESTRICT x_save = xp_save.dataPtr();
        ParticleReal* AMREX_RESTRICT y_save = yp_save.dataPtr();
        ParticleReal* AMREX_RESTRICT z_save = zp_save.dataPtr();
        ParticleReal* AMREX_RESTRICT ux_save = uxp_save.dataPtr();
        ParticleReal* AMREX_RESTRICT uy_save = uyp_save.dataPtr();
        ParticleReal* AMREX_RESTRICT uz_save = uzp_save.dataPtr();

        // Undo the push for particles not injected yet.
        // The zp are advanced a fixed amount.
        const Real z_plane_lev = zinject_plane_lev;
        const Real vz_ave_boosted = vzbeam_ave_boosted;
        const bool rigid = rigid_advance;
        const Real inv_csq = 1./(PhysConst::c*PhysConst::c);
        amrex::ParallelFor( pti.numParticles(),
                            [=] AMREX_GPU_DEVICE (long i) {
                                ParticleReal xp, yp, zp;
                                GetPosition(i, xp, yp, zp);
                                if (zp <= z_plane_lev) {
                                    ux[i] = ux_save[i];
                                    uy[i] = uy_save[i];
                                    uz[i] = uz_save[i];
                                    xp = x_save[i];
                                    yp = y_save[i];
                                    if (rigid) {
                                        zp = z_save[i] + dt*vz_ave_boosted;
                                    }
                                    else {
                                        const Real gi = 1./std::sqrt(1. + (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i])*inv_csq);
                                        zp = z_save[i] + dt*uz[i]*gi;
                                    }
                                    SetPosition(i, xp, yp, zp);
                                }
                            });
    }
}

void
RigidInjectedParticleContainer::Evolve (int lev,
                                        const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                        const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                        const MultiFab& Ex_avg, const MultiFab& Ey_avg, const MultiFab& Ez_avg,
                                        const MultiFab& Bx_avg, const MultiFab& By_avg, const MultiFab& Bz_avg,
                                        MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                        MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                        MultiFab* rho, MultiFab* crho,
                                        const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                        const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                        Real t, Real dt, DtType a_dt_type)
{

    // Update location of injection plane in the boosted frame
    zinject_plane_lev_previous = zinject_plane_levels[lev];
    zinject_plane_levels[lev] -= dt*WarpX::beta_boost*PhysConst::c;
    zinject_plane_lev = zinject_plane_levels[lev];

    // Set the done injecting flag whan the inject plane moves out of the
    // simulation domain.
    // It is much easier to do this check, rather than checking if all of the
    // particles have crossed the inject plane.
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();
    const int zdir = AMREX_SPACEDIM-1;
    done_injecting[lev] = ((zinject_plane_levels[lev] < plo[zdir] && WarpX::moving_window_v + WarpX::beta_boost*PhysConst::c >= 0.) ||
                           (zinject_plane_levels[lev] > phi[zdir] && WarpX::moving_window_v + WarpX::beta_boost*PhysConst::c <= 0.));
    done_injecting_lev = done_injecting[lev];

    PhysicalParticleContainer::Evolve (lev,
                                       Ex, Ey, Ez,
                                       Bx, By, Bz,
                                       Ex_avg, Ey_avg, Ez_avg,
                                       Bx_avg, By_avg, Bz_avg,
                                       jx, jy, jz,
                                       cjx, cjy, cjz,
                                       rho, crho,
                                       cEx, cEy, cEz,
                                       cBx, cBy, cBz,
                                       t, dt, a_dt_type);
}

void
RigidInjectedParticleContainer::PushP (int lev, Real dt,
                                       const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                       const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    WARPX_PROFILE("RigidInjectedParticleContainer::PushP");

    if (do_not_push) return;

    const std::array<Real,3>& dx = WarpX::CellSize(std::max(lev,0));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            amrex::Box box = pti.tilebox();
            box.grow(Ex.nGrow());

            const long np = pti.numParticles();

            // Data on the grid
            const FArrayBox& exfab = Ex[pti];
            const FArrayBox& eyfab = Ey[pti];
            const FArrayBox& ezfab = Ez[pti];
            const FArrayBox& bxfab = Bx[pti];
            const FArrayBox& byfab = By[pti];
            const FArrayBox& bzfab = Bz[pti];

            const auto getPosition = GetParticlePosition(pti);

            const auto getExternalE = GetExternalEField(pti);
            const auto getExternalB = GetExternalBField(pti);

            const auto& xyzmin = WarpX::GetInstance().LowerCornerWithGalilean(box,m_v_galilean,lev);

            const Dim3 lo = lbound(box);

            bool galerkin_interpolation = WarpX::galerkin_interpolation;
            int nox = WarpX::nox;
            int n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;

            amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
            amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};

            amrex::Array4<const amrex::Real> const& ex_arr = exfab.array();
            amrex::Array4<const amrex::Real> const& ey_arr = eyfab.array();
            amrex::Array4<const amrex::Real> const& ez_arr = ezfab.array();
            amrex::Array4<const amrex::Real> const& bx_arr = bxfab.array();
            amrex::Array4<const amrex::Real> const& by_arr = byfab.array();
            amrex::Array4<const amrex::Real> const& bz_arr = bzfab.array();

            amrex::IndexType const ex_type = exfab.box().ixType();
            amrex::IndexType const ey_type = eyfab.box().ixType();
            amrex::IndexType const ez_type = ezfab.box().ixType();
            amrex::IndexType const bx_type = bxfab.box().ixType();
            amrex::IndexType const by_type = byfab.box().ixType();
            amrex::IndexType const bz_type = bzfab.box().ixType();

            auto& attribs = pti.GetAttribs();
            amrex::ParticleReal* const AMREX_RESTRICT uxpp = attribs[PIdx::ux].dataPtr();
            amrex::ParticleReal* const AMREX_RESTRICT uypp = attribs[PIdx::uy].dataPtr();
            amrex::ParticleReal* const AMREX_RESTRICT uzpp = attribs[PIdx::uz].dataPtr();

            int* AMREX_RESTRICT ion_lev = nullptr;
            if (do_field_ionization) {
                ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
            }

            // Save the position and momenta, making copies
            amrex::Gpu::DeviceVector<ParticleReal> uxp_save(np);
            amrex::Gpu::DeviceVector<ParticleReal> uyp_save(np);
            amrex::Gpu::DeviceVector<ParticleReal> uzp_save(np);
            ParticleReal* const AMREX_RESTRICT ux_save = uxp_save.dataPtr();
            ParticleReal* const AMREX_RESTRICT uy_save = uyp_save.dataPtr();
            ParticleReal* const AMREX_RESTRICT uz_save = uzp_save.dataPtr();

            // Loop over the particles and update their momentum
            const amrex::Real q = this->charge;
            const amrex::Real m = this-> mass;

            const auto pusher_algo = WarpX::particle_pusher_algo;
            const auto do_crr = do_classical_radiation_reaction;

            amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
            {
                ux_save[ip] = uxpp[ip];
                uy_save[ip] = uypp[ip];
                uz_save[ip] = uzpp[ip];

                amrex::ParticleReal xp, yp, zp;
                getPosition(ip, xp, yp, zp);

                amrex::ParticleReal Exp = 0._rt, Eyp = 0._rt, Ezp = 0._rt;
                amrex::ParticleReal Bxp = 0._rt, Byp = 0._rt, Bzp = 0._rt;

                // first gather E and B to the particle positions
                doGatherShapeN(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                               ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                               ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                               dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes,
                               nox, galerkin_interpolation);
                getExternalE(ip, Exp, Eyp, Ezp);
                getExternalB(ip, Bxp, Byp, Bzp);

                if (do_crr) {
                    amrex::Real qp = q;
                    if (ion_lev) { qp *= ion_lev[ip]; }
                    UpdateMomentumBorisWithRadiationReaction(uxpp[ip], uypp[ip], uzpp[ip],
                                                             Exp, Eyp, Ezp, Bxp,
                                                             Byp, Bzp, qp, m, dt);
                } else if (pusher_algo == ParticlePusherAlgo::Boris) {
                    amrex::Real qp = q;
                    if (ion_lev) { qp *= ion_lev[ip]; }
                    UpdateMomentumBoris( uxpp[ip], uypp[ip], uzpp[ip],
                                         Exp, Eyp, Ezp, Bxp,
                                         Byp, Bzp, qp, m, dt);
                } else if (pusher_algo == ParticlePusherAlgo::Vay) {
                    amrex::Real qp = q;
                    if (ion_lev){ qp *= ion_lev[ip]; }
                    UpdateMomentumVay( uxpp[ip], uypp[ip], uzpp[ip],
                                       Exp, Eyp, Ezp, Bxp,
                                       Byp, Bzp, qp, m, dt);
                } else if (pusher_algo == ParticlePusherAlgo::HigueraCary) {
                    amrex::Real qp = q;
                    if (ion_lev){ qp *= ion_lev[ip]; }
                    UpdateMomentumHigueraCary( uxpp[ip], uypp[ip], uzpp[ip],
                                               Exp, Eyp, Ezp, Bxp,
                                               Byp, Bzp, qp, m, dt);
                } else {
                    amrex::Abort("Unknown particle pusher");
                }
            });

            // Undo the push for particles not injected yet.
            // It is assumed that PushP will only be called on the first and last steps
            // and that no particles will cross zinject_plane.
            const ParticleReal zz = zinject_plane_levels[lev];
            amrex::ParallelFor( pti.numParticles(), [=] AMREX_GPU_DEVICE (long i)
            {
                ParticleReal xp, yp, zp;
                getPosition(i, xp, yp, zp);
                if (zp <= zz) {
                    uxpp[i] = ux_save[i];
                    uypp[i] = uy_save[i];
                    uzpp[i] = uz_save[i];
                }
            });

            amrex::Gpu::synchronize();
        }
    }
}
