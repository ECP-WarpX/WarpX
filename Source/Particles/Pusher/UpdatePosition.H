#ifndef WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_
#define WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_

#include <AMReX_FArrayBox.H>
#include <WarpXConst.H>
#include <AMReX_REAL.H>

/** \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdatePosition(
    amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
    const amrex::ParticleReal ux, const amrex::ParticleReal uy, const amrex::ParticleReal uz,
    const amrex::Real dt )
{

    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);

    // Compute inverse Lorentz factor
    const amrex::Real inv_gamma = 1./std::sqrt(1. + (ux*ux + uy*uy + uz*uz)*inv_c2);
    // Update positions over one time step
    x += ux * inv_gamma * dt;
#if (AMREX_SPACEDIM == 3) || (defined WARPX_DIM_RZ) // RZ pushes particles in 3D
    y += uy * inv_gamma * dt;
#endif
    z += uz * inv_gamma * dt;

}

/* \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdatePositionGalilean(
    amrex::Real& x, amrex::Real& y, amrex::Real& z,
    const amrex::Real ux, const amrex::Real uy, const amrex::Real uz,
    const amrex::Real vx, const amrex::Real vy, const amrex::Real vz,
    const amrex::Real dt )
{

    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);

    // Compute inverse Lorentz factor
    const amrex::Real inv_gamma = 1./std::sqrt(1. + (ux*ux + uy*uy + uz*uz)*inv_c2);
    // Update positions over one time step
    x += ux * inv_gamma * dt - vx * dt;  //oshapoval
#if (AMREX_SPACEDIM == 3) || (defined WARPX_RZ) // RZ pushes particles in 3D
    y += uy * inv_gamma * dt - vy * dt;
#endif
    z += uz * inv_gamma * dt - vz * dt;

}


#endif // WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_
