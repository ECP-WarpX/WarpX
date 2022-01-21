#include "Particles/WarpXParticleContainer_fwd.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Gather/ScaleFields.H"
#include "Particles/ElementaryProcess/QEDPhotonEmission.H"

#include <AMReX_BaseFwd.H>

#include <AMReX_REAL.H>

void doPushPXQEDSync(
    const amrex::Real dt,
    const long particle_pusher_algo,
    WarpXParIter& pti,
    const long offset,
    const long np_to_push,
    const amrex::Real q,
    const amrex::Real m,
    int* AMREX_RESTRICT ion_lev,
    const GetParticlePosition& getPosition,
    amrex::FArrayBox const* exfab,
    amrex::FArrayBox const* eyfab,
    amrex::FArrayBox const* ezfab,
    amrex::FArrayBox const* bxfab,
    amrex::FArrayBox const* byfab,
    amrex::FArrayBox const* bzfab,
    ScaleFields scaleFields,
    amrex::GpuArray<amrex::Real, 3> dx_arr,
    amrex::GpuArray<amrex::Real, 3> xyzmin_arr,
    const amrex::Dim3 lo,
    const int n_rz_azimuthal_modes,
    const int nox,
    const bool galerkin_interpolation,
    const bool do_not_gather,
    const bool do_classical_radiation_reaction,
    QuantumSynchrotronEvolveOpticalDepth& evolve_opt,
    amrex::ParticleReal* AMREX_RESTRICT p_optical_depth_QSR,
    const amrex::Real chi_max);
