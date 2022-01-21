#include "PushPX.H"

#include "Particles/ElementaryProcess/QEDInternals/QedChiFunctions.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Pusher/PushSelector.H"
#include "Particles/Gather/FieldGather.H"
#include "Particles/Gather/ScaleFields.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_GpuLaunch.H>
#include <AMReX_Particle.H>

using namespace amrex;

template <
    bool GalerkinInterpolation,
    bool DoNotGather,
    int NOX,
    bool DoCrr,
    long ParticlePusherAlgo
#ifdef WARPX_QED
    , bool DoSync
#endif
>
void doPushPX(
    const Real dt,
    const int np_to_push,
    ParticleReal* const AMREX_RESTRICT ux,
    ParticleReal* const AMREX_RESTRICT uy,
    ParticleReal* const AMREX_RESTRICT uz,
    const Real q,
    const Real m,
    int* AMREX_RESTRICT ion_lev,
    const GetParticlePosition& getPosition,
    const GetExternalEBField& getExternalEB,
    SetParticlePosition & setPosition,
    ScaleFields& scaleFields,
    Array4<const Real> const& ex_arr,
    Array4<const Real> const& ey_arr,
    Array4<const Real> const& ez_arr,
    Array4<const Real> const& bx_arr,
    Array4<const Real> const& by_arr,
    Array4<const Real> const& bz_arr,
    IndexType const ex_type,
    IndexType const ey_type,
    IndexType const ez_type,
    IndexType const bx_type,
    IndexType const by_type,
    IndexType const bz_type,
    GpuArray<Real, 3> dx_arr,
    GpuArray<Real, 3> xyzmin_arr,
    const Dim3 lo,
    const int n_rz_azimuthal_modes
#ifdef WARPX_QED
    , QuantumSynchrotronEvolveOpticalDepth& evolve_opt,
    amrex::ParticleReal* AMREX_RESTRICT p_optical_depth_QSR,
    const amrex::Real chi_max
#endif
    )
{
    ParallelFor( np_to_push, [=] AMREX_GPU_DEVICE (long ip)
    {
        ParticleReal xp, yp, zp;
        getPosition(ip, xp, yp, zp);

        ParticleReal Exp = 0._rt, Eyp = 0._rt, Ezp = 0._rt;
        ParticleReal Bxp = 0._rt, Byp = 0._rt, Bzp = 0._rt;

        if constexpr (!DoNotGather){
            // first gather E and B to the particle positions
            constexpr int tGal = (GalerkinInterpolation)? 1 : 0;
            doGatherShapeN<NOX,tGal>(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                           ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                           ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                           dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes);
        }
        // Externally applied E and B-field in Cartesian co-ordinates
        getExternalEB(ip, Exp, Eyp, Ezp, Bxp, Byp, Bzp);

        scaleFields(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp);

        doParticlePush<
            ParticlePusherAlgo,
            DoCrr, DoSync>(
                getPosition, setPosition, ip,
                ux[ip], uy[ip], uz[ip],
                Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                ion_lev ? ion_lev[ip] : 0,
                m, q,
#ifdef WARPX_QED
                chi_max,
#endif
                dt);

#ifdef WARPX_QED
        if constexpr (DoSync) {
            evolve_opt(ux[ip], uy[ip], uz[ip],
                       Exp, Eyp, Ezp,Bxp, Byp, Bzp,
                       dt, p_optical_depth_QSR[ip]);
        }
        else{
          ignore_unused(evolve_opt);
        }

#endif

    });
}

#ifdef WARPX_QED
template <
    bool GalerkinInterpolation,
    bool DoNotGather,
    int NOX,
    bool DoCrr,
    long ParticlePusherAlgo,
    typename ...Params
>
void doPushPX(bool do_sync, Params&&... params)
{
    if (do_sync)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            NOX,
            DoCrr,
            ParticlePusherAlgo,
            true>(std::forward<Params>(params)...);
    else
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            NOX,
            DoCrr,
            ParticlePusherAlgo,
            false>(std::forward<Params>(params)...);
}

#endif

template <
    bool GalerkinInterpolation,
    bool DoNotGather,
    int NOX,
    bool DoCrr,
    typename ...Params
>
void doPushPX(const long pusher_algo, Params&&... params)
{
    if (pusher_algo == 0)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            NOX,
            DoCrr,
            0>(std::forward<Params>(params)...);
    else if (pusher_algo == 1)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            NOX,
            DoCrr,
            1>(std::forward<Params>(params)...);
    else if (pusher_algo == 2)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            NOX,
            DoCrr,
            2>(std::forward<Params>(params)...);
    else
        Abort("ERR");
}

template <
    bool GalerkinInterpolation,
    bool DoNotGather,
    int NOX,
    typename ...Params
>
void doPushPX(const bool do_crr, Params&&... params)
{
    if (do_crr)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            NOX,
            true>(std::forward<Params>(params)...);
    else
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            NOX,
            false>(std::forward<Params>(params)...);
}

template <
    bool GalerkinInterpolation,
    bool DoNotGather,
    typename ...Params
>
void doPushPX(const int nox, Params&&... params)
{
    if (nox == 1)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            1>(std::forward<Params>(params)...);
    else if (nox == 2)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            2>(std::forward<Params>(params)...);
    else if (nox == 3)
        doPushPX<GalerkinInterpolation,
            DoNotGather,
            3>(std::forward<Params>(params)...);
    else
        Abort("ERR");
}



template <
    bool GalerkinInterpolation,
    typename ...Params
>
void doPushPX(const bool do_not_gather, Params&&... params)
{
    if (do_not_gather)
        doPushPX<GalerkinInterpolation,
            true>(std::forward<Params>(params)...);
    else
        doPushPX<GalerkinInterpolation,
            false>(std::forward<Params>(params)...);
}

void doPushPX(
    const Real dt,
    const long particle_pusher_algo,
    WarpXParIter& pti,
    const long offset,
    const long np_to_push,
    const Real q,
    const Real m,
    int* AMREX_RESTRICT ion_lev,
    const GetParticlePosition& getPosition,
    FArrayBox const* exfab,
    FArrayBox const* eyfab,
    FArrayBox const* ezfab,
    FArrayBox const* bxfab,
    FArrayBox const* byfab,
    FArrayBox const* bzfab,
    ScaleFields scaleFields,
    GpuArray<amrex::Real, 3> dx_arr,
    GpuArray<amrex::Real, 3> xyzmin_arr,
    const Dim3 lo,
    const int n_rz_azimuthal_modes,
    const int nox,
    const bool galerkin_interpolation,
    const bool do_not_gather,
    const bool do_classical_radiation_reaction
#ifdef WARPX_QED
        , bool do_sync,
        QuantumSynchrotronEvolveOpticalDepth& evolve_opt,
        amrex::ParticleReal* AMREX_RESTRICT p_optical_depth_QSR,
        const amrex::Real chi_max
#endif
    )
{

    Array4<const Real> const& ex_arr = exfab->array();
    Array4<const Real> const& ey_arr = eyfab->array();
    Array4<const Real> const& ez_arr = ezfab->array();
    Array4<const Real> const& bx_arr = bxfab->array();
    Array4<const Real> const& by_arr = byfab->array();
    Array4<const Real> const& bz_arr = bzfab->array();

    IndexType const ex_type = exfab->box().ixType();
    IndexType const ey_type = eyfab->box().ixType();
    IndexType const ez_type = ezfab->box().ixType();
    IndexType const bx_type = bxfab->box().ixType();
    IndexType const by_type = byfab->box().ixType();
    IndexType const bz_type = bzfab->box().ixType();

    auto setPosition = SetParticlePosition(pti, offset);

    const auto getExternalEB = GetExternalEBField(pti, offset);

    auto& attribs = pti.GetAttribs();
    ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr() + offset;
    ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr() + offset;
    ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr() + offset;

    if(galerkin_interpolation)
        doPushPX<true>(
            do_not_gather,
            nox,
            do_classical_radiation_reaction,
            particle_pusher_algo,
#ifdef WARPX_QED
            do_sync,
#endif
            dt,
            np_to_push,
            ux, uy, uz,
            q, m,
            ion_lev,
            getPosition,
            getExternalEB,
            setPosition,
            scaleFields,
            ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
            ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
            dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes
#ifdef WARPX_QED
            , evolve_opt,
            p_optical_depth_QSR,
            chi_max
#endif
        );
    else
        doPushPX<false>(
            do_not_gather,
            nox,
            do_classical_radiation_reaction,
            particle_pusher_algo,
#ifdef WARPX_QED
            do_sync,
#endif
            dt,
            np_to_push,
            ux, uy, uz,
            q, m,
            ion_lev,
            getPosition,
            getExternalEB,
            setPosition,
            scaleFields,
            ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
            ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
            dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes
#ifdef WARPX_QED
            , evolve_opt,
            p_optical_depth_QSR,
            chi_max
#endif
        );
}
