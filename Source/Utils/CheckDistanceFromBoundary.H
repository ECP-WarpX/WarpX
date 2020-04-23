#ifndef WARPX_CHECKDISTANCEFROMBOUNDARY_H_
#define WARPX_CHECKDISTANCEFROMBOUNDARY_H_

#include "WarpXUtil.H"

/** \brief This function checks if 1 particle is too close to a domain boundary.
 *
 * \param[in] xp Particle position in x
 * \param[in] yp Particle position in y
 * \param[in] zp Particle position in z
 * \param[in] prob_lo Lower end of the domain in all directions
 * \param[in] prob_hi Higher end of the domain in all directions
 */
AMREX_FORCE_INLINE
void
CheckDistanceFromBoundary (
    const amrex::ParticleReal xp, const amrex::ParticleReal yp, const amrex::ParticleReal zp,
    const amrex::Real* const prob_lo, const amrex::Real* const prob_hi)
{
    amrex::ParticleReal machine_epsilon = std::numeric_limits<amrex::ParticleReal>::epsilon();
    WarpXUtilMsg::AlwaysAssert(
#if (AMREX_SPACEDIM == 3)
        abs(xp-prob_lo[0]) > machine_epsilon &&
        abs(xp-prob_hi[0]) > machine_epsilon &&
        abs(yp-prob_lo[1]) > machine_epsilon &&
        abs(yp-prob_hi[1]) > machine_epsilon &&
        abs(zp-prob_lo[2]) > machine_epsilon &&
        abs(zp-prob_hi[2]) > machine_epsilon,
#else // AMREX_SPACEDIM == 2
        abs(xp-prob_lo[0]) > machine_epsilon &&
        abs(xp-prob_hi[0]) > machine_epsilon &&
        abs(zp-prob_lo[1]) > machine_epsilon &&
        abs(zp-prob_hi[1]) > machine_epsilon,
#endif
        "ERROR: A particle was detected too close to a domain boundary (at a distance smaller than "
        "the machine epsilon), which may cause an error in amrex::Redistribute(). Please consider "
        "extending the domain size a tiny bit. This error easily occurs with the laser antenna, in "
        "which case extending the domain size by 1% in the dimension of laser propagation could fix "
        "the issue"
        );
}

#endif // WARPX_CHECKDISTANCEFROMBOUNDARY_H_
