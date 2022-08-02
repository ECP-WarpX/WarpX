/* Copyright 2019 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
#define WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_

#include "Utils/WarpXConst.H"

#include <AMReX_Math.H>
#include <AMReX_Random.H>

#include <cmath>  // isnan() isinf()
#include <limits> // numeric_limits<float>::min()

/* \brief Update particle velocities according to
 *        F. Perez et al., Phys.Plasmas.19.083104 (2012),
 *        which is based on Nanbu's method, PhysRevE.55.4642 (1997).
 *        @param[in] LmdD is max(Debye length, minimal interparticle distance).
 *        @param[in] L is the Coulomb log. A fixed L will be used if L > 0,
 *        otherwise L will be calculated based on the algorithm.
 *        To see if there are nan or inf updated velocities,
 *        compile with USE_ASSERTION=TRUE.
*/

template <typename T_PR, typename T_R>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumPerezElastic (
    T_PR& u1x, T_PR& u1y, T_PR& u1z, T_PR& u2x, T_PR& u2y, T_PR& u2z,
    T_PR const n1, T_PR const n2, T_PR const n12,
    T_PR const q1, T_PR const m1, T_PR const w1,
    T_PR const q2, T_PR const m2, T_PR const w2,
    T_R const dt, T_PR const L, T_PR const lmdD,
    amrex::RandomEngine const& engine)
{

    T_PR const diffx = amrex::Math::abs(u1x-u2x);
    T_PR const diffy = amrex::Math::abs(u1y-u2y);
    T_PR const diffz = amrex::Math::abs(u1z-u2z);
    T_PR const diffm = std::sqrt(diffx*diffx+diffy*diffy+diffz*diffz);
    T_PR const summm = std::sqrt(u1x*u1x+u1y*u1y+u1z*u1z) + std::sqrt(u2x*u2x+u2y*u2y+u2z*u2z);
    // If g = u1 - u2 = 0, do not collide.
    // Or if the relative difference is less than 1.0e-10.
    if ( diffm < std::numeric_limits<T_PR>::min() || diffm/summm < 1.0e-10 ) { return; }

    T_PR constexpr inv_c2 = T_PR(1.0) / ( PhysConst::c * PhysConst::c );

    // Compute Lorentz factor gamma
    T_PR const g1 = std::sqrt( T_PR(1.0) + (u1x*u1x+u1y*u1y+u1z*u1z)*inv_c2 );
    T_PR const g2 = std::sqrt( T_PR(1.0) + (u2x*u2x+u2y*u2y+u2z*u2z)*inv_c2 );

    // Compute momenta
    T_PR const p1x = u1x * m1;
    T_PR const p1y = u1y * m1;
    T_PR const p1z = u1z * m1;
    T_PR const p2x = u2x * m2;
    T_PR const p2y = u2y * m2;
    T_PR const p2z = u2z * m2;

    // Compute center-of-mass (COM) velocity and gamma
    T_PR const mass_g = m1 * g1 + m2 * g2;
    T_PR const vcx    = (p1x+p2x) / mass_g;
    T_PR const vcy    = (p1y+p2y) / mass_g;
    T_PR const vcz    = (p1z+p2z) / mass_g;
    T_PR const vcms   = vcx*vcx + vcy*vcy + vcz*vcz;
    T_PR const gc     = T_PR(1.0) / std::sqrt( T_PR(1.0) - vcms*inv_c2 );

    // Compute vc dot v1 and v2
    T_PR const vcDv1 = (vcx*u1x + vcy*u1y + vcz*u1z) / g1;
    T_PR const vcDv2 = (vcx*u2x + vcy*u2y + vcz*u2z) / g2;

    // Compute p1 star
    T_PR p1sx;
    T_PR p1sy;
    T_PR p1sz;
    if ( vcms > std::numeric_limits<T_PR>::min() )
    {
        T_PR const lorentz_tansform_factor =
            ( (gc-T_PR(1.0))/vcms*vcDv1 - gc )*m1*g1;
        p1sx = p1x + vcx*lorentz_tansform_factor;
        p1sy = p1y + vcy*lorentz_tansform_factor;
        p1sz = p1z + vcz*lorentz_tansform_factor;
    }
    else // If vcms = 0, don't do Lorentz-transform.
    {
        p1sx = p1x;
        p1sy = p1y;
        p1sz = p1z;
    }
    T_PR const p1sm = std::sqrt( p1sx*p1sx + p1sy*p1sy + p1sz*p1sz );

    // Compute gamma star
    T_PR const g1s = ( T_PR(1.0) - vcDv1*inv_c2 )*gc*g1;
    T_PR const g2s = ( T_PR(1.0) - vcDv2*inv_c2 )*gc*g2;

    // Compute the Coulomb log lnLmd
    T_PR lnLmd;
    if ( L > T_PR(0.0) ) { lnLmd = L; }
    else
    {
        // Compute b0
        T_PR const b0 = amrex::Math::abs(q1*q2) * inv_c2 /
               (T_PR(4.0)*MathConst::pi*PhysConst::ep0) * gc/mass_g *
               ( m1*g1s*m2*g2s/(p1sm*p1sm*inv_c2) + T_PR(1.0) );

        // Compute the minimal impact parameter
        constexpr T_PR hbar_pi = static_cast<T_PR>(PhysConst::hbar*MathConst::pi);
        const T_PR bmin = amrex::max(hbar_pi/p1sm, b0);

        // Compute the Coulomb log lnLmd
        lnLmd = amrex::max( T_PR(2.0),
                T_PR(0.5)*std::log(T_PR(1.0)+lmdD*lmdD/(bmin*bmin)) );
    }

    // Compute s
    const auto tts = m1*g1s*m2*g2s/(inv_c2*p1sm*p1sm) + T_PR(1.0);
    const auto tts2 = tts*tts;
    T_PR s = n1*n2/n12 * dt*lnLmd*q1*q1*q2*q2 /
          ( T_PR(4.0) * MathConst::pi * PhysConst::ep0 * PhysConst::ep0 *
            m1*g1*m2*g2/(inv_c2*inv_c2) ) * gc*p1sm/mass_g * tts2;

    // Compute s'
    const auto cbrt_n1 = std::cbrt(n1);
    const auto cbrt_n2 = std::cbrt(n2);
    const auto coeff = static_cast<T_PR>(
        std::pow(4.0*MathConst::pi/3.0,1.0/3.0));
    T_PR const vrel = mass_g*p1sm/(m1*g1s*m2*g2s*gc);
    T_PR const sp = coeff * n1*n2/n12 * dt * vrel * (m1+m2) /
          amrex::max( m1*cbrt_n1*cbrt_n1,
                      m2*cbrt_n2*cbrt_n2);

    // Determine s
    s = amrex::min(s,sp);

    // Get random numbers
    T_PR r = amrex::Random(engine);

    // Compute scattering angle
    T_PR cosXs;
    T_PR sinXs;
    if ( s <= T_PR(0.1) )
    {
        while ( true )
        {
            cosXs = T_PR(1.0) + s * std::log(r);
            // Avoid the bug when r is too small such that cosXs < -1
            if ( cosXs >= T_PR(-1.0) ) { break; }
            r = amrex::Random(engine);
        }
    }
    else if ( s > T_PR(0.1) && s <= T_PR(3.0) )
    {
        T_PR const Ainv = static_cast<T_PR>(
            0.0056958 + 0.9560202*s - 0.508139*s*s +
            0.47913906*s*s*s - 0.12788975*s*s*s*s + 0.02389567*s*s*s*s*s);
        cosXs = Ainv * std::log( std::exp(T_PR(-1.0)/Ainv) +
                T_PR(2.0) * r * std::sinh(T_PR(1.0)/Ainv) );
    }
    else if ( s > T_PR(3.0) && s <= T_PR(6.0) )
    {
        T_PR const A = T_PR(3.0) * std::exp(-s);
        cosXs = T_PR(1.0)/A * std::log( std::exp(-A) +
                T_PR(2.0) * r * std::sinh(A) );
    }
    else
    {
        cosXs = T_PR(2.0) * r - T_PR(1.0);
    }
    sinXs = std::sqrt(T_PR(1.0) - cosXs*cosXs);

    // Get random azimuthal angle
    T_PR const phis = amrex::Random(engine) * T_PR(2.0) * MathConst::pi;
    T_PR const cosphis = std::cos(phis);
    T_PR const sinphis = std::sin(phis);

    // Compute post-collision momenta pfs in COM
    T_PR p1fsx;
    T_PR p1fsy;
    T_PR p1fsz;
    // p1sp is the p1s perpendicular
    T_PR p1sp = std::sqrt( p1sx*p1sx + p1sy*p1sy );
    // Make sure p1sp is not almost zero
    if ( p1sp > std::numeric_limits<T_PR>::min() )
    {
        p1fsx = ( p1sx*p1sz/p1sp ) * sinXs*cosphis +
                ( p1sy*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sx           ) * cosXs;
        p1fsy = ( p1sy*p1sz/p1sp ) * sinXs*cosphis +
                (-p1sx*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sy           ) * cosXs;
        p1fsz = (-p1sp           ) * sinXs*cosphis +
                ( T_PR(0.0)      ) * sinXs*sinphis +
                ( p1sz           ) * cosXs;
        // Note a negative sign is different from
        // Eq. (12) in Perez's paper,
        // but they are the same due to the random nature of phis.
    }
    else
    {
        // If the previous p1sp is almost zero
        // x->y  y->z  z->x
        // This set is equivalent to the one in Nanbu's paper
        p1sp = std::sqrt( p1sy*p1sy + p1sz*p1sz );
        p1fsy = ( p1sy*p1sx/p1sp ) * sinXs*cosphis +
                ( p1sz*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sy           ) * cosXs;
        p1fsz = ( p1sz*p1sx/p1sp ) * sinXs*cosphis +
                (-p1sy*p1sm/p1sp ) * sinXs*sinphis +
                ( p1sz           ) * cosXs;
        p1fsx = (-p1sp           ) * sinXs*cosphis +
                ( T_PR(0.0)      ) * sinXs*sinphis +
                ( p1sx           ) * cosXs;
    }

    T_PR const p2fsx = -p1fsx;
    T_PR const p2fsy = -p1fsy;
    T_PR const p2fsz = -p1fsz;

    // Transform from COM to lab frame
    T_PR p1fx;    T_PR p2fx;
    T_PR p1fy;    T_PR p2fy;
    T_PR p1fz;    T_PR p2fz;
    if ( vcms > std::numeric_limits<T_PR>::min() )
    {
        T_PR const vcDp1fs = vcx*p1fsx + vcy*p1fsy + vcz*p1fsz;
        T_PR const vcDp2fs = vcx*p2fsx + vcy*p2fsy + vcz*p2fsz;
        T_PR const factor = (gc-T_PR(1.0))/vcms;
        T_PR const factor1 = factor*vcDp1fs + m1*g1s*gc;
        T_PR const factor2 = factor*vcDp2fs + m2*g2s*gc;
        p1fx = p1fsx + vcx * factor1;
        p1fy = p1fsy + vcy * factor1;
        p1fz = p1fsz + vcz * factor1;
        p2fx = p2fsx + vcx * factor2;
        p2fy = p2fsy + vcy * factor2;
        p2fz = p2fsz + vcz * factor2;
    }
    else // If vcms = 0, don't do Lorentz-transform.
    {
        p1fx = p1fsx;
        p1fy = p1fsy;
        p1fz = p1fsz;
        p2fx = p2fsx;
        p2fy = p2fsy;
        p2fz = p2fsz;
    }

    // Rejection method
    r = amrex::Random(engine);
    if ( w2 > r*amrex::max(w1, w2) )
    {
        u1x  = p1fx / m1;
        u1y  = p1fy / m1;
        u1z  = p1fz / m1;
#ifndef AMREX_USE_DPCPP
        AMREX_ASSERT(!std::isnan(u1x+u1y+u1z+u2x+u2y+u2z));
        AMREX_ASSERT(!std::isinf(u1x+u1y+u1z+u2x+u2y+u2z));
#endif
    }
    r = amrex::Random(engine);
    if ( w1 > r*amrex::max(w1, w2) )
    {
        u2x  = p2fx / m2;
        u2y  = p2fy / m2;
        u2z  = p2fz / m2;
#ifndef AMREX_USE_DPCPP
        AMREX_ASSERT(!std::isnan(u1x+u1y+u1z+u2x+u2y+u2z));
        AMREX_ASSERT(!std::isinf(u1x+u1y+u1z+u2x+u2y+u2z));
#endif
    }

}

#endif // WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
