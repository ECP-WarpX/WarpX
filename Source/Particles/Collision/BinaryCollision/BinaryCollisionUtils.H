/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef BINARY_COLLISION_UTILS_H_
#define BINARY_COLLISION_UTILS_H_

#include <string>

#include "Particles/MultiParticleContainer.H"

#include <AMReX_Math.H>

enum struct CollisionType { DeuteriumTritiumToNeutronHeliumFusion,
                            DeuteriumDeuteriumToProtonTritiumFusion,
                            DeuteriumDeuteriumToNeutronHeliumFusion,
                            DeuteriumHeliumToProtonHeliumFusion,
                            ProtonBoronToAlphasFusion,
                            DSMC,
                            PairwiseCoulomb,
                            Undefined };

enum struct NuclearFusionType {
                                DeuteriumTritiumToNeutronHelium,
                                DeuteriumDeuteriumToProtonTritium,
                                DeuteriumDeuteriumToNeutronHelium,
                                DeuteriumHeliumToProtonHelium,
                                ProtonBoronToAlphas,
                                Undefined };

namespace BinaryCollisionUtils{

    NuclearFusionType get_nuclear_fusion_type (std::string collision_name,
                                               MultiParticleContainer const * mypc);

    CollisionType get_collision_type (std::string collision_name,
                                      MultiParticleContainer const * mypc);

    CollisionType nuclear_fusion_type_to_collision_type (NuclearFusionType fusion_type);

    /**
     * \brief Return (relativistic) collision energy, collision speed and
     * Lorentz factor for transforming between the lab and center-of-momentum
     * frames.
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void get_collision_parameters (
        const amrex::ParticleReal& u1x, const amrex::ParticleReal& u1y,
        const amrex::ParticleReal& u1z, const amrex::ParticleReal& u2x,
        const amrex::ParticleReal& u2y, const amrex::ParticleReal& u2z,
        const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
        amrex::ParticleReal& E_kin_COM, amrex::ParticleReal& v_rel_COM,
        amrex::ParticleReal& lab_to_COM_lorentz_factor )
    {
        // General notations in this function:
        //     x_sq denotes the square of x
        //     x_star denotes the value of x in the center of mass frame

        using namespace amrex::literals;
        using namespace amrex::Math;

        constexpr auto one_pr = amrex::ParticleReal(1.);
        constexpr auto inv_four_pr = amrex::ParticleReal(1./4.);
        constexpr double c_sq = PhysConst::c * PhysConst::c;
        constexpr double inv_csq = 1.0 / c_sq;

        const amrex::ParticleReal m1_sq = m1*m1;
        const amrex::ParticleReal m2_sq = m2*m2;

        // Compute Lorentz factor gamma in the lab frame
        const double g1 = std::sqrt( 1.0 + static_cast<double>(u1x*u1x+u1y*u1y+u1z*u1z)*inv_csq );
        const double g2 = std::sqrt( 1.0 + static_cast<double>(u2x*u2x+u2y*u2y+u2z*u2z)*inv_csq );

        // Compute momenta
        const amrex::ParticleReal p1x = u1x * m1;
        const amrex::ParticleReal p1y = u1y * m1;
        const amrex::ParticleReal p1z = u1z * m1;
        const amrex::ParticleReal p2x = u2x * m2;
        const amrex::ParticleReal p2y = u2y * m2;
        const amrex::ParticleReal p2z = u2z * m2;
        // Square norm of the total (sum between the two particles) momenta in the lab frame
        const auto p_total_sq = static_cast<double>(
            powi<2>(p1x + p2x) + powi<2>(p1y + p2y) + powi<2>(p1z + p2z)
        );

        // Total energy in the lab frame
        // Note the use of `double` for energy since this calculation is
        // prone to error with single precision.
        const auto m1_dbl = static_cast<double>(m1);
        const auto m2_dbl = static_cast<double>(m2);
        const double E_lab = (m1_dbl * g1 + m2_dbl * g2) * c_sq;
        // Total energy squared in the center of mass frame, calculated using the Lorentz invariance
        // of the four-momentum norm
        const double E_star_sq = E_lab*E_lab - c_sq*p_total_sq;

        // Kinetic energy in the center of mass frame
        const double E_star = std::sqrt(E_star_sq);
        E_kin_COM = static_cast<amrex::ParticleReal>(E_star - (m1_dbl + m2_dbl)*c_sq);

        // Square of the norm of the momentum of one of the particles in the center of mass frame
        // Formula obtained by inverting E^2 = p^2*c^2 + m^2*c^4 in the COM frame for each particle
        // The expression below is specifically written in a form that avoids returning
        // small negative numbers due to machine precision errors, for low-energy particles
        const auto E_ratio = static_cast<amrex::ParticleReal>(E_star/((m1 + m2)*c_sq));
        const auto p_star_sq = static_cast<amrex::ParticleReal>(
            m1*m2*c_sq * ( powi<2>(E_ratio) - one_pr )
            + powi<2>(m1 - m2)*c_sq*inv_four_pr * powi<2>( E_ratio - 1._prt/E_ratio)
        );

        // Lorentz factors in the center of mass frame
        const auto g1_star = std::sqrt(one_pr + p_star_sq / static_cast<amrex::ParticleReal>(m1_sq*c_sq));
        const auto g2_star = std::sqrt(one_pr + p_star_sq / static_cast<amrex::ParticleReal>(m2_sq*c_sq));

        // relative velocity in the center of mass frame
        v_rel_COM = std::sqrt(p_star_sq) * (one_pr/(m1*g1_star) + one_pr/(m2*g2_star));

        // Cross sections and relative velocity are computed in the center of mass frame.
        // On the other hand, the particle densities (weight over volume) in the lab frame are used.
        // To take this disrepancy into account, it is needed to multiply the
        // collision probability by the ratio between the Lorentz factors in the
        // COM frame and the Lorentz factors in the lab frame (see
        // Perez et al., Phys.Plasmas.19.083104 (2012)). The correction factor
        // is calculated here.
        lab_to_COM_lorentz_factor = g1_star*g2_star/static_cast<amrex::ParticleReal>(g1*g2);
    }
}

#endif // BINARY_COLLISION_UTILS_H_
