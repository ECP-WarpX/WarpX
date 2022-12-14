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

enum struct CollisionType { DeuteriumTritiumToNeutronHeliumFusion,
                            DeuteriumDeuteriumToProtonTritiumFusion,
                            DeuteriumDeuteriumToNeutronHeliumFusion,
                            DeuteriumHeliumToProtonHeliumFusion,
                            ProtonBoronToAlphasFusion,
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

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void get_collision_parameters (
        const amrex::ParticleReal& u1x, const amrex::ParticleReal& u1y,
        const amrex::ParticleReal& u1z, const amrex::ParticleReal& u2x,
        const amrex::ParticleReal& u2y, const amrex::ParticleReal& u2z,
        const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
        amrex::ParticleReal& E_coll, amrex::ParticleReal& v_coll,
        amrex::ParticleReal& lab_to_COM_factor )
    {
        // General notations in this function:
        //     x_sq denotes the square of x
        //     x_star denotes the value of x in the center of mass frame

        using namespace amrex::literals;

        constexpr auto one_pr = amrex::ParticleReal(1.);
        constexpr auto inv_four_pr = amrex::ParticleReal(1./4.);
        constexpr amrex::ParticleReal c_sq = PhysConst::c * PhysConst::c;
        constexpr amrex::ParticleReal inv_csq = one_pr / ( c_sq );

        const amrex::ParticleReal m1_sq = m1*m1;
        const amrex::ParticleReal m2_sq = m2*m2;

        // Compute Lorentz factor gamma in the lab frame
        const amrex::ParticleReal g1 = std::sqrt( one_pr + (u1x*u1x+u1y*u1y+u1z*u1z)*inv_csq );
        const amrex::ParticleReal g2 = std::sqrt( one_pr + (u2x*u2x+u2y*u2y+u2z*u2z)*inv_csq );

        // Compute momenta
        const amrex::ParticleReal p1x = u1x * m1;
        const amrex::ParticleReal p1y = u1y * m1;
        const amrex::ParticleReal p1z = u1z * m1;
        const amrex::ParticleReal p2x = u2x * m2;
        const amrex::ParticleReal p2y = u2y * m2;
        const amrex::ParticleReal p2z = u2z * m2;
        // Square norm of the total (sum between the two particles) momenta in the lab frame
        auto constexpr pow2 = [](double const x) { return x*x; };
        const amrex::ParticleReal p_total_sq = pow2(p1x + p2x) +
                                            pow2(p1y+p2y) +
                                            pow2(p1z+p2z);

        // Total energy in the lab frame
        const amrex::ParticleReal E_lab = (m1 * g1 + m2 * g2) * c_sq;
        // Total energy squared in the center of mass frame, calculated using the Lorentz invariance
        // of the four-momentum norm
        const amrex::ParticleReal E_star_sq = E_lab*E_lab - c_sq*p_total_sq;

        // Kinetic energy in the center of mass frame
        const amrex::ParticleReal E_star = std::sqrt(E_star_sq);
        E_coll = E_star - (m1 + m2)*c_sq;

        // Square of the norm of the momentum of one of the particles in the center of mass frame
        // Formula obtained by inverting E^2 = p^2*c^2 + m^2*c^4 in the COM frame for each particle
        // The expression below is specifically written in a form that avoids returning
        // small negative numbers due to machine precision errors, for low-energy particles
        const amrex::ParticleReal E_ratio = E_star/((m1 + m2)*c_sq);
        const amrex::ParticleReal p_star_sq = m1*m2*c_sq * ( pow2(E_ratio) - one_pr )
                + pow2(m1 - m2)*c_sq*inv_four_pr * pow2( E_ratio - 1._prt/E_ratio );

        // Lorentz factors in the center of mass frame
        const amrex::ParticleReal g1_star = std::sqrt(one_pr + p_star_sq / (m1_sq*c_sq));
        const amrex::ParticleReal g2_star = std::sqrt(one_pr + p_star_sq / (m2_sq*c_sq));

        // relative velocity in the center of mass frame
        v_coll = std::sqrt(p_star_sq) * (one_pr/(m1*g1_star) + one_pr/(m2*g2_star));

        // Cross sections and relative velocity are computed in the center of mass frame.
        // On the other hand, the particle densities (weight over volume) in the lab frame are used. To
        // take into account this discrepancy, we need to multiply the fusion probability by the ratio
        // between the Lorentz factors in the COM frame and the Lorentz factors in the lab frame
        // (see Perez et al., Phys.Plasmas.19.083104 (2012))
        lab_to_COM_factor = g1_star*g2_star/(g1*g2);
    }
}

#endif // BINARY_COLLISION_UTILS_H_
