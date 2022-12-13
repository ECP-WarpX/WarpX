/* Copyright 2021 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "BinaryCollisionUtils.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <string>

namespace BinaryCollisionUtils{

    NuclearFusionType get_nuclear_fusion_type (const std::string collision_name,
                                               MultiParticleContainer const * const mypc)
        {
            const amrex::ParmParse pp_collision_name(collision_name);
            amrex::Vector<std::string> species_names;
            pp_collision_name.getarr("species", species_names);
            auto& species1 = mypc->GetParticleContainerFromName(species_names[0]);
            auto& species2 = mypc->GetParticleContainerFromName(species_names[1]);
            amrex::Vector<std::string> product_species_name;
            pp_collision_name.getarr("product_species", product_species_name);

            if ((species1.AmIA<PhysicalSpecies::hydrogen2>() && species2.AmIA<PhysicalSpecies::hydrogen3>())
                ||
                (species1.AmIA<PhysicalSpecies::hydrogen3>() && species2.AmIA<PhysicalSpecies::hydrogen2>())
                )
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    product_species_name.size() == 2u,
                    "ERROR: Deuterium-tritium fusion must contain exactly two product species");
                auto& product_species1 = mypc->GetParticleContainerFromName(product_species_name[0]);
                auto& product_species2 = mypc->GetParticleContainerFromName(product_species_name[1]);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    (product_species1.AmIA<PhysicalSpecies::helium4>() && product_species2.AmIA<PhysicalSpecies::neutron>())
                    ||
                    (product_species1.AmIA<PhysicalSpecies::neutron>() && product_species2.AmIA<PhysicalSpecies::helium4>()),
                    "ERROR: Product species of deuterium-tritium fusion must be of type neutron and helium4");
                return NuclearFusionType::DeuteriumTritiumToNeutronHelium;
            }
            else if (species1.AmIA<PhysicalSpecies::hydrogen2>() && species2.AmIA<PhysicalSpecies::hydrogen2>())
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    product_species_name.size() == 2u,
                    "ERROR: Deuterium-deuterium fusion must contain exactly two product species");
                auto& product_species1 = mypc->GetParticleContainerFromName(product_species_name[0]);
                auto& product_species2 = mypc->GetParticleContainerFromName(product_species_name[1]);
                if (
                    (product_species1.AmIA<PhysicalSpecies::helium3>() && product_species2.AmIA<PhysicalSpecies::neutron>())
                  ||(product_species1.AmIA<PhysicalSpecies::neutron>() && product_species2.AmIA<PhysicalSpecies::helium3>())){
                    return NuclearFusionType::DeuteriumDeuteriumToNeutronHelium;
                } else if (
                    (product_species1.AmIA<PhysicalSpecies::hydrogen3>() && product_species2.AmIA<PhysicalSpecies::hydrogen1>())
                  ||(product_species1.AmIA<PhysicalSpecies::hydrogen1>() && product_species2.AmIA<PhysicalSpecies::hydrogen3>())){
                    return NuclearFusionType::DeuteriumDeuteriumToProtonTritium;
                } else {
                    WARPX_ABORT_WITH_MESSAGE("ERROR: Product species of deuterium-deuterium fusion must be of type helium3 and neutron, or hydrogen3 and hydrogen1");
                }
            }
            else if ((species1.AmIA<PhysicalSpecies::hydrogen2>() && species2.AmIA<PhysicalSpecies::helium3>())
                ||
                (species1.AmIA<PhysicalSpecies::helium3>() && species2.AmIA<PhysicalSpecies::hydrogen2>())
                )
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    product_species_name.size() == 2u,
                    "ERROR: Deuterium-helium fusion must contain exactly two product species");
                auto& product_species1 = mypc->GetParticleContainerFromName(product_species_name[0]);
                auto& product_species2 = mypc->GetParticleContainerFromName(product_species_name[1]);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    (product_species1.AmIA<PhysicalSpecies::helium4>() && product_species2.AmIA<PhysicalSpecies::hydrogen1>())
                    ||
                    (product_species1.AmIA<PhysicalSpecies::hydrogen1>() && product_species2.AmIA<PhysicalSpecies::helium4>()),
                    "ERROR: Product species of deuterium-helium fusion must be of type hydrogen1 and helium4");
                return NuclearFusionType::DeuteriumHeliumToProtonHelium;
            }
            else if ((species1.AmIA<PhysicalSpecies::hydrogen1>() && species2.AmIA<PhysicalSpecies::boron11>())
                ||
                (species1.AmIA<PhysicalSpecies::boron11>() && species2.AmIA<PhysicalSpecies::hydrogen1>())
                )
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    product_species_name.size() == 1,
                    "ERROR: Proton-boron must contain exactly one product species");
                auto& product_species = mypc->GetParticleContainerFromName(product_species_name[0]);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    product_species.AmIA<PhysicalSpecies::helium4>(),
                    "ERROR: Product species of proton-boron fusion must be of type helium4");
                return NuclearFusionType::ProtonBoronToAlphas;
            }
            WARPX_ABORT_WITH_MESSAGE("Binary nuclear fusion not implemented between species " +
                         species_names[0] + " of type " + species1.getSpeciesTypeName() +
                         " and species " + species_names[1] + " of type " +
                         species2.getSpeciesTypeName());
            return NuclearFusionType::Undefined;
        }

    CollisionType get_collision_type (const std::string collision_name,
                                      MultiParticleContainer const * const mypc)
        {
            const amrex::ParmParse pp_collision_name(collision_name);
            std::string type;
            pp_collision_name.get("type", type);
            if (type == "nuclearfusion") {
                const NuclearFusionType fusion_type = get_nuclear_fusion_type(collision_name, mypc);
                return nuclear_fusion_type_to_collision_type(fusion_type);
            }
            WARPX_ABORT_WITH_MESSAGE(type + " is not a valid type of collision that creates new particles");
            return CollisionType::Undefined;
        }

    CollisionType nuclear_fusion_type_to_collision_type (const NuclearFusionType fusion_type)
        {
            if (fusion_type == NuclearFusionType::DeuteriumTritiumToNeutronHelium) {
                return CollisionType::DeuteriumTritiumToNeutronHeliumFusion;
            }
            if (fusion_type == NuclearFusionType::DeuteriumDeuteriumToProtonTritium) {
                return CollisionType::DeuteriumDeuteriumToProtonTritiumFusion;
            }
            if (fusion_type == NuclearFusionType::DeuteriumDeuteriumToNeutronHelium) {
                return CollisionType::DeuteriumDeuteriumToNeutronHeliumFusion;
            }
            if (fusion_type == NuclearFusionType::DeuteriumHeliumToProtonHelium) {
                return CollisionType::DeuteriumHeliumToProtonHeliumFusion;
            }
            if (fusion_type == NuclearFusionType::ProtonBoronToAlphas) {
                return CollisionType::ProtonBoronToAlphasFusion;
            }
            WARPX_ABORT_WITH_MESSAGE("Invalid nuclear fusion type");
            return CollisionType::Undefined;
        }

    void get_collision_parameters(
        const amrex::ParticleReal& u1x, const amrex::ParticleReal& u1y,
        const amrex::ParticleReal& u1z, const amrex::ParticleReal& u2x,
        const amrex::ParticleReal& u2y, const amrex::ParticleReal& u2z,
        const amrex::ParticleReal& m1, const amrex::ParticleReal& m2,
        amrex::ParticleReal& E_coll, amrex::ParticleReal& v_coll,
        amrex::ParticleReal& lab_to_COM_factor)
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
