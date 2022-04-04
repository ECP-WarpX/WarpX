/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENT_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENT_H_

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

class LatticeElement
{
public:
    LatticeElement (std::string const & element_name);

    LatticeElement(LatticeElement const &) = delete;
    LatticeElement(LatticeElement &&) = delete;
    LatticeElement & operator=(LatticeElement const &) = delete;

    virtual ~LatticeElement () = default;

    bool isDefined () {return static_cast<int>(m_zstarts.size()) > 0;}

    AMREX_GPU_HOST_DEVICE
    const amrex::Vector<amrex::Real> & get_zstarts () const& {return m_zstarts;}
    AMREX_GPU_HOST_DEVICE
    const amrex::Vector<amrex::Real> & get_zends () const& {return m_zends;}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    virtual void get_field (const int /*ielement*/,
                            const amrex::ParticleReal /*x*/,
                            const amrex::ParticleReal /*y*/,
                            const amrex::ParticleReal /*z*/,
                            const amrex::ParticleReal /*ux*/,
                            const amrex::ParticleReal /*uy*/,
                            const amrex::ParticleReal /*uz*/,
                            const amrex::ParticleReal /*dt*/,
                            amrex::ParticleReal& /*Ex*/,
                            amrex::ParticleReal& /*Ey*/,
                            amrex::ParticleReal& /*Ez*/,
                            amrex::ParticleReal& /*Bx*/,
                            amrex::ParticleReal& /*By*/,
                            amrex::ParticleReal& /*Bz*/) {}

protected:

    amrex::Vector<amrex::Real> m_zstarts;
    amrex::Vector<amrex::Real> m_zends;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real hard_edged_fraction(const amrex::Real z,
                                    const amrex::Real ux,
                                    const amrex::Real uy,
                                    const amrex::Real uz,
                                    const amrex::Real dt,
                                    const int ielement)
    {

        using namespace amrex::literals;

        constexpr amrex::Real inv_c2 = 1._rt/(PhysConst::c*PhysConst::c);
        const amrex::ParticleReal gamma = std::sqrt(1._rt + (ux*ux + uy*uy + uz*uz)*inv_c2);
        const amrex::ParticleReal vz = uz/gamma;

        amrex::ParticleReal const z1 = z;
        amrex::ParticleReal const z2 = z + vz*dt;

        // This allows vz to be positive or negative
        amrex::ParticleReal const zl = std::min(z1, z2);
        amrex::ParticleReal const zr = std::max(z1, z2);

        amrex::Real const zs = m_zstarts[ielement];
        amrex::Real const ze = m_zends[ielement];

        // Calculate the residence correction
        // frac will be 1 if the step is completely inside the lens, between 0 and 1
        // when entering or leaving the lens, and otherwise 0.
        amrex::Real const fl = ( (zl >= zs && zl < ze) ? 1. : 0.);
        amrex::Real const fr = ( (zr >= zs && zr < ze) ? 1. : 0.);
        amrex::Real frac = fl;
        if (fl > fr) frac = (ze - zl)/(zr - zl);
        if (fr > fl) frac = (zr - zs)/(zr - zl);

        return frac;
    }

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENT_H_
