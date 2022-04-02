/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENT_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENT_H_

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

    const amrex::Vector<amrex::Real> & get_zstarts () const {return m_zstarts;}
    const amrex::Vector<amrex::Real> & get_zends () const {return m_zends;}

    // AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    virtual void get_field (const int /*ielement*/,
                            const amrex::ParticleReal /*x*/,
                            const amrex::ParticleReal /*y*/,
                            const amrex::ParticleReal /*z*/,
                            amrex::ParticleReal& /*Ex*/,
                            amrex::ParticleReal& /*Ey*/,
                            amrex::ParticleReal& /*Ez*/,
                            amrex::ParticleReal& /*Bx*/,
                            amrex::ParticleReal& /*By*/,
                            amrex::ParticleReal& /*Bz*/) {}

protected:

    amrex::Vector<amrex::Real> m_zstarts;
    amrex::Vector<amrex::Real> m_zends;

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_LATTICEELEMENT_H_
