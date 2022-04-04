/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDQUADRUPOLE_H_
#define WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDQUADRUPOLE_H_

#include "LatticeElement.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

class HardEdgedQuadrupole
    : public LatticeElement
{
public:

    HardEdgedQuadrupole ();

    ~HardEdgedQuadrupole () = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void get_field (const int ielement,
                    const amrex::ParticleReal x,
                    const amrex::ParticleReal y,
                    const amrex::ParticleReal z,
                    const amrex::ParticleReal ux,
                    const amrex::ParticleReal uy,
                    const amrex::ParticleReal uz,
                    const amrex::ParticleReal dt,
                    amrex::ParticleReal& Ex,
                    amrex::ParticleReal& Ey,
                    amrex::ParticleReal& /*Ez*/,
                    amrex::ParticleReal& /*Bx*/,
                    amrex::ParticleReal& /*By*/,
                    amrex::ParticleReal& /*Bz*/)
    {

        amrex::Real const frac = hard_edged_fraction(z, ux, uy, uz, dt, ielement);

        amrex::Real dEdx = frac*m_dEdx[ielement];

        Ex = +x*dEdx;
        Ey = -y*dEdx;

    }

private:

    std::string const element_name = "quad";

    amrex::Vector<amrex::Real> m_dEdx;

};

#endif // WARPX_ACCELERATORLATTICE_LATTICEELEMENTS_HARDEDGEDQUADRUPOLE_H_
