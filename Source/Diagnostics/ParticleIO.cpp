/* Copyright 2019 Andrew Myers, Axel Huebl, David Grote
 * Luca Fedeli, Maxence Thevenet, Revathi Jambunathan
 * Weiqun Zhang, levinem, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Particles/MultiParticleContainer.H"
#include "Particles/PhysicalParticleContainer.H"
#include "Particles/RigidInjectedParticleContainer.H"
#include "Particles/SpeciesPhysicalProperties.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX_BLassert.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParticleIO.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <istream>
#include <memory>
#include <string>
#include <vector>

using namespace amrex;

void
RigidInjectedParticleContainer::ReadHeader (std::istream& is)
{
    is >> charge >> mass;
    WarpX::GotoNextLine(is);

    int nlevs;
    is >> nlevs;
    WarpX::GotoNextLine(is);

    AMREX_ASSERT(zinject_plane_levels.size() == 0);
    AMREX_ASSERT(done_injecting.size() == 0);

    for (int i = 0; i < nlevs; ++i)
    {
        int zinject_plane_tmp;
        is >> zinject_plane_tmp;
        zinject_plane_levels.push_back(zinject_plane_tmp);
        WarpX::GotoNextLine(is);
    }

    for (int i = 0; i < nlevs; ++i)
    {
        int done_injecting_tmp;
        is >> done_injecting_tmp;
        done_injecting.push_back(done_injecting_tmp);
        WarpX::GotoNextLine(is);
    }
}

void
RigidInjectedParticleContainer::WriteHeader (std::ostream& os) const
{
    // no need to write species_id
    os << charge << " " << mass << "\n";
    int nlevs = zinject_plane_levels.size();
    os << nlevs << "\n";
    for (int i = 0; i < nlevs; ++i)
    {
        os << zinject_plane_levels[i] << "\n";
    }
    for (int i = 0; i < nlevs; ++i)
    {
        os << done_injecting[i] << "\n";
    }
}

void
WarpXParticleContainer::ReadHeader (std::istream& is)
{
    is >> charge >> mass;
    WarpX::GotoNextLine(is);
}

void
WarpXParticleContainer::WriteHeader (std::ostream& os) const
{
    // no need to write species_id
    os << charge << " " << mass << "\n";
}

void
MultiParticleContainer::Restart (const std::string& dir)
{
    for (unsigned i = 0, n = species_names.size(); i < n; ++i) {
        allcontainers[i]->Restart(dir, species_names[i]);
    }
}

void
MultiParticleContainer::ReadHeader (std::istream& is)
{
    for (auto& pc : allcontainers) {
        pc->ReadHeader(is);
    }
}

void
MultiParticleContainer::WriteHeader (std::ostream& os) const
{
    for (const auto& pc : allcontainers) {
        pc->WriteHeader(os);
    }
}

// Particle momentum is defined as gamma*velocity, which is neither
// SI mass*gamma*velocity nor normalized gamma*velocity/c.
// This converts momentum to SI units (or vice-versa) to write SI data
// to file.
// Photons are a special case, since particle momentum is defined as
// (photon_energy/(m_e * c) ) * u, where u is the photon direction (a
// unit vector).
void
PhysicalParticleContainer::ConvertUnits(ConvertDirection convert_direction)
{
    WARPX_PROFILE("PhysicalParticleContainer::ConvertUnits()");

    // Compute conversion factor
    auto factor = 1_rt;

    // Account for the special case of photons
    const auto t_mass =
        AmIA<PhysicalSpecies::photon>() ? PhysConst::m_e : mass;

    if (convert_direction == ConvertDirection::WarpX_to_SI){
        factor = t_mass;
    } else if (convert_direction == ConvertDirection::SI_to_WarpX){
        factor = 1._rt/t_mass;
    }

    const int nLevels = finestLevel();
    for (int lev=0; lev<=nLevels; lev++){
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // - momenta are stored as a struct of array, in `attribs`
            auto& attribs = pti.GetAttribs();
            ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
            // Loop over the particles and convert momentum
            const long np = pti.numParticles();
            ParallelFor( np,
                [=] AMREX_GPU_DEVICE (long i) {
                    ux[i] *= factor;
                    uy[i] *= factor;
                    uz[i] *= factor;
                }
            );
        }
    }
}
