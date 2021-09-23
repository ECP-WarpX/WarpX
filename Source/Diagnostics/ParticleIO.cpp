/* Copyright 2019 Andrew Myers, Axel Huebl, David Grote
 * Luca Fedeli, Maxence Thevenet, Revathi Jambunathan
 * Weiqun Zhang, levinem, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Particles/ParticleIO.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBuffer.H"
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

#include <algorithm>
#include <array>
#include <istream>
#include <memory>
#include <string>
#include <sstream>
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
        WarpXParticleContainer* pc = allcontainers[i].get();
        std::string header_fn = dir + "/" + species_names[i] + "/Header";

        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(header_fn, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        std::string line, word;

        std::getline(is, line); // Version
        std::getline(is, line); // SpaceDim

        int nr;
        is >> nr;

        std::vector<std::string> real_comp_names;
        for (int j = 0; j < nr; ++j) {
            std::string comp_name;
            is >> comp_name;
            real_comp_names.push_back(comp_name);
        }

        for (auto const& comp : pc->getParticleRuntimeComps()) {
            auto search = std::find(real_comp_names.begin(), real_comp_names.end(), comp.first);
            if (search == real_comp_names.end()) {
                std::stringstream ss;
                ss << "Species " << species_names[i] << "needs runtime real component " << comp.first;
                ss << ", but it was not found in the checkpoint file. \n";
                amrex::Abort(ss.str());
            }
        }

        for (int j = PIdx::nattribs; j < nr; ++j) {
            const auto& comp_name = real_comp_names[j];
            auto current_comp_names = pc->getParticleComps();
            auto search = current_comp_names.find(comp_name);
            if (search == current_comp_names.end()) {
                amrex::Print() << "Runtime real component " << comp_name
                               << " was found in the checkpoint file, but it has not been added yet. "
                               << " Adding it now. \n";
                pc->AddRealComp(comp_name);
            }
        }

        int ni;
        is >> ni;

        std::vector<std::string> int_comp_names;
        for (int j = 0; j < ni; ++j) {
            std::string comp_name;
            is >> comp_name;
            int_comp_names.push_back(comp_name);
        }

        for (auto const& comp : pc->getParticleRuntimeiComps()) {
            auto search = std::find(int_comp_names.begin(), int_comp_names.end(), comp.first);
            if (search == int_comp_names.end()) {
                std::stringstream ss;
                ss << "Species " << species_names[i] << "needs runtime int component " << comp.first;
                ss << ", but it was not found in the checkpoint file. \n";
                amrex::Abort(ss.str());
            }
        }

        for (int j = 0; j < ni; ++j) {
            const auto& comp_name = int_comp_names[j];
            auto current_comp_names = pc->getParticleiComps();
            auto search = current_comp_names.find(comp_name);
            if (search == current_comp_names.end()) {
                amrex::Print() << "Runtime int component " << comp_name
                               << " was found in the checkpoint file, but it has not been added yet. "
                               << " Adding it now. \n";
                pc->AddIntComp(comp_name);
            }
        }

        pc->Restart(dir, species_names[i]);
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

void
PhysicalParticleContainer::ConvertUnits (ConvertDirection convert_direction)
{
    WARPX_PROFILE("PhysicalParticleContainer::ConvertUnits()");

    // Account for the special case of photons
    const auto t_mass =
        this->AmIA<PhysicalSpecies::photon>() ? PhysConst::m_e : this->getMass();

    particlesConvertUnits(convert_direction, this, t_mass);
}
