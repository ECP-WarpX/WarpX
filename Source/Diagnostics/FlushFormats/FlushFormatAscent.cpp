#include "FlushFormatAscent.H"
#include "WarpX.H"

#include <AMReX.H>

using namespace amrex;

void
FlushFormatAscent::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<amrex::MultiFab>& mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
    const std::string prefix, bool plot_raw_fields,
    bool plot_raw_fields_guards, bool /*plot_raw_rho*/, bool plot_raw_F) const
{
#ifdef AMREX_USE_ASCENT
    auto & warpx = WarpX::GetInstance();

    amrex::Real convert_factor = 1._rt/(warpx.gamma_boost * (1._rt - warpx.beta_boost) );

    // hack mesh labels to make box square
    auto geom_hack = geom;
    for (int i = 0; i < geom_hack.size(); ++i) {
        Real new_lo[AMREX_SPACEDIM];
        Real new_hi[AMREX_SPACEDIM];
        const Real* current_lo = geom_hack[i].ProbLo();
        const Real* current_hi = geom_hack[i].ProbHi();

        new_lo[0] = current_lo[0];// / convert_factor;
        new_lo[1] = current_lo[1];// / convert_factor;
        new_lo[2] = current_lo[2] / convert_factor;
        new_hi[0] = current_hi[0];// / convert_factor;
        new_hi[1] = current_hi[1];// / convert_factor;
        new_hi[2] = current_hi[2] / convert_factor;

        geom_hack[i].ProbDomain(RealBox(new_lo, new_hi));

        std::cout << "buffer_i = " << i << std::endl;
        std::cout << "    domain.lo = "
            << geom[i].ProbLo()[0] << " " << geom_hack[i].ProbLo()[1] << " " << geom_hack[i].ProbLo()[2]
            << std::endl;
        std::cout << "    domain.hi = "
            << geom[i].ProbHi()[0] << " " << geom_hack[i].ProbHi()[1] << " " << geom_hack[i].ProbHi()[2]
            << std::endl;
    }

    // wrap mesh data
    conduit::Node bp_mesh;
    amrex::MultiLevelToBlueprint(
        nlev, amrex::GetVecOfConstPtrs(mf), varnames, geom_hack, time, iteration, warpx.refRatio(), bp_mesh);

    // hack particle positions to make box square
    for (unsigned k = 0, n = particle_diags.size(); k < n; ++k) {
        WarpXParticleContainer* pc = particle_diags[k].getParticleContainer();
        for (auto lev = 0; lev <= pc->finestLevel(); lev++)
        {
            for (WarpXParIter pti(*pc, lev); pti.isValid(); ++pti)
            {
                auto const numParticleOnTile = pti.numParticles();
                auto& aos = pti.GetArrayOfStructs();
                for (auto i = 0; i < numParticleOnTile; i++)
                {
                    //aos[i].pos(0) /= convert_factor;
                    //aos[i].pos(1) /= convert_factor;
                    aos[i].pos(2) /= convert_factor;
                }
            }
        }
    }

    WriteParticles(particle_diags, bp_mesh);

    // If you want to save blueprint HDF5 files w/o using an Ascent
    // extract, you can call the following AMReX helper:
    // const auto step = istep[0];
    // WriteBlueprintFiles(bp_mesh,"bp_export",step,"hdf5");

    ascent::Ascent ascent;
    conduit::Node opts;
    opts["exceptions"] = "catch";
    opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
    ascent.open(opts);
    ascent.publish(bp_mesh);
    conduit::Node actions;
    ascent.execute(actions);
    ascent.close();

    // revert particle position hack
    for (unsigned k = 0, n = particle_diags.size(); k < n; ++k) {
        WarpXParticleContainer* pc = particle_diags[k].getParticleContainer();
        for (auto lev = 0; lev <= pc->finestLevel(); lev++)
        {
            for (WarpXParIter pti(*pc, lev); pti.isValid(); ++pti)
            {
                auto const numParticleOnTile = pti.numParticles();
                auto& aos = pti.GetArrayOfStructs();
                for (auto i = 0; i < numParticleOnTile; i++)
                {
                    //aos[i].pos(0) *= convert_factor;
                    //aos[i].pos(1) *= convert_factor;
                    aos[i].pos(2) *= convert_factor;
                }
            }
        }
    }

#else
    amrex::ignore_unused(varnames, mf, geom, iteration, time,
        particle_diags, nlev);
#endif // AMREX_USE_ASCENT
    amrex::ignore_unused(prefix, plot_raw_fields, plot_raw_fields_guards,
        plot_raw_F);
}

#ifdef AMREX_USE_ASCENT
void
FlushFormatAscent::WriteParticles(const amrex::Vector<ParticleDiag>& particle_diags, conduit::Node& a_bp_mesh) const
{
    // wrap particle data for each species
    // we prefix the fields with "particle_{species_name}" b/c we
    // want to to uniquely name all the fields that can be plotted

    for (unsigned i = 0, n = particle_diags.size(); i < n; ++i) {
        Vector<std::string> particle_varnames;
        Vector<std::string> particle_int_varnames;
        std::string prefix = "particle_" + particle_diags[i].getSpeciesName();

        // Get pc for species
        // auto& pc = mypc->GetParticleContainer(i);
        WarpXParticleContainer* pc = particle_diags[i].getParticleContainer();

        // get names of real comps
        std::map<std::string, int> real_comps_map = pc->getParticleComps();

        // WarpXParticleContainer compile-time extra AoS attributes (Real): 0
        // WarpXParticleContainer compile-time extra AoS attributes (int): 0

        // WarpXParticleContainer compile-time extra SoA attributes (Real): PIdx::nattribs
        // not an efficient search, but N is small...
        for(int i = 0; i < PIdx::nattribs; ++i)
        {
            auto rvn_it = real_comps_map.begin();
            for (; rvn_it != real_comps_map.end(); ++rvn_it)
                if (rvn_it->second == i)
                    break;
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                rvn_it != real_comps_map.end(),
                "Ascent: SoA real attribute not found");
            std::string varname = rvn_it->first;
            particle_varnames.push_back(prefix + "_" + varname);
        }
        // WarpXParticleContainer compile-time extra SoA attributes (int): 0

        // WarpXParticleContainer "runtime" SoA attributes (Real), e.g QED: to do

        // wrap pc for current species into a blueprint topology
        amrex::ParticleContainerToBlueprint(*pc,
                                            particle_varnames,
                                            particle_int_varnames,
                                            a_bp_mesh,
                                            prefix);
    }
}
#endif // AMREX_USE_ASCENT
