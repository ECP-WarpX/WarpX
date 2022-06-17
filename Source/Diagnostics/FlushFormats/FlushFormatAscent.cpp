#include "FlushFormatAscent.H"

#include "WarpX.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX.H>
#include <AMReX_REAL.H>

using namespace amrex;

void
FlushFormatAscent::WriteToFile (
    const amrex::Vector<std::string> varnames,
    const amrex::Vector<amrex::MultiFab>& mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
    const std::string prefix, int /*file_min_digits*/, bool plot_raw_fields,
    bool plot_raw_fields_guards,
    bool /*isBTD*/, int /*snapshotID*/, const amrex::Geometry& /*full_BTD_snapshot*/,
    bool /*isLastBTDFlush*/, const amrex::Vector<int>& /* totalParticlesFlushedAlready*/) const
{
#ifdef AMREX_USE_ASCENT
    WARPX_PROFILE("FlushFormatAscent::WriteToFile()");

    auto & warpx = WarpX::GetInstance();

    // wrap mesh data
    WARPX_PROFILE_VAR("FlushFormatAscent::WriteToFile::MultiLevelToBlueprint", prof_ascent_mesh_blueprint);
    conduit::Node bp_mesh;
    amrex::MultiLevelToBlueprint(
        nlev, amrex::GetVecOfConstPtrs(mf), varnames, geom, time, iteration, warpx.refRatio(), bp_mesh);
    WARPX_PROFILE_VAR_STOP(prof_ascent_mesh_blueprint);

    WARPX_PROFILE_VAR("FlushFormatAscent::WriteToFile::WriteParticles", prof_ascent_particles);
    WriteParticles(particle_diags, bp_mesh);
    WARPX_PROFILE_VAR_STOP(prof_ascent_particles);

    // If you want to save blueprint HDF5 files w/o using an Ascent
    // extract, you can call the following AMReX helper:
    // const auto step = istep[0];
    // WriteBlueprintFiles(bp_mesh,"bp_export",step,"hdf5");

    WARPX_PROFILE_VAR("FlushFormatAscent::WriteToFile::publish", prof_ascent_publish);
    ascent::Ascent ascent;
    conduit::Node opts;
    opts["exceptions"] = "catch";
    opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
    ascent.open(opts);
    ascent.publish(bp_mesh);
    WARPX_PROFILE_VAR_STOP(prof_ascent_publish);

    WARPX_PROFILE_VAR("FlushFormatAscent::WriteToFile::execute", prof_ascent_execute);
    conduit::Node actions;
    ascent.execute(actions);
    ascent.close();
    WARPX_PROFILE_VAR_STOP(prof_ascent_execute);

#else
    amrex::ignore_unused(varnames, mf, geom, iteration, time,
        particle_diags, nlev);
#endif // AMREX_USE_ASCENT
    amrex::ignore_unused(prefix, plot_raw_fields, plot_raw_fields_guards);
}

#ifdef AMREX_USE_ASCENT
void
FlushFormatAscent::WriteParticles(const amrex::Vector<ParticleDiag>& particle_diags, conduit::Node& a_bp_mesh) const
{
    WARPX_PROFILE("FlushFormatAscent::WriteParticles()");

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
        for(int j = 0; j < PIdx::nattribs; ++j)
        {
            auto rvn_it = real_comps_map.begin();
            for (; rvn_it != real_comps_map.end(); ++rvn_it)
                if (rvn_it->second == j)
                    break;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
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
