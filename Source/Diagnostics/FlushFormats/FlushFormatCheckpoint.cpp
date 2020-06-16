#include "FlushFormatCheckpoint.H"
#include "WarpX.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX_buildInfo.H>

using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

void
FlushFormatCheckpoint::WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev, const std::string prefix,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool plot_raw_rho, bool plot_raw_F) const
{
    WARPX_PROFILE("FlushFormatCheckpoint::WriteToFile()");

    auto & warpx = WarpX::GetInstance();

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    const std::string& checkpointname = amrex::Concatenate(prefix,iteration[0]);

    amrex::Print() << "  Writing checkpoint " << checkpointname << "\n";

    // const int nlevels = finestLevel()+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlev, true);

    WriteWarpXHeader(checkpointname, particle_diags);

    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlev; ++lev)
    {
        VisMF::Write(warpx.getEfield_fp(lev, 0),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ex_fp"));
        VisMF::Write(warpx.getEfield_fp(lev, 1),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ey_fp"));
        VisMF::Write(warpx.getEfield_fp(lev, 2),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ez_fp"));
        VisMF::Write(warpx.getBfield_fp(lev, 0),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bx_fp"));
        VisMF::Write(warpx.getBfield_fp(lev, 1),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "By_fp"));
        VisMF::Write(warpx.getBfield_fp(lev, 2),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bz_fp"));
        if (warpx.getis_synchronized()) {
            // Need to save j if synchronized because after restart we need j to evolve E by dt/2.
            VisMF::Write(warpx.getcurrent_fp(lev, 0),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jx_fp"));
            VisMF::Write(warpx.getcurrent_fp(lev, 1),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jy_fp"));
            VisMF::Write(warpx.getcurrent_fp(lev, 2),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jz_fp"));
        }

        if (lev > 0)
        {
            VisMF::Write(warpx.getEfield_cp(lev, 0),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ex_cp"));
            VisMF::Write(warpx.getEfield_cp(lev, 1),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ey_cp"));
            VisMF::Write(warpx.getEfield_cp(lev, 2),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Ez_cp"));
            VisMF::Write(warpx.getBfield_cp(lev, 0),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bx_cp"));
            VisMF::Write(warpx.getBfield_cp(lev, 1),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "By_cp"));
            VisMF::Write(warpx.getBfield_cp(lev, 2),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "Bz_cp"));
            if (warpx.getis_synchronized()) {
                // Need to save j if synchronized because after restart we need j to evolve E by dt/2.
                VisMF::Write(warpx.getcurrent_cp(lev, 0),
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jx_cp"));
                VisMF::Write(warpx.getcurrent_cp(lev, 1),
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jy_cp"));
                VisMF::Write(warpx.getcurrent_cp(lev, 2),
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "jz_cp"));
            }
        }

        if (warpx.DoPML() && warpx.GetPML(lev)) {
            warpx.GetPML(lev)->CheckPoint(
                amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "pml"));
        }
    }

    CheckpointParticles(checkpointname, particle_diags);

    VisMF::SetHeaderVersion(current_version);

}

void
FlushFormatCheckpoint::CheckpointParticles(
    const std::string& dir,
    const amrex::Vector<ParticleDiag>& particle_diags) const
{
    for (unsigned i = 0, n = particle_diags.size(); i < n; ++i) {
        particle_diags[i].getParticleContainer()->Checkpoint(
            dir, particle_diags[i].getSpeciesName());
    }
}
