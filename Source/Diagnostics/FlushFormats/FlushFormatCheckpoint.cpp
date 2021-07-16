#include "FlushFormatCheckpoint.H"

#include "BoundaryConditions/PML.H"
#include "Diagnostics/ParticleDiag/ParticleDiag.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParticleIO.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

namespace
{
    const std::string default_level_prefix {"Level_"};
}

void
FlushFormatCheckpoint::WriteToFile (
        const amrex::Vector<std::string> /*varnames*/,
        const amrex::Vector<amrex::MultiFab>& /*mf*/,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double /*time*/,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
        const std::string prefix, int file_min_digits,
        bool /*plot_raw_fields*/,
        bool /*plot_raw_fields_guards*/,
        bool /*plot_raw_rho*/, bool /*plot_raw_F*/,
        bool /*isBTD*/, int /*snapshotID*/,
        const amrex::Geometry& /*full_BTD_snapshot*/,
        bool /*isLastBTDFlush*/) const
{
    WARPX_PROFILE("FlushFormatCheckpoint::WriteToFile()");

    auto & warpx = WarpX::GetInstance();

    VisMF::Header::Version current_version = VisMF::GetHeaderVersion();
    VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    const std::string& checkpointname = amrex::Concatenate(prefix, iteration[0], file_min_digits);

    amrex::Print() << "  Writing checkpoint " << checkpointname << "\n";

    // const int nlevels = finestLevel()+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, default_level_prefix, nlev, true);

    WriteWarpXHeader(checkpointname, particle_diags, geom);

    WriteJobInfo(checkpointname);

    for (int lev = 0; lev < nlev; ++lev)
    {
        VisMF::Write(warpx.getEfield_fp(lev, 0),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Ex_fp"));
        VisMF::Write(warpx.getEfield_fp(lev, 1),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Ey_fp"));
        VisMF::Write(warpx.getEfield_fp(lev, 2),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Ez_fp"));
        VisMF::Write(warpx.getBfield_fp(lev, 0),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Bx_fp"));
        VisMF::Write(warpx.getBfield_fp(lev, 1),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "By_fp"));
        VisMF::Write(warpx.getBfield_fp(lev, 2),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Bz_fp"));
        if (warpx.getis_synchronized()) {
            // Need to save j if synchronized because after restart we need j to evolve E by dt/2.
            VisMF::Write(warpx.getcurrent_fp(lev, 0),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "jx_fp"));
            VisMF::Write(warpx.getcurrent_fp(lev, 1),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "jy_fp"));
            VisMF::Write(warpx.getcurrent_fp(lev, 2),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "jz_fp"));
        }

        if (lev > 0)
        {
            VisMF::Write(warpx.getEfield_cp(lev, 0),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Ex_cp"));
            VisMF::Write(warpx.getEfield_cp(lev, 1),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Ey_cp"));
            VisMF::Write(warpx.getEfield_cp(lev, 2),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Ez_cp"));
            VisMF::Write(warpx.getBfield_cp(lev, 0),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Bx_cp"));
            VisMF::Write(warpx.getBfield_cp(lev, 1),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "By_cp"));
            VisMF::Write(warpx.getBfield_cp(lev, 2),
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "Bz_cp"));
            if (warpx.getis_synchronized()) {
                // Need to save j if synchronized because after restart we need j to evolve E by dt/2.
                VisMF::Write(warpx.getcurrent_cp(lev, 0),
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "jx_cp"));
                VisMF::Write(warpx.getcurrent_cp(lev, 1),
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "jy_cp"));
                VisMF::Write(warpx.getcurrent_cp(lev, 2),
                             amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "jz_cp"));
            }
        }

        if (warpx.DoPML() && warpx.GetPML(lev)) {
            warpx.GetPML(lev)->CheckPoint(
                amrex::MultiFabFileFullPrefix(lev, checkpointname, default_level_prefix, "pml"));
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
