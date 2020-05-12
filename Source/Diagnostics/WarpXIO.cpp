/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Axel Huebl
 * Burlen Loring, David Grote, Gunther H. Weber
 * Junmin Gu, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "FieldIO.H"
#include "SliceDiagnostic.H"
#include "Utils/CoarsenIO.H"

#ifdef WARPX_USE_OPENPMD
#   include "Diagnostics/WarpXOpenPMD.H"
#endif

#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil_F.H>
#include <AMReX_buildInfo.H>

#ifdef BL_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif

#ifdef AMREX_USE_ASCENT
#   include <ascent.hpp>
#   include <AMReX_Conduit_Blueprint.H>
#endif


using namespace amrex;

namespace
{
    const std::string level_prefix {"Level_"};
}

void
WarpX::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void
WarpX::InitFromCheckpoint ()
{
    WARPX_PROFILE("WarpX::InitFromCheckpoint()");

    amrex::Print() << "  Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    {
        std::string File(restart_chkfile + "/WarpXHeader");

        VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        std::string line, word;

        std::getline(is, line);

        int nlevs;
        is >> nlevs;
        GotoNextLine(is);
        finest_level = nlevs-1;

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
                istep[i++] = std::stoi(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
                nsubsteps[i++] = std::stoi(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
                t_new[i++] = std::stod(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
                t_old[i++] = std::stod(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
                dt[i++] = std::stod(word);
            }
        }

        is >> moving_window_x;
        GotoNextLine(is);

        is >> is_synchronized;
        GotoNextLine(is);

        Real prob_lo[AMREX_SPACEDIM];
        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
                prob_lo[i++] = std::stod(word);
            }
        }

        Real prob_hi[AMREX_SPACEDIM];
        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
                prob_hi[i++] = std::stod(word);
            }
        }

        ResetProbDomain(RealBox(prob_lo,prob_hi));

        for (int lev = 0; lev < nlevs; ++lev) {
            BoxArray ba;
            ba.readFrom(is);
            GotoNextLine(is);
            DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
            SetBoxArray(lev, ba);
            SetDistributionMap(lev, dm);
            AllocLevelData(lev, ba, dm);
        }

        mypc->ReadHeader(is);
    }

    const int nlevs = finestLevel()+1;

    // Initialize the field data
    for (int lev = 0; lev < nlevs; ++lev)
    {
        for (int i = 0; i < 3; ++i) {
            current_fp[lev][i]->setVal(0.0);
            Efield_fp[lev][i]->setVal(0.0);
            Bfield_fp[lev][i]->setVal(0.0);
        }

        if (lev > 0) {
            for (int i = 0; i < 3; ++i) {
                Efield_aux[lev][i]->setVal(0.0);
                Bfield_aux[lev][i]->setVal(0.0);

                current_cp[lev][i]->setVal(0.0);
                Efield_cp[lev][i]->setVal(0.0);
                Bfield_cp[lev][i]->setVal(0.0);
            }
        }

        VisMF::Read(*Efield_fp[lev][0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex_fp"));
        VisMF::Read(*Efield_fp[lev][1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey_fp"));
        VisMF::Read(*Efield_fp[lev][2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez_fp"));

        VisMF::Read(*Bfield_fp[lev][0],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx_fp"));
        VisMF::Read(*Bfield_fp[lev][1],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By_fp"));
        VisMF::Read(*Bfield_fp[lev][2],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz_fp"));

        if (is_synchronized) {
            VisMF::Read(*current_fp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jx_fp"));
            VisMF::Read(*current_fp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jy_fp"));
            VisMF::Read(*current_fp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jz_fp"));
        }

        if (lev > 0)
        {
            VisMF::Read(*Efield_cp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex_cp"));
            VisMF::Read(*Efield_cp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey_cp"));
            VisMF::Read(*Efield_cp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez_cp"));

            VisMF::Read(*Bfield_cp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx_cp"));
            VisMF::Read(*Bfield_cp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By_cp"));
            VisMF::Read(*Bfield_cp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz_cp"));

            if (is_synchronized) {
                VisMF::Read(*current_cp[lev][0],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jx_cp"));
                VisMF::Read(*current_cp[lev][1],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jy_cp"));
                VisMF::Read(*current_cp[lev][2],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "jz_cp"));
            }
        }
    }

    if (do_pml)
    {
        InitPML();
        for (int lev = 0; lev < nlevs; ++lev) {
            pml[lev]->Restart(amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml"));
        }
    }

    // Initilize particles
    mypc->AllocData();
    mypc->Restart(restart_chkfile);

}


std::unique_ptr<MultiFab>
WarpX::GetCellCenteredData() {

    WARPX_PROFILE("WarpX::GetCellCenteredData");

    const int ng =  1;
    const int nc = 10;

    Vector<std::unique_ptr<MultiFab> > cc(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        cc[lev].reset( new MultiFab(grids[lev], dmap[lev], nc, ng) );

        int dcomp = 0;
        // first the electric field
        AverageAndPackVectorField( *cc[lev], Efield_aux[lev], dmap[lev], dcomp, ng );
        dcomp += 3;
        // then the magnetic field
        AverageAndPackVectorField( *cc[lev], Bfield_aux[lev], dmap[lev], dcomp, ng );
        dcomp += 3;
        // then the current density
        AverageAndPackVectorField( *cc[lev], current_fp[lev], dmap[lev], dcomp, ng );
        dcomp += 3;
        // then the charge density
        const std::unique_ptr<MultiFab>& charge_density = mypc->GetChargeDensity(lev);
        AverageAndPackScalarField( *cc[lev], *charge_density, dmap[lev], dcomp, ng );

        cc[lev]->FillBoundary(geom[lev].periodicity());
    }

    for (int lev = finest_level; lev > 0; --lev)
    {
        CoarsenIO::Coarsen( *cc[lev-1], *cc[lev], 0, 0, nc, 0, refRatio(lev-1) );
    }

    return std::move(cc[0]);
}

void
WarpX::UpdateInSitu () const
{
#if defined(BL_USE_SENSEI_INSITU) || defined(AMREX_USE_ASCENT)
    WARPX_PROFILE("WarpX::UpdateInSitu()");

    // Average the fields from the simulation to the cell centers
    const int ngrow = 1;
    Vector<std::string> varnames; // Name of the written fields
    // mf_avg will contain the averaged, cell-centered fields
    Vector<MultiFab> mf_avg;
    WarpX::AverageAndPackFields( varnames, mf_avg, ngrow );

#ifdef BL_USE_SENSEI_INSITU
    if (insitu_bridge->update(istep[0], t_new[0],
        dynamic_cast<amrex::AmrMesh*>(const_cast<WarpX*>(this)),
        {&mf_avg}, {varnames}))
    {
        amrex::ErrorStream()
            << "WarpXIO::UpdateInSitu : Failed to update the in situ bridge."
            << std::endl;

        amrex::Abort();
    }
#endif

#ifdef AMREX_USE_ASCENT
    // wrap mesh data
    conduit::Node bp_mesh;
    MultiLevelToBlueprint(finest_level+1,
            amrex::GetVecOfConstPtrs(mf_avg),
            varnames,
            Geom(),
            t_new[0],
            istep,
            refRatio(),
            bp_mesh);

    // wrap particle data for each species
    // we prefix the fields with "particle_{species_name}" b/c we
    // want to to uniquely name all the fields that can be plotted

    std::vector<std::string> species_names = mypc->GetSpeciesNames();

    for (unsigned i = 0, n = species_names.size(); i < n; ++i)
    {

        Vector<std::string> particle_varnames;
        Vector<std::string> particle_int_varnames;
        std::string prefix = "particle_" + species_names[i];

        // Get pc for species
        auto& pc = mypc->GetParticleContainer(i);

        // get names of real comps
        std::map<std::string, int> real_comps_map = pc.getParticleComps();
        std::map<std::string, int>::const_iterator r_itr = real_comps_map.begin();

        // TODO: Looking at other code paths, I am not sure compile time
        //  QED field is included in getParticleComps()?
        while (r_itr != real_comps_map.end())
        {
            // get next real particle name
            std::string varname = r_itr->first;
            particle_varnames.push_back(prefix + "_" + varname);
            r_itr++;
        }

        // get names of int comps
        std::map<std::string, int> int_comps_map = pc.getParticleiComps();
        std::map<std::string, int>::const_iterator i_itr = int_comps_map.begin();

        while (i_itr != int_comps_map.end())
        {
            // get next real particle name
            std::string varname = i_itr->first;
            particle_int_varnames.push_back(prefix + "_" + varname);
            i_itr++;
        }

        // wrap pc for current species into a blueprint topology
        amrex::ParticleContainerToBlueprint(pc,
                                            particle_varnames,
                                            particle_int_varnames,
                                            bp_mesh,
                                            prefix);
    }

    // // If you want to save blueprint HDF5 files w/o using an Ascent
    // // extract, you can call the following AMReX helper:
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
#endif

#endif
}
