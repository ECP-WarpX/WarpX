/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Axel Huebl
 * Burlen Loring, David Grote, Gunther H. Weber
 * Junmin Gu, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BoundaryConditions/PML.H"
#include "FieldIO.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/CoarsenIO.H"
#include "Parallelization/WarpXCommUtil.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#ifdef AMREX_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

#include <array>
#include <istream>
#include <memory>
#include <string>
#include <utility>

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

amrex::DistributionMapping
WarpX::GetRestartDMap (const std::string& chkfile, const amrex::BoxArray& ba, int lev) const {
    std::string DMFileName = chkfile;
    if (!DMFileName.empty() && DMFileName[DMFileName.size()-1] != '/') {DMFileName += '/';}
    DMFileName = amrex::Concatenate(DMFileName + "Level_", lev, 1);
    DMFileName += "/DM";

    if (!amrex::FileExists(DMFileName)) {
        return amrex::DistributionMapping{ba, ParallelDescriptor::NProcs()};
    }

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(DMFileName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream DMFile(fileCharPtrString, std::istringstream::in);
    if ( ! DMFile.good()) amrex::FileOpenFailed(DMFileName);
    DMFile.exceptions(std::ios_base::failbit | std::ios_base::badbit);

    int nprocs_in_checkpoint;
    DMFile >> nprocs_in_checkpoint;
    if (nprocs_in_checkpoint != ParallelDescriptor::NProcs()) {
        return amrex::DistributionMapping{ba, ParallelDescriptor::NProcs()};
    }

    amrex::DistributionMapping dm;
    dm.readFrom(DMFile);
    if (dm.size() != ba.size()) {
        return amrex::DistributionMapping{ba, ParallelDescriptor::NProcs()};
    }

    return dm;
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
        is.exceptions(std::ios_base::failbit | std::ios_base::badbit);

        std::string line, word;

        std::getline(is, line);

        int nlevs;
        is >> nlevs;
        GotoNextLine(is);
        finest_level = nlevs-1;

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (int i = 0; i < istep.size(); ++i) {
                lis >> word;
                istep.at(i) = std::stoi(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (int i = 0; i < nsubsteps.size(); ++i) {
                lis >> word;
                nsubsteps.at(i) = std::stoi(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (int i = 0; i < t_new.size(); ++i) {
                lis >> word;
                t_new.at(i) = static_cast<Real>(std::stod(word));
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (int i = 0; i < t_old.size(); ++i) {
                lis >> word;
                t_old.at(i) = static_cast<Real>(std::stod(word));
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (int i = 0; i < dt.size(); ++i) {
                lis >> word;
                dt.at(i) = static_cast<Real>(std::stod(word));
            }
        }

        amrex::Real moving_window_x_checkpoint;
        is >> moving_window_x_checkpoint;
        GotoNextLine(is);

        is >> is_synchronized;
        GotoNextLine(is);

        amrex::Vector<amrex::Real> prob_lo( AMREX_SPACEDIM );
        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (int i = 0; i < prob_lo.size(); ++i) {
                lis >> word;
                prob_lo.at(i) = static_cast<Real>(std::stod(word));
            }
        }

        amrex::Vector<amrex::Real> prob_hi( AMREX_SPACEDIM );
        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (int i = 0; i < prob_hi.size(); ++i) {
                lis >> word;
                prob_hi.at(i) = static_cast<Real>(std::stod(word));
            }
        }

        ResetProbDomain(RealBox(prob_lo.data(),prob_hi.data()));

        for (int lev = 0; lev < nlevs; ++lev) {
            BoxArray ba;
            ba.readFrom(is);
            GotoNextLine(is);
            DistributionMapping dm = GetRestartDMap(restart_chkfile, ba, lev);
            SetBoxArray(lev, ba);
            SetDistributionMap(lev, dm);
            AllocLevelData(lev, ba, dm);
        }

        mypc->ReadHeader(is);
        is >> current_injection_position;
        GotoNextLine(is);

        int do_moving_window_before_restart;
        is >> do_moving_window_before_restart;
        GotoNextLine(is);

        if (do_moving_window_before_restart) {
            moving_window_x = moving_window_x_checkpoint;
        }

        is >> time_of_last_gal_shift;
        GotoNextLine(is);
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

        if (WarpX::fft_do_time_averaging)
        {
            VisMF::Read(*Efield_avg_fp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex_avg_fp"));
            VisMF::Read(*Efield_avg_fp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey_avg_fp"));
            VisMF::Read(*Efield_avg_fp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez_avg_fp"));

            VisMF::Read(*Bfield_avg_fp[lev][0],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx_avg_fp"));
            VisMF::Read(*Bfield_avg_fp[lev][1],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By_avg_fp"));
            VisMF::Read(*Bfield_avg_fp[lev][2],
                        amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz_avg_fp"));
        }

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

            if (WarpX::fft_do_time_averaging)
            {
                VisMF::Read(*Efield_avg_cp[lev][0],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ex_avg_cp"));
                VisMF::Read(*Efield_avg_cp[lev][1],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ey_avg_cp"));
                VisMF::Read(*Efield_avg_cp[lev][2],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Ez_avg_cp"));

                VisMF::Read(*Bfield_avg_cp[lev][0],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bx_avg_cp"));
                VisMF::Read(*Bfield_avg_cp[lev][1],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "By_avg_cp"));
                VisMF::Read(*Bfield_avg_cp[lev][2],
                            amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "Bz_avg_cp"));
            }

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

    InitPML();
    if (do_pml)
    {
        for (int lev = 0; lev < nlevs; ++lev) {
            pml[lev]->Restart(amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml"));
        }
    }

    // Initialize particles
    mypc->AllocData();
    mypc->Restart(restart_chkfile);

}


std::unique_ptr<MultiFab>
WarpX::GetCellCenteredData() {

    WARPX_PROFILE("WarpX::GetCellCenteredData()");

    const amrex::IntVect ng(1);
    const int nc = 10;

    Vector<std::unique_ptr<MultiFab> > cc(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        cc[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev], nc, ng );

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

        WarpXCommUtil::FillBoundary(*cc[lev], geom[lev].periodicity());
    }

    for (int lev = finest_level; lev > 0; --lev)
    {
        CoarsenIO::Coarsen( *cc[lev-1], *cc[lev], 0, 0, nc, 0, refRatio(lev-1) );
    }

    return std::move(cc[0]);
}
