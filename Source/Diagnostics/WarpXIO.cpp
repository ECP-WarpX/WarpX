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
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
#    include "BoundaryConditions/PML_RZ.H"
#endif
#include "FieldIO.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"
#include "Diagnostics/MultiDiagnostics.H"

#include <ablastr/utils/Communication.H>

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

    amrex::Print()<< Utils::TextMsg::Info(
        "restart from checkpoint " + restart_chkfile);

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
            for (auto& istep_lev : istep) {
                lis >> word;
                istep_lev = std::stoi(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (auto& nsub : nsubsteps) {
                lis >> word;
                nsub = std::stoi(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (auto& t_new_lev : t_new) {
                lis >> word;
                t_new_lev = static_cast<Real>(std::stod(word));
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (auto& t_old_lev : t_old) {
                lis >> word;
                t_old_lev = static_cast<Real>(std::stod(word));
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (auto& dt_lev : dt) {
                lis >> word;
                dt_lev = static_cast<Real>(std::stod(word));
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
            for (auto& prob_lo_comp : prob_lo) {
                lis >> word;
                prob_lo_comp = static_cast<Real>(std::stod(word));
            }
        }

        amrex::Vector<amrex::Real> prob_hi( AMREX_SPACEDIM );
        std::getline(is, line);
        {
            std::istringstream lis(line);
            lis.exceptions(std::ios_base::failbit | std::ios_base::badbit);
            for (auto& prob_hi_comp : prob_hi) {
                lis >> word;
                prob_hi_comp = static_cast<Real>(std::stod(word));
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


        auto & warpx = WarpX::GetInstance();
        for (int idiag = 0; idiag < warpx.GetMultiDiags().GetTotalDiags(); ++idiag)
        {
            if( warpx.GetMultiDiags().diagstypes(idiag) == DiagTypes::BackTransformed )
            {
                auto& diag = warpx.GetMultiDiags().GetDiag(idiag);
                if (diag.getnumbuffers() > 0) {
                    diag.InitDataBeforeRestart();
                    for (int i_buffer=0; i_buffer<diag.getnumbuffers(); ++i_buffer){
                        amrex::Real tlab;
                        is >> tlab;
                        diag.settlab(i_buffer, tlab);
                        int kindex_hi;
                        is >> kindex_hi;
                        diag.set_buffer_k_index_hi(i_buffer, kindex_hi);

                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            amrex::Real snapshot_lo;
                            is >> snapshot_lo;
                            diag.setSnapshotDomainLo(i_buffer, idim, snapshot_lo);
                        }
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            amrex::Real snapshot_hi;
                            is >> snapshot_hi;
                            diag.setSnapshotDomainHi(i_buffer, idim, snapshot_hi);
                        }

                        int flush_counter;
                        is >> flush_counter;
                        diag.set_flush_counter(i_buffer, flush_counter);

                        int last_valid_Zslice;
                        is >> last_valid_Zslice;
                        diag.set_last_valid_Zslice(i_buffer, last_valid_Zslice);

                        int snapshot_full_flag;
                        is >> snapshot_full_flag;
                        diag.set_snapshot_full(i_buffer, snapshot_full_flag);

                    }
                    diag.InitDataAfterRestart();
                } else {
                    diag.InitData();
                }
            } else {
                warpx.GetMultiDiags().GetDiag(idiag).InitData();
            }
        }
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
            if (pml[lev])
                pml[lev]->Restart(amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml"));
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
            if (pml_rz[lev])
                pml_rz[lev]->Restart(amrex::MultiFabFileFullPrefix(lev, restart_chkfile, level_prefix, "pml_rz"));
#endif
        }
    }

    InitializeEBGridData(maxLevel());

    // Initialize particles
    mypc->AllocData();
    mypc->Restart(restart_chkfile);

}
