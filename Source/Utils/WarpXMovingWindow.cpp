/* Copyright 2019-2020 Andrew Myers, Axel Huebl, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
#   include "BoundaryConditions/PML_RZ.H"
#endif
#include "Particles/MultiParticleContainer.H"
#include "Parallelization/WarpXCommUtil.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_Dim3.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_INT.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <array>
#include <cmath>
#include <memory>
#include <string>

using namespace amrex;

void
WarpX::UpdatePlasmaInjectionPosition (amrex::Real a_dt)
{
    int dir = moving_window_dir;
    // Continuously inject plasma in new cells (by default only on level 0)
    if (WarpX::warpx_do_continuous_injection and (WarpX::gamma_boost > 1)){
        // In boosted-frame simulations, the plasma has moved since the last
        // call to this function, and injection position needs to be updated
        current_injection_position -= WarpX::beta_boost *
#if defined(WARPX_DIM_3D)
            WarpX::boost_direction[dir] * PhysConst::c * a_dt;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            // In 2D, dir=0 corresponds to x and dir=1 corresponds to z
            // This needs to be converted in order to index `boost_direction`
            // which has 3 components, for both 2D and 3D simulations.
            WarpX::boost_direction[2*dir] * PhysConst::c * a_dt;
#elif defined(WARPX_DIM_1D_Z)
            // In 1D, dir=0 corresponds to z
            // This needs to be converted in order to index `boost_direction`
            // which has 3 components, for 1D, 2D, and 3D simulations.
            WarpX::boost_direction[2] * PhysConst::c * a_dt;
            amrex::ignore_unused(dir);
#endif
    }
}

int
WarpX::MoveWindow (const int step, bool move_j)
{
    WARPX_PROFILE("WarpX::MoveWindow");

    if (step == start_moving_window_step) {
        amrex::Print() << Utils::TextMsg::Info("Starting moving window");
    }
    if (step == end_moving_window_step) {
        amrex::Print() << Utils::TextMsg::Info("Stopping moving window");
    }
    if (moving_window_active(step) == false) return 0;

    // Update the continuous position of the moving window,
    // and of the plasma injection
    moving_window_x += (moving_window_v - WarpX::beta_boost * PhysConst::c)/(1 - moving_window_v * WarpX::beta_boost / PhysConst::c) * dt[0];
    int dir = moving_window_dir;

    // Update warpx.current_injection_position
    // PhysicalParticleContainer uses this injection position
    UpdatePlasmaInjectionPosition( dt[0] );
    if (WarpX::warpx_do_continuous_injection){
        // Update injection position for WarpXParticleContainer in mypc.
        // Nothing to do for PhysicalParticleContainers
        // For LaserParticleContainer, need to update the antenna position.
        mypc->UpdateContinuousInjectionPosition( dt[0] );
    }

    // compute the number of cells to shift on the base level
    amrex::Real new_lo[AMREX_SPACEDIM];
    amrex::Real new_hi[AMREX_SPACEDIM];
    const amrex::Real* current_lo = geom[0].ProbLo();
    const amrex::Real* current_hi = geom[0].ProbHi();
    const amrex::Real* cdx = geom[0].CellSize();
    int num_shift_base = static_cast<int>((moving_window_x - current_lo[dir]) / cdx[dir]);

    if (num_shift_base == 0) return 0;

    // update the problem domain. Note the we only do this on the base level because
    // amrex::Geometry objects share the same, static RealBox.
    for (int i=0; i<AMREX_SPACEDIM; i++) {
        new_lo[i] = current_lo[i];
        new_hi[i] = current_hi[i];
    }
    new_lo[dir] = current_lo[dir] + num_shift_base * cdx[dir];
    new_hi[dir] = current_hi[dir] + num_shift_base * cdx[dir];

    ResetProbDomain(amrex::RealBox(new_lo, new_hi));

    // Moving slice coordinates - lo and hi - with moving window //
    // slice box is modified only if slice diagnostics is initialized in input //
    if ( slice_plot_int > 0 )
    {
       amrex::Real new_slice_lo[AMREX_SPACEDIM];
       amrex::Real new_slice_hi[AMREX_SPACEDIM];
       const amrex::Real* current_slice_lo = slice_realbox.lo();
       const amrex::Real* current_slice_hi = slice_realbox.hi();
       for ( int i = 0; i < AMREX_SPACEDIM; i++) {
           new_slice_lo[i] = current_slice_lo[i];
           new_slice_hi[i] = current_slice_hi[i];
       }
       int num_shift_base_slice = static_cast<int> ((moving_window_x -
                                  current_slice_lo[dir]) / cdx[dir]);
       new_slice_lo[dir] = current_slice_lo[dir] + num_shift_base_slice*cdx[dir];
       new_slice_hi[dir] = current_slice_hi[dir] + num_shift_base_slice*cdx[dir];
       slice_realbox.setLo(new_slice_lo);
       slice_realbox.setHi(new_slice_hi);
    }

    int num_shift      = num_shift_base;
    int num_shift_crse = num_shift;

    // Shift the mesh fields
    for (int lev = 0; lev <= finest_level; ++lev) {

        if (lev > 0) {
            num_shift_crse = num_shift;
            num_shift *= refRatio(lev-1)[dir];
        }

        // Shift each component of vector fields (E, B, j)
        for (int dim = 0; dim < 3; ++dim) {
            // Fine grid
            amrex::ParserExecutor<3> Bfield_parser;
            amrex::ParserExecutor<3> Efield_parser;
            bool use_Bparser = false;
            bool use_Eparser = false;
            if (B_ext_grid_s == "parse_b_ext_grid_function") {
                use_Bparser = true;
                if (dim == 0) Bfield_parser = Bxfield_parser->compile<3>();
                if (dim == 1) Bfield_parser = Byfield_parser->compile<3>();
                if (dim == 2) Bfield_parser = Bzfield_parser->compile<3>();
            }
            if (E_ext_grid_s == "parse_e_ext_grid_function") {
                use_Eparser = true;
                if (dim == 0) Efield_parser = Exfield_parser->compile<3>();
                if (dim == 1) Efield_parser = Eyfield_parser->compile<3>();
                if (dim == 2) Efield_parser = Ezfield_parser->compile<3>();
            }
            shiftMF(*Bfield_fp[lev][dim], geom[lev], num_shift, dir, lev, B_external_grid[dim], use_Bparser, Bfield_parser);
            shiftMF(*Efield_fp[lev][dim], geom[lev], num_shift, dir, lev, E_external_grid[dim], use_Eparser, Efield_parser);
            if (fft_do_time_averaging) {
                shiftMF(*Bfield_avg_fp[lev][dim], geom[lev], num_shift, dir, lev, B_external_grid[dim], use_Bparser, Bfield_parser);
                shiftMF(*Efield_avg_fp[lev][dim], geom[lev], num_shift, dir, lev, E_external_grid[dim], use_Eparser, Efield_parser);
            }
            if (move_j) {
                shiftMF(*current_fp[lev][dim], geom[lev], num_shift, dir, lev);
            }
            if (pml[lev] && pml[lev]->ok()) {
                const std::array<amrex::MultiFab*, 3>& pml_B = pml[lev]->GetB_fp();
                const std::array<amrex::MultiFab*, 3>& pml_E = pml[lev]->GetE_fp();
                shiftMF(*pml_B[dim], geom[lev], num_shift, dir, lev);
                shiftMF(*pml_E[dim], geom[lev], num_shift, dir, lev);
            }
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
            if (pml_rz[lev] && dim < 2) {
                const std::array<amrex::MultiFab*, 2>& pml_rz_B = pml_rz[lev]->GetB_fp();
                const std::array<amrex::MultiFab*, 2>& pml_rz_E = pml_rz[lev]->GetE_fp();
                shiftMF(*pml_rz_B[dim], geom[lev], num_shift, dir, lev);
                shiftMF(*pml_rz_E[dim], geom[lev], num_shift, dir, lev);
            }
#endif
            if (lev > 0) {
                // coarse grid
                shiftMF(*Bfield_cp[lev][dim], geom[lev-1], num_shift_crse, dir, lev, B_external_grid[dim], use_Bparser, Bfield_parser);
                shiftMF(*Efield_cp[lev][dim], geom[lev-1], num_shift_crse, dir, lev, E_external_grid[dim], use_Eparser, Efield_parser);
                shiftMF(*Bfield_aux[lev][dim], geom[lev], num_shift, dir, lev);
                shiftMF(*Efield_aux[lev][dim], geom[lev], num_shift, dir, lev);
                if (fft_do_time_averaging) {
                    shiftMF(*Bfield_avg_cp[lev][dim], geom[lev-1], num_shift_crse, dir, lev, B_external_grid[dim], use_Bparser, Bfield_parser);
                    shiftMF(*Efield_avg_cp[lev][dim], geom[lev-1], num_shift_crse, dir, lev, E_external_grid[dim], use_Eparser, Efield_parser);
                }
                if (move_j) {
                    shiftMF(*current_cp[lev][dim], geom[lev-1], num_shift_crse, dir, lev);
                }
                if (do_pml && pml[lev]->ok()) {
                    const std::array<amrex::MultiFab*, 3>& pml_B = pml[lev]->GetB_cp();
                    const std::array<amrex::MultiFab*, 3>& pml_E = pml[lev]->GetE_cp();
                    shiftMF(*pml_B[dim], geom[lev-1], num_shift_crse, dir, lev);
                    shiftMF(*pml_E[dim], geom[lev-1], num_shift_crse, dir, lev);
                }
            }
        }

        // Shift scalar component F for dive cleaning
        if (do_dive_cleaning) {
            // Fine grid
            shiftMF(*F_fp[lev], geom[lev], num_shift, dir, lev);
            if (do_pml && pml[lev]->ok()) {
                amrex::MultiFab* pml_F = pml[lev]->GetF_fp();
                shiftMF(*pml_F, geom[lev], num_shift, dir, lev);
            }
            if (lev > 0) {
                // Coarse grid
                shiftMF(*F_cp[lev], geom[lev-1], num_shift_crse, dir, lev);
                if (do_pml && pml[lev]->ok()) {
                    amrex::MultiFab* pml_F = pml[lev]->GetF_cp();
                    shiftMF(*pml_F, geom[lev-1], num_shift_crse, dir, lev);
                }
            }
        }

        // Shift scalar component rho
        if (move_j) {
            if (rho_fp[lev]){
                // Fine grid
                shiftMF(*rho_fp[lev],   geom[lev], num_shift, dir, lev);
                if (lev > 0){
                    // Coarse grid
                    shiftMF(*rho_cp[lev], geom[lev-1], num_shift_crse, dir, lev);
                }
            }
        }
    }

    // Continuously inject plasma in new cells (by default only on level 0)
    if (WarpX::warpx_do_continuous_injection) {

        const int lev = 0;

        // particleBox encloses the cells where we generate particles
        // (only injects particles in an integer number of cells,
        // for correct particle spacing)
        amrex::RealBox particleBox = geom[lev].ProbDomain();
        amrex::Real new_injection_position;
        if (moving_window_v >= 0){
            // Forward-moving window
            amrex::Real dx = geom[lev].CellSize(dir);
            new_injection_position = current_injection_position +
                std::floor( (geom[lev].ProbHi(dir) - current_injection_position)/dx ) * dx;
        } else {
            // Backward-moving window
            amrex::Real dx = geom[lev].CellSize(dir);
            new_injection_position = current_injection_position -
                std::floor( (current_injection_position - geom[lev].ProbLo(dir))/dx) * dx;
        }
        // Modify the corresponding bounds of the particleBox
        if (moving_window_v >= 0) {
            particleBox.setLo( dir, current_injection_position );
            particleBox.setHi( dir, new_injection_position );
        } else {
            particleBox.setLo( dir, new_injection_position );
            particleBox.setHi( dir, current_injection_position );
        }

        if (particleBox.ok() and (current_injection_position != new_injection_position)){
            // Performs continuous injection of all WarpXParticleContainer
            // in mypc.
            mypc->ContinuousInjection(particleBox);
            current_injection_position = new_injection_position;
        }
    }

    return num_shift_base;
}

void
WarpX::shiftMF (amrex::MultiFab& mf, const amrex::Geometry& geom,
                int num_shift, int dir, const int lev,
                amrex::Real external_field, bool useparser,
                amrex::ParserExecutor<3> const& field_parser)
{
    using namespace amrex::literals;
    WARPX_PROFILE("WarpX::shiftMF()");
    const amrex::BoxArray& ba = mf.boxArray();
    const amrex::DistributionMapping& dm = mf.DistributionMap();
    const int nc = mf.nComp();
    const amrex::IntVect& ng = mf.nGrowVect();

    AMREX_ALWAYS_ASSERT(ng[dir] >= num_shift);

    amrex::MultiFab tmpmf(ba, dm, nc, ng);
    amrex::MultiFab::Copy(tmpmf, mf, 0, 0, nc, ng);

    if ( WarpX::safe_guard_cells ) {
        // Fill guard cells.
        WarpXCommUtil::FillBoundary(tmpmf, geom.periodicity());
    } else {
        amrex::IntVect ng_mw = amrex::IntVect::TheUnitVector();
        // Enough guard cells in the MW direction
        ng_mw[dir] = num_shift;
        // Make sure we don't exceed number of guard cells allocated
        ng_mw = ng_mw.min(ng);
        // Fill guard cells.
        WarpXCommUtil::FillBoundary(tmpmf, ng_mw, geom.periodicity());
    }

    // Make a box that covers the region that the window moved into
    const amrex::IndexType& typ = ba.ixType();
    const amrex::Box& domainBox = geom.Domain();
    amrex::Box adjBox;
    if (num_shift > 0) {
        adjBox = adjCellHi(domainBox, dir, ng[dir]);
    } else {
        adjBox = adjCellLo(domainBox, dir, ng[dir]);
    }
    adjBox = amrex::convert(adjBox, typ);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim == dir and typ.nodeCentered(dir)) {
            if (num_shift > 0) {
                adjBox.growLo(idim, -1);
            } else {
                adjBox.growHi(idim, -1);
            }
        } else if (idim != dir) {
            adjBox.growLo(idim, ng[idim]);
            adjBox.growHi(idim, ng[idim]);
        }
    }

    amrex::IntVect shiftiv(0);
    shiftiv[dir] = num_shift;
    amrex::Dim3 shift = shiftiv.dim3();

    const amrex::RealBox& real_box = geom.ProbDomain();
    const auto dx = geom.CellSizeArray();

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

    for (amrex::MFIter mfi(tmpmf); mfi.isValid(); ++mfi )
    {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
        }
        amrex::Real wt = amrex::second();

        auto const& dstfab = mf.array(mfi);
        auto const& srcfab = tmpmf.array(mfi);

        const amrex::Box& outbox = mfi.fabbox() & adjBox;

        if (outbox.ok()) {
            if (useparser == false) {
                AMREX_PARALLEL_FOR_4D ( outbox, nc, i, j, k, n,
                {
                    srcfab(i,j,k,n) = external_field;
                })
            } else if (useparser == true) {
                // index type of the src mf
                auto const& mf_IndexType = (tmpmf).ixType();
                amrex::IntVect mf_type(AMREX_D_DECL(0,0,0));
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    mf_type[idim] = mf_IndexType.nodeCentered(idim);
                }

                amrex::ParallelFor (outbox, nc,
                      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                      // Compute x,y,z co-ordinates based on index type of mf
#if defined(WARPX_DIM_1D_Z)
                      amrex::Real x = 0.0_rt;
                      amrex::Real y = 0.0_rt;
                      amrex::Real fac_z = (1.0_rt - mf_type[0]) * dx[0]*0.5_rt;
                      amrex::Real z = i*dx[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                      amrex::Real fac_x = (1.0_rt - mf_type[0]) * dx[0]*0.5_rt;
                      amrex::Real x = i*dx[0] + real_box.lo(0) + fac_x;
                      amrex::Real y = 0.0;
                      amrex::Real fac_z = (1.0_rt - mf_type[1]) * dx[1]*0.5_rt;
                      amrex::Real z = j*dx[1] + real_box.lo(1) + fac_z;
#else
                      amrex::Real fac_x = (1.0_rt - mf_type[0]) * dx[0]*0.5_rt;
                      amrex::Real x = i*dx[0] + real_box.lo(0) + fac_x;
                      amrex::Real fac_y = (1.0_rt - mf_type[1]) * dx[1]*0.5_rt;
                      amrex::Real y = j*dx[1] + real_box.lo(1) + fac_y;
                      amrex::Real fac_z = (1.0_rt - mf_type[2]) * dx[2]*0.5_rt;
                      amrex::Real z = k*dx[2] + real_box.lo(2) + fac_z;
#endif
                      srcfab(i,j,k,n) = field_parser(x,y,z);
                });
            }

        }

        amrex::Box dstBox = mf[mfi].box();
        if (num_shift > 0) {
            dstBox.growHi(dir, -num_shift);
        } else {
            dstBox.growLo(dir,  num_shift);
        }
        AMREX_PARALLEL_FOR_4D ( dstBox, nc, i, j, k, n,
        {
            dstfab(i,j,k,n) = srcfab(i+shift.x,j+shift.y,k+shift.z,n);
        })

        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            amrex::Gpu::synchronize();
            wt = amrex::second() - wt;
            amrex::HostDevice::Atomic::Add( &(*cost)[mfi.index()], wt);
        }
    }

#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
    if (WarpX::GetInstance().getPMLRZ()) {
        // This does the exchange of data in the corner guard cells, the cells that are in the
        // guard region both radially and longitudinally. These are the PML cells in the overlapping
        // longitudinal region. FillBoundary normally does not update these cells.
        // This update is needed so that the cells at the end of the FABs are updated appropriately
        // with the data shifted from the nieghboring FAB. Without this update, the RZ PML becomes
        // unstable with the moving grid.
        // This code creates a temporary MultiFab using a BoxList where the radial size of all of
        // its boxes is increased so that the radial guard cells are included in the boxes valid domain.
        // The temporary MultiFab is setup to refer to the data of the original Multifab (this can
        // be done since the shape of the data is all the same, just the indexing is different).
        amrex::BoxList bl;
        for (int i = 0, N=ba.size(); i < N; ++i) {
            bl.push_back(amrex::grow(ba[i], 0, mf.nGrowVect()[0]));
        }
        amrex::BoxArray rba(std::move(bl));
        amrex::MultiFab rmf(rba, dm, mf.nComp(), IntVect(0,mf.nGrowVect()[1]), MFInfo().SetAlloc(false));

        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            rmf.setFab(mfi, FArrayBox(mf[mfi], amrex::make_alias, 0, mf.nComp()));
        }
        rmf.FillBoundary(false);
    }
#endif

}

void
WarpX::ShiftGalileanBoundary ()
{
    amrex::Real cur_time = t_new[0];
    amrex::Real new_lo[AMREX_SPACEDIM];
    amrex::Real new_hi[AMREX_SPACEDIM];
    const amrex::Real* current_lo = geom[0].ProbLo();
    const amrex::Real* current_hi = geom[0].ProbHi();

    amrex::Real time_shift = (cur_time - time_of_last_gal_shift);

#if defined(WARPX_DIM_3D)
        m_galilean_shift = {
            m_v_galilean[0]*time_shift,
            m_v_galilean[1]*time_shift,
            m_v_galilean[2]*time_shift };
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        m_galilean_shift = {
            m_v_galilean[0]*time_shift,
            std::numeric_limits<amrex::Real>::quiet_NaN(),
            m_v_galilean[2]*time_shift };
#elif defined(WARPX_DIM_1D_Z)
        m_galilean_shift = {
            std::numeric_limits<Real>::quiet_NaN(),
            std::numeric_limits<Real>::quiet_NaN(),
            m_v_galilean[2]*time_shift };
#endif

#if defined(WARPX_DIM_3D)
        for (int i=0; i<AMREX_SPACEDIM; i++) {
            new_lo[i] = current_lo[i] + m_galilean_shift[i];
            new_hi[i] = current_hi[i] + m_galilean_shift[i];
        }
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    {
        new_lo[0] = current_lo[0] + m_galilean_shift[0];
        new_hi[0] = current_hi[0] + m_galilean_shift[0];
        new_lo[1] = current_lo[1] + m_galilean_shift[2];
        new_hi[1] = current_hi[1] + m_galilean_shift[2];
    }
#elif defined(WARPX_DIM_1D_Z)
    {
        new_lo[0] = current_lo[0] + m_galilean_shift[2];
        new_hi[0] = current_hi[0] + m_galilean_shift[2];
    }
#endif
    time_of_last_gal_shift = cur_time;

    ResetProbDomain(amrex::RealBox(new_lo, new_hi));
}

void
WarpX::ResetProbDomain (const amrex::RealBox& rb)
{
    amrex::Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        amrex::Geometry g = Geom(lev);
        g.ProbDomain(rb);
        SetGeometry(lev, g);
    }
}
