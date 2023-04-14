/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Aurore Blelly
 * Axel Huebl, Burlen Loring, Maxence Thevenet
 * Michael Rowan, Remi Lehe, Revathi Jambunathan
 * Weiqun Zhang
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
#include "Diagnostics/MultiDiagnostics.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "FieldSolver/FiniteDifferenceSolver/MacroscopicProperties/MacroscopicProperties.H"
#include "Filter/BilinearFilter.H"
#include "Filter/NCIGodfreyFilter.H"
#include "Particles/MultiParticleContainer.H"
#include "Utils/Logo/GetLogo.H"
#include "Utils/MPIInitHelpers.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"
#include "Python/WarpX_py.H"

#include <ablastr/utils/Communication.H>
#include <ablastr/utils/UsedInputsFile.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#ifdef AMREX_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_INT.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"

using namespace amrex;

void
WarpX::PostProcessBaseGrids (BoxArray& ba0) const
{
    if (numprocs != 0) {
        const Box& dom = Geom(0).Domain();
        const IntVect& domlo = dom.smallEnd();
        const IntVect& domlen = dom.size();
        const IntVect sz = domlen / numprocs;
        const IntVect extra = domlen - sz*numprocs;
        BoxList bl;
#if defined(WARPX_DIM_3D)
        for (int k = 0; k < numprocs[2]; ++k) {
            // The first extra[2] blocks get one extra cell with a total of
            // sz[2]+1.  The rest get sz[2] cells.  The docomposition in y
            // and x directions are similar.
            int klo = (k < extra[2]) ? k*(sz[2]+1) : (k*sz[2]+extra[2]);
            int khi = (k < extra[2]) ? klo+(sz[2]+1)-1 : klo+sz[2]-1;
            klo += domlo[2];
            khi += domlo[2];
#endif
#if (AMREX_SPACEDIM >= 2)
            for (int j = 0; j < numprocs[1]; ++j) {
                int jlo = (j < extra[1]) ? j*(sz[1]+1) : (j*sz[1]+extra[1]);
                int jhi = (j < extra[1]) ? jlo+(sz[1]+1)-1 : jlo+sz[1]-1;
                jlo += domlo[1];
                jhi += domlo[1];
#endif
                for (int i = 0; i < numprocs[0]; ++i) {
                    int ilo = (i < extra[0]) ? i*(sz[0]+1) : (i*sz[0]+extra[0]);
                    int ihi = (i < extra[0]) ? ilo+(sz[0]+1)-1 : ilo+sz[0]-1;
                    ilo += domlo[0];
                    ihi += domlo[0];
                    bl.push_back(Box(IntVect(AMREX_D_DECL(ilo,jlo,klo)),
                                     IntVect(AMREX_D_DECL(ihi,jhi,khi))));
        AMREX_D_TERM(},},})
        ba0 = BoxArray(std::move(bl));
    }
}

void
WarpX::PrintMainPICparameters ()
{
    amrex::Print() << "-------------------------------------------------------------------------------\n";
    amrex::Print() << "--------------------------- MAIN EM PIC PARAMETERS ----------------------------\n";
    amrex::Print() << "-------------------------------------------------------------------------------\n";

    // print warpx build information
    if constexpr (std::is_same<Real, float>::value) {
      amrex::Print() << "Precision:            | SINGLE" << "\n";
    }
    else {
      amrex::Print() << "Precision:            | DOUBLE" << "\n";
    }
    if constexpr (std::is_same<ParticleReal, float>::value) {
      amrex::Print() << "Particle precision:   | SINGLE" << "\n";
    }
    else {
      amrex::Print() << "Particle precision:   | DOUBLE" << "\n";
    }

    // Print geometry dimensionality
    amrex::ParmParse pp_geometry("geometry");
    std::string dims;
    pp_geometry.query( "dims", dims );
    if (dims=="1") {
      amrex::Print() << "Geometry:             | 1D (Z)" << "\n";
    }
    else if (dims=="2") {
      amrex::Print() << "Geometry:             | 2D (XZ)" << "\n";
    }
    else if (dims=="3") {
      amrex::Print() << "Geometry:             | 3D (XYZ)" << "\n";
    }
    else if (dims=="RZ") {
      amrex::Print() << "Geometry:             | 2D (RZ)" << "\n";
    }

    #ifdef WARPX_DIM_RZ
      amrex::Print() << "                      | - n_rz_azimuthal_modes = " <<
                     WarpX::n_rz_azimuthal_modes << "\n";
    #endif // WARPX_USE_RZ
    //Print solver's operation mode (e.g., EM or electrostatic)
    if (electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame) {
      amrex::Print() << "Operation mode:       | Electrostatic" << "\n";
      amrex::Print() << "                      | - laboratory frame" << "\n";
    }
    else if (electrostatic_solver_id == ElectrostaticSolverAlgo::Relativistic){
      amrex::Print() << "Operation mode:       | Electrostatic" << "\n";
      amrex::Print() << "                      | - relativistic" << "\n";
    }
    else if (electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic){
      amrex::Print() << "Operation mode:       | Electrostatic" << "\n";
      amrex::Print() << "                      | - laboratory frame, electrostatic + magnetostatic" << "\n";
    }
    else{
      amrex::Print() << "Operation mode:       | Electromagnetic" << "\n";
    }
    if (em_solver_medium == MediumForEM::Vacuum ){
      amrex::Print() << "                      | - vacuum" << "\n";
    }
    else if (em_solver_medium == MediumForEM::Macroscopic ){
      amrex::Print() << "                      | - macroscopic" << "\n";
    }
    if ( (em_solver_medium == MediumForEM::Macroscopic) &&
       (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::LaxWendroff)){
      amrex::Print() << "                      |  - Lax-Wendroff algorithm\n";
      }
    else if ((em_solver_medium == MediumForEM::Macroscopic) &&
            (WarpX::macroscopic_solver_algo == MacroscopicSolverAlgo::BackwardEuler)){
      amrex::Print() << "                      |  - Backward Euler algorithm\n";
      }
    amrex::Print() << "-------------------------------------------------------------------------------\n";
    // Print type of current deposition
    if (current_deposition_algo == CurrentDepositionAlgo::Direct){
      amrex::Print() << "Current Deposition:   | direct \n";
    }
    else if (current_deposition_algo == CurrentDepositionAlgo::Vay){
      amrex::Print() << "Current Deposition:   | Vay \n";
    }
    else if (current_deposition_algo == CurrentDepositionAlgo::Esirkepov){
      amrex::Print() << "Current Deposition:   | Esirkepov \n";
    }
    // Print type of particle pusher
    if (particle_pusher_algo == ParticlePusherAlgo::Vay){
      amrex::Print() << "Particle Pusher:      | Vay \n";
    }
    else if (particle_pusher_algo == ParticlePusherAlgo::HigueraCary){
      amrex::Print() << "Particle Pusher:      | Higuera-Cary \n";
    }
    else if (particle_pusher_algo == ParticlePusherAlgo::Boris){
      amrex::Print() << "Particle Pusher:      | Boris \n";
    }
    // Print type of charge deposition
    if (charge_deposition_algo == ChargeDepositionAlgo::Standard){
      amrex::Print() << "Charge Deposition:    | standard \n";
    }
    // Print field gathering algorithm
    if (field_gathering_algo == GatheringAlgo::MomentumConserving){
      amrex::Print() << "Field Gathering:      | momentum-conserving \n";
    }
    else{
      amrex::Print() << "Field Gathering:      | energy-conserving \n";
    }
    // Print particle's shape factors
    amrex::Print() << "Particle Shape Factor:| " << WarpX::nox << "\n";
    amrex::Print() << "-------------------------------------------------------------------------------\n";
    // Print solver's type: Yee, CKC, ECT
    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::Yee){
      amrex::Print() << "Maxwell Solver:       | Yee \n";
    }
    else if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::CKC){
      amrex::Print() << "Maxwell Solver:       | CKC \n";
    }
    else if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT){
      amrex::Print() << "Maxwell Solver:       | ECT \n";
    }
  #ifdef WARPX_USE_PSATD
    // Print PSATD solver's configuration
    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD){
      amrex::Print() << "Maxwell Solver:       | PSATD \n";
      }
    if ((m_v_galilean[0]!=0) or (m_v_galilean[1]!=0) or (m_v_galilean[2]!=0)) {
      amrex::Print() << "                      | - Galilean \n" <<
      "                      |  - v_galilean = (" << m_v_galilean[0] << "," <<
                              m_v_galilean[1] << "," << m_v_galilean[2] << ")\n";
      }
    if ((m_v_comoving[0]!=0) or (m_v_comoving[1]!=0) or (m_v_comoving[2]!=0)) {
      amrex::Print() << "                      | - comoving \n" <<
      "                      |  - v_comoving = (" << m_v_comoving[0] << "," <<
                              m_v_comoving[1] << "," << m_v_comoving[2] << ")\n";
      }
    if (WarpX::update_with_rho==1) {
      amrex::Print() << "                      | - update with rho is ON \n";
      }
    if (current_correction==1) {
      amrex::Print() << "                      | - current correction is ON \n";
        }
    if (WarpX::do_dive_cleaning==1) {
      amrex::Print() << "                      | - div(E) cleaning is ON \n";
      }
    if (WarpX::do_divb_cleaning==1) {
      amrex::Print() << "                      | - div(B) cleaning is ON \n";
      }
    if (do_multi_J == 1){
      amrex::Print() << "                      | - multi-J deposition is ON \n";
      amrex::Print() << "                      |   - do_multi_J_n_depositions = "
                                        << WarpX::do_multi_J_n_depositions << "\n";
    }
    if (fft_do_time_averaging == 1){
      amrex::Print()<<"                      | - time-averaged is ON \n";
    }
  #endif // WARPX_USE_PSATD

  if (grid_type == GridType::Collocated){
    amrex::Print() << "                      | - collocated grid \n";
  }
  #ifdef WARPX_USE_PSATD
    if ( (grid_type == GridType::Staggered) && (field_gathering_algo == GatheringAlgo::EnergyConserving) ){
      amrex::Print()<<"                      | - staggered grid " << "\n";
    }
    else if ( (grid_type == GridType::Hybrid) && (field_gathering_algo == GatheringAlgo::MomentumConserving) ){
    amrex::Print()<<"                      | - hybrid grid " << "\n";
    if (dims=="3"){
      amrex::Print() << "                      |   - field_centering_nox = " << WarpX::field_centering_nox << "\n";
      amrex::Print() << "                      |   - field_centering_noy = " << WarpX::field_centering_noy << "\n";
      amrex::Print() << "                      |   - field_centering_noz = " << WarpX::field_centering_noz << "\n";
      amrex::Print() << "                      |   - current_centering_nox = " << WarpX::current_centering_nox << "\n";
      amrex::Print() << "                      |   - current_centering_noy = " << WarpX::current_centering_noy << "\n";
      amrex::Print() << "                      |   - current_centering_noz = " << WarpX::current_centering_noz << "\n";
    }
    else if (dims=="2"){
      amrex::Print() << "                      |   - field_centering_nox = " << WarpX::field_centering_nox << "\n";
      amrex::Print() << "                      |   - field_centering_noz = " << WarpX::field_centering_noz << "\n";
      amrex::Print() << "                      |   - current_centering_nox = " << WarpX::current_centering_nox << "\n";
      amrex::Print() << "                      |   - current_centering_noz = " << WarpX::current_centering_noz << "\n";
     }
    else if (dims=="1"){
      amrex::Print() << "                      |   - field_centering_noz = " << WarpX::field_centering_noz << "\n";
      amrex::Print() << "                      |   - current_centering_noz = " << WarpX::current_centering_noz << "\n";
     }
    }
    if (WarpX::use_hybrid_QED == true){
      amrex::Print() << "                      | - use_hybrid_QED = true \n";
    }

    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD){
    // Print solver's order
      std::string psatd_nox_fft, psatd_noy_fft, psatd_noz_fft;
      psatd_nox_fft = (nox_fft == -1) ? "inf" : std::to_string(nox_fft);
      psatd_noy_fft = (noy_fft == -1) ? "inf" : std::to_string(noy_fft);
      psatd_noz_fft = (noz_fft == -1) ? "inf" : std::to_string(noz_fft);

      if (dims=="3" ){
        amrex::Print() << "Spectral order:       | - psatd.nox = " << psatd_nox_fft << "\n";
        amrex::Print() << "                      | - psatd.noy = " << psatd_noy_fft << "\n";
        amrex::Print() << "                      | - psatd.noz = " << psatd_noz_fft << "\n";
      }
      else if (dims=="2" and WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD){
        amrex::Print() << "Spectral order:       | - psatd.nox = " << psatd_nox_fft << "\n";
        amrex::Print() << "                      | - psatd.noz = " << psatd_noz_fft << "\n";
      }
      else if (dims=="1" and WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD){
        amrex::Print() << "Spectral order:       | - psatd.noz = " << psatd_noz_fft << "\n";
      }
    }
    // Print guard cells number
    amrex::Print() << "Guard cells           | - ng_alloc_EB = " << guard_cells.ng_alloc_EB << "\n";
    amrex::Print() << " (allocated for E/B)  | \n";

    #endif // WARPX_USE_PSATD
    amrex::Print() << "-------------------------------------------------------------------------------" << "\n";
    //Print main boosted frame algorithm's parameters
    if (WarpX::gamma_boost!=1){
    amrex::Print() << "Boosted Frame:        |    ON  \n";
    amrex::Print() << "                      |  - gamma_boost = " << WarpX::gamma_boost << "\n";
    amrex::Print() << "                      |  - boost_direction = (" << WarpX::boost_direction[0] <<
                             "," << WarpX::boost_direction[1] << "," << WarpX::boost_direction[2] << ")\n";
    amrex::Print() << "------------------------------------------------------------------------------- \n";
    }
    //Print moving window details
    if (WarpX::do_moving_window == 1){
      amrex::Print() << "Moving window:        |    ON  \n";
      if (WarpX::moving_window_dir == 0){
        amrex::Print() << "                      |  - moving_window_dir = x \n";
      }
      #if defined(WARPX_DIM_3D)
      else if (WarpX::moving_window_dir == 1){
        amrex::Print() << "                      |  - moving_window_dir = y \n";
      }
      #endif
      else if (WarpX::moving_window_dir == WARPX_ZINDEX) {
        amrex::Print() << "                      |  - moving_window_dir = z \n";
      }
      amrex::Print() << "                      |  - moving_window_v = " << WarpX::moving_window_v << "\n";
      amrex::Print() << "------------------------------------------------------------------------------- \n";
    }
}

void
WarpX::WriteUsedInputsFile () const
{
    std::string filename = "warpx_used_inputs";
    ParmParse pp_warpx("warpx");
    pp_warpx.queryAdd("used_inputs_file", filename);

    ablastr::utils::write_used_inputs_file(filename);
}

void
WarpX::InitData ()
{
    WARPX_PROFILE("WarpX::InitData()");
    utils::warpx_check_mpi_thread_level();

#ifdef WARPX_QED
    Print() << "PICSAR (" << WarpX::PicsarVersion() << ")\n";
#endif

    Print() << "WarpX (" << WarpX::Version() << ")\n";

    Print() << utils::logo::get_logo();

    // WarpX::computeMaxStepBoostAccelerator
    // needs to start from the initial zmin_domain_boost,
    // even if restarting from a checkpoint file
    if (do_compute_max_step_from_zmax) {
        zmin_domain_boost_step_0 = geom[0].ProbLo(WARPX_ZINDEX);
    }
    if (restart_chkfile.empty())
    {
        ComputeDt();
        WarpX::PrintDtDxDyDz();
        InitFromScratch();
        InitDiagnostics();
    }
    else
    {
        InitFromCheckpoint();
        WarpX::PrintDtDxDyDz();
        PostRestart();
        reduced_diags->InitData();
    }

    ComputeMaxStep();

    ComputePMLFactors();

    if (WarpX::use_fdtd_nci_corr) {
        WarpX::InitNCICorrector();
    }

    BuildBufferMasks();

    if (WarpX::em_solver_medium==1) {
        m_macroscopic_properties->InitData();
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\nGrids Summary:\n";
        printGridSummary(std::cout, 0, finestLevel());
    }

    // Check that the number of guard cells is smaller than the number of valid cells for all MultiFabs
    // (example: a box with 16 valid cells and 32 guard cells in z will not be considered valid)
    CheckGuardCells();

    PrintMainPICparameters();
    WriteUsedInputsFile();

    if (restart_chkfile.empty())
    {
        // Loop through species and calculate their space-charge field
        bool const reset_fields = false; // Do not erase previous user-specified values on the grid
        ExecutePythonCallback("beforeInitEsolve");
        ComputeSpaceChargeField(reset_fields);
        ExecutePythonCallback("afterInitEsolve");
        if (electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic)
            ComputeMagnetostaticField();

        // Write full diagnostics before the first iteration.
        multi_diags->FilterComputePackFlush( -1 );

        // Write reduced diagnostics before the first iteration.
        if (reduced_diags->m_plot_rd != 0)
        {
            reduced_diags->ComputeDiags(-1);
            reduced_diags->WriteToFile(-1);
        }
    }

    PerformanceHints();

    CheckKnownIssues();
}

void
WarpX::InitDiagnostics () {
    multi_diags->InitData();
    reduced_diags->InitData();
}

void
WarpX::InitFromScratch ()
{
    const Real time = 0.0;

    AmrCore::InitFromScratch(time);  // This will call MakeNewLevelFromScratch

    mypc->AllocData();
    mypc->InitData();

    InitPML();
}

void
WarpX::InitPML ()
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (WarpX::field_boundary_lo[idim] == FieldBoundaryType::PML) {
            do_pml = 1;
            do_pml_Lo[0][idim] = 1; // on level 0
        }
        if (WarpX::field_boundary_hi[idim] == FieldBoundaryType::PML) {
            do_pml = 1;
            do_pml_Hi[0][idim] = 1; // on level 0
        }
    }
    if (finest_level > 0) do_pml = 1;
    if (do_pml)
    {
#if (defined WARPX_DIM_RZ) && (defined WARPX_USE_PSATD)
        do_pml_Lo[0][0] = 0; // no PML at r=0, in cylindrical geometry
        pml_rz[0] = std::make_unique<PML_RZ>(0, boxArray(0), DistributionMap(0), &Geom(0), pml_ncell, do_pml_in_domain);
#else
        // Note: fill_guards_fields and fill_guards_current are both set to
        // zero (amrex::IntVect(0)) (what we do with damping BCs does not apply
        // to the PML, for example in the presence of mesh refinement patches)
        pml[0] = std::make_unique<PML>(0, boxArray(0), DistributionMap(0), &Geom(0), nullptr,
                             pml_ncell, pml_delta, amrex::IntVect::TheZeroVector(),
                             dt[0], nox_fft, noy_fft, noz_fft, grid_type,
                             do_moving_window, pml_has_particles, do_pml_in_domain,
                             psatd_solution_type, J_in_time, rho_in_time,
                             do_pml_dive_cleaning, do_pml_divb_cleaning,
                             amrex::IntVect(0), amrex::IntVect(0),
                             guard_cells.ng_FieldSolver.max(),
                             v_particle_pml,
                             do_pml_Lo[0], do_pml_Hi[0]);
#endif

        for (int lev = 1; lev <= finest_level; ++lev)
        {
            do_pml_Lo[lev] = amrex::IntVect::TheUnitVector();
            do_pml_Hi[lev] = amrex::IntVect::TheUnitVector();
            // check if fine patch edges co-incide with domain boundary
            amrex::Box levelBox = boxArray(lev).minimalBox();
            // Domain box at level, lev
            amrex::Box DomainBox = Geom(lev).Domain();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (levelBox.smallEnd(idim) == DomainBox.smallEnd(idim))
                    do_pml_Lo[lev][idim] = do_pml_Lo[0][idim];
                if (levelBox.bigEnd(idim) == DomainBox.bigEnd(idim))
                    do_pml_Hi[lev][idim] = do_pml_Hi[0][idim];
            }

#ifdef WARPX_DIM_RZ
            //In cylindrical geometry, if the edge of the patch is at r=0, do not add PML
            if ((max_level > 0) && (fine_tag_lo[0]==0.)) {
                do_pml_Lo[lev][0] = 0;
            }
#endif
            // Note: fill_guards_fields and fill_guards_current are both set to
            // zero (amrex::IntVect(0)) (what we do with damping BCs does not apply
            // to the PML, for example in the presence of mesh refinement patches)
            pml[lev] = std::make_unique<PML>(lev, boxArray(lev), DistributionMap(lev),
                                   &Geom(lev), &Geom(lev-1),
                                   pml_ncell, pml_delta, refRatio(lev-1),
                                   dt[lev], nox_fft, noy_fft, noz_fft, grid_type,
                                   do_moving_window, pml_has_particles, do_pml_in_domain,
                                   psatd_solution_type, J_in_time, rho_in_time, do_pml_dive_cleaning, do_pml_divb_cleaning,
                                   amrex::IntVect(0), amrex::IntVect(0),
                                   guard_cells.ng_FieldSolver.max(),
                                   v_particle_pml,
                                   do_pml_Lo[lev], do_pml_Hi[lev]);
        }
    }
}

void
WarpX::ComputePMLFactors ()
{
    if (do_pml)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            if (pml[lev])
                pml[lev]->ComputePMLFactors(dt[lev]);
        }
    }
}

void
WarpX::ComputeMaxStep ()
{
    if (do_compute_max_step_from_zmax) {
        computeMaxStepBoostAccelerator();
    }
}

/* \brief computes max_step for wakefield simulation in boosted frame.
 * \param geom: Geometry object that contains simulation domain.
 *
 * max_step is set so that the simulation stop when the lower corner of the
 * simulation box passes input parameter zmax_plasma_to_compute_max_step.
 */
void
WarpX::computeMaxStepBoostAccelerator() {
    // Sanity checks: can use zmax_plasma_to_compute_max_step only if
    // the moving window and the boost are all in z direction.
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::moving_window_dir == WARPX_ZINDEX,
        "Can use zmax_plasma_to_compute_max_step only if "
        "moving window along z. TODO: all directions.");
    if (gamma_boost > 1){
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            (WarpX::boost_direction[0]-0)*(WarpX::boost_direction[0]-0) +
            (WarpX::boost_direction[1]-0)*(WarpX::boost_direction[1]-0) +
            (WarpX::boost_direction[2]-1)*(WarpX::boost_direction[2]-1) < 1.e-12,
            "Can use zmax_plasma_to_compute_max_step in boosted frame only if "
            "warpx.boost_direction = z. TODO: all directions.");
    }

    // Lower end of the simulation domain. All quantities are given in boosted
    // frame except zmax_plasma_to_compute_max_step.

    // End of the plasma: Transform input argument
    // zmax_plasma_to_compute_max_step to boosted frame.
    const Real len_plasma_boost = zmax_plasma_to_compute_max_step/gamma_boost;
    // Plasma velocity
    const Real v_plasma_boost = -beta_boost * PhysConst::c;
    // Get time at which the lower end of the simulation domain passes the
    // upper end of the plasma (in the z direction).
    const Real interaction_time_boost = (len_plasma_boost-zmin_domain_boost_step_0)/
        (moving_window_v-v_plasma_boost);
    // Divide by dt, and update value of max_step.
    int computed_max_step;
    if (do_subcycling){
        computed_max_step = static_cast<int>(interaction_time_boost/dt[0]);
    } else {
        computed_max_step =
            static_cast<int>(interaction_time_boost/dt[maxLevel()]);
    }
    max_step = computed_max_step;
    Print()<<"max_step computed in computeMaxStepBoostAccelerator: "
           <<max_step<<std::endl;
}

void
WarpX::InitNCICorrector ()
{
#if !(defined WARPX_DIM_1D_Z)
    if (WarpX::use_fdtd_nci_corr)
    {
        for (int lev = 0; lev <= max_level; ++lev)
        {
            const Geometry& gm = Geom(lev);
            const Real* dx = gm.CellSize();
            amrex::Real dz, cdtodz;
#if defined(WARPX_DIM_3D)
                dz = dx[2];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                dz = dx[1];
#else
                dz = dx[0];
#endif
            cdtodz = PhysConst::c * dt[lev] / dz;

            // Initialize Godfrey filters
            // Same filter for fields Ex, Ey and Bz
            const bool nodal_gather = !galerkin_interpolation;
            nci_godfrey_filter_exeybz[lev] = std::make_unique<NCIGodfreyFilter>(
                godfrey_coeff_set::Ex_Ey_Bz, cdtodz, nodal_gather);
            // Same filter for fields Bx, By and Ez
            nci_godfrey_filter_bxbyez[lev] = std::make_unique<NCIGodfreyFilter>(
                godfrey_coeff_set::Bx_By_Ez, cdtodz, nodal_gather);
            // Compute Godfrey filters stencils
            nci_godfrey_filter_exeybz[lev]->ComputeStencils();
            nci_godfrey_filter_bxbyez[lev]->ComputeStencils();
        }
    }
#endif
}

void
WarpX::InitFilter (){
    if (WarpX::use_filter){
        WarpX::bilinear_filter.npass_each_dir = WarpX::filter_npass_each_dir.toArray<unsigned int>();
        WarpX::bilinear_filter.ComputeStencils();
    }
}

void
WarpX::PostRestart ()
{
    mypc->PostRestart();
}


void
WarpX::InitLevelData (int lev, Real /*time*/)
{

    ParmParse pp_warpx("warpx");

    // default values of E_external_grid and B_external_grid
    // are used to set the E and B field when "constant" or
    // "parser" is not explicitly used in the input.
    pp_warpx.query("B_ext_grid_init_style", B_ext_grid_s);
    std::transform(B_ext_grid_s.begin(),
                   B_ext_grid_s.end(),
                   B_ext_grid_s.begin(),
                   ::tolower);

    pp_warpx.query("E_ext_grid_init_style", E_ext_grid_s);
    std::transform(E_ext_grid_s.begin(),
                   E_ext_grid_s.end(),
                   E_ext_grid_s.begin(),
                   ::tolower);

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    if (B_ext_grid_s == "constant")
        utils::parser::getArrWithParser(pp_warpx, "B_external_grid", B_external_grid);

    // if the input string is "constant", the values for the
    // external grid must be provided in the input.
    if (E_ext_grid_s == "constant")
        utils::parser::getArrWithParser(pp_warpx, "E_external_grid", E_external_grid);

    // initialize the averaged fields only if the averaged algorithm
    // is activated ('psatd.do_time_averaging=1')
    ParmParse pp_psatd("psatd");
    pp_psatd.query("do_time_averaging", fft_do_time_averaging );

    for (int i = 0; i < 3; ++i) {

        if (B_ext_grid_s == "constant" || B_ext_grid_s == "default") {
           Bfield_fp[lev][i]->setVal(B_external_grid[i]);
           if (fft_do_time_averaging) {
                Bfield_avg_fp[lev][i]->setVal(B_external_grid[i]);
           }

           if (lev > 0) {
              Bfield_aux[lev][i]->setVal(B_external_grid[i]);
              Bfield_cp[lev][i]->setVal(B_external_grid[i]);
              if (fft_do_time_averaging) {
                  Bfield_avg_cp[lev][i]->setVal(B_external_grid[i]);
              }
           }
        }
        if (E_ext_grid_s == "constant" || E_ext_grid_s == "default") {
           Efield_fp[lev][i]->setVal(E_external_grid[i]);
           if (fft_do_time_averaging) {
               Efield_avg_fp[lev][i]->setVal(E_external_grid[i]);
            }

           if (lev > 0) {
              Efield_aux[lev][i]->setVal(E_external_grid[i]);
              Efield_cp[lev][i]->setVal(E_external_grid[i]);
              if (fft_do_time_averaging) {
                  Efield_avg_cp[lev][i]->setVal(E_external_grid[i]);
              }
           }
        }
    }

#ifdef AMREX_USE_EB
    InitializeEBGridData(lev);
#endif

    // if the input string for the B-field is "parse_b_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (B_ext_grid_s == "parse_b_ext_grid_function") {

#ifdef WARPX_DIM_RZ
       amrex::Abort(Utils::TextMsg::Err(
           "E and B parser for external fields does not work with RZ -- TO DO"));
#endif
       utils::parser::Store_parserString(pp_warpx, "Bx_external_grid_function(x,y,z)",
          str_Bx_ext_grid_function);
       utils::parser::Store_parserString(pp_warpx, "By_external_grid_function(x,y,z)",
          str_By_ext_grid_function);
       utils::parser::Store_parserString(pp_warpx, "Bz_external_grid_function(x,y,z)",
          str_Bz_ext_grid_function);
       Bxfield_parser = std::make_unique<amrex::Parser>(
       utils::parser::makeParser(str_Bx_ext_grid_function,{"x","y","z"}));
       Byfield_parser = std::make_unique<amrex::Parser>(
          utils::parser::makeParser(str_By_ext_grid_function,{"x","y","z"}));
       Bzfield_parser = std::make_unique<amrex::Parser>(
          utils::parser::makeParser(str_Bz_ext_grid_function,{"x","y","z"}));

       // Initialize Bfield_fp with external function
       InitializeExternalFieldsOnGridUsingParser(Bfield_fp[lev][0].get(),
                                                 Bfield_fp[lev][1].get(),
                                                 Bfield_fp[lev][2].get(),
                                                 Bxfield_parser->compile<3>(),
                                                 Byfield_parser->compile<3>(),
                                                 Bzfield_parser->compile<3>(),
                                                 m_edge_lengths[lev],
                                                 m_face_areas[lev],
                                                 'B',
                                                 lev, PatchType::fine);
       if (lev > 0) {
          InitializeExternalFieldsOnGridUsingParser(Bfield_aux[lev][0].get(),
                                                    Bfield_aux[lev][1].get(),
                                                    Bfield_aux[lev][2].get(),
                                                    Bxfield_parser->compile<3>(),
                                                    Byfield_parser->compile<3>(),
                                                    Bzfield_parser->compile<3>(),
                                                    m_edge_lengths[lev],
                                                    m_face_areas[lev],
                                                    'B',
                                                    lev, PatchType::fine);

          InitializeExternalFieldsOnGridUsingParser(Bfield_cp[lev][0].get(),
                                                    Bfield_cp[lev][1].get(),
                                                    Bfield_cp[lev][2].get(),
                                                    Bxfield_parser->compile<3>(),
                                                    Byfield_parser->compile<3>(),
                                                    Bzfield_parser->compile<3>(),
                                                    m_edge_lengths[lev],
                                                    m_face_areas[lev],
                                                    'B',
                                                    lev, PatchType::coarse);
       }
    }

    // if the input string for the E-field is "parse_e_ext_grid_function",
    // then the analytical expression or function must be
    // provided in the input file.
    if (E_ext_grid_s == "parse_e_ext_grid_function") {

#ifdef WARPX_DIM_RZ
       amrex::Abort(Utils::TextMsg::Err(
           "E and B parser for external fields does not work with RZ -- TO DO"));
#endif
       utils::parser::Store_parserString(pp_warpx, "Ex_external_grid_function(x,y,z)",
           str_Ex_ext_grid_function);
       utils::parser::Store_parserString(pp_warpx, "Ey_external_grid_function(x,y,z)",
           str_Ey_ext_grid_function);
       utils::parser::Store_parserString(pp_warpx, "Ez_external_grid_function(x,y,z)",
           str_Ez_ext_grid_function);

       Exfield_parser = std::make_unique<amrex::Parser>(
           utils::parser::makeParser(str_Ex_ext_grid_function,{"x","y","z"}));
       Eyfield_parser = std::make_unique<amrex::Parser>(
           utils::parser::makeParser(str_Ey_ext_grid_function,{"x","y","z"}));
       Ezfield_parser = std::make_unique<amrex::Parser>(
           utils::parser::makeParser(str_Ez_ext_grid_function,{"x","y","z"}));

       // Initialize Efield_fp with external function
       InitializeExternalFieldsOnGridUsingParser(Efield_fp[lev][0].get(),
                                                 Efield_fp[lev][1].get(),
                                                 Efield_fp[lev][2].get(),
                                                 Exfield_parser->compile<3>(),
                                                 Eyfield_parser->compile<3>(),
                                                 Ezfield_parser->compile<3>(),
                                                 m_edge_lengths[lev],
                                                 m_face_areas[lev],
                                                 'E',
                                                 lev, PatchType::fine);

#ifdef AMREX_USE_EB
        // We initialize ECTRhofield consistently with the Efield
        if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT) {
            m_fdtd_solver_fp[lev]->EvolveECTRho(Efield_fp[lev], m_edge_lengths[lev],
                                                m_face_areas[lev], ECTRhofield[lev], lev);

        }
#endif

       if (lev > 0) {
          InitializeExternalFieldsOnGridUsingParser(Efield_aux[lev][0].get(),
                                                    Efield_aux[lev][1].get(),
                                                    Efield_aux[lev][2].get(),
                                                    Exfield_parser->compile<3>(),
                                                    Eyfield_parser->compile<3>(),
                                                    Ezfield_parser->compile<3>(),
                                                    m_edge_lengths[lev],
                                                    m_face_areas[lev],
                                                    'E',
                                                    lev, PatchType::fine);

          InitializeExternalFieldsOnGridUsingParser(Efield_cp[lev][0].get(),
                                                    Efield_cp[lev][1].get(),
                                                    Efield_cp[lev][2].get(),
                                                    Exfield_parser->compile<3>(),
                                                    Eyfield_parser->compile<3>(),
                                                    Ezfield_parser->compile<3>(),
                                                    m_edge_lengths[lev],
                                                    m_face_areas[lev],
                                                    'E',
                                                    lev, PatchType::coarse);
#ifdef AMREX_USE_EB
           if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT) {
               // We initialize ECTRhofield consistently with the Efield
               m_fdtd_solver_cp[lev]->EvolveECTRho(Efield_cp[lev], m_edge_lengths[lev],
                                                   m_face_areas[lev], ECTRhofield[lev], lev);

           }
#endif
       }
    }

    if (costs[lev]) {
        const auto iarr = costs[lev]->IndexArray();
        for (int i : iarr) {
            (*costs[lev])[i] = 0.0;
            WarpX::setLoadBalanceEfficiency(lev, -1);
        }
    }
}

void
WarpX::InitializeExternalFieldsOnGridUsingParser (
       MultiFab *mfx, MultiFab *mfy, MultiFab *mfz,
       ParserExecutor<3> const& xfield_parser, ParserExecutor<3> const& yfield_parser,
       ParserExecutor<3> const& zfield_parser,
       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths,
       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& face_areas,
       const char field,
       const int lev, PatchType patch_type)
{

    auto dx_lev = geom[lev].CellSizeArray();
    amrex::IntVect refratio = (lev > 0 ) ? WarpX::RefRatio(lev-1) : amrex::IntVect(1);
    if (patch_type == PatchType::coarse) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            dx_lev[idim] = dx_lev[idim] * refratio[idim];
        }
    }
    const RealBox& real_box = geom[lev].ProbDomain();
    amrex::IntVect x_nodal_flag = mfx->ixType().toIntVect();
    amrex::IntVect y_nodal_flag = mfy->ixType().toIntVect();
    amrex::IntVect z_nodal_flag = mfz->ixType().toIntVect();

    for ( MFIter mfi(*mfx, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const amrex::Box& tbx = mfi.tilebox( x_nodal_flag, mfx->nGrowVect() );
       const amrex::Box& tby = mfi.tilebox( y_nodal_flag, mfy->nGrowVect() );
       const amrex::Box& tbz = mfi.tilebox( z_nodal_flag, mfz->nGrowVect() );

       auto const& mfxfab = mfx->array(mfi);
       auto const& mfyfab = mfy->array(mfi);
       auto const& mfzfab = mfz->array(mfi);

#ifdef AMREX_USE_EB
       amrex::Array4<amrex::Real> const& lx = edge_lengths[0]->array(mfi);
       amrex::Array4<amrex::Real> const& ly = edge_lengths[1]->array(mfi);
       amrex::Array4<amrex::Real> const& lz = edge_lengths[2]->array(mfi);
       amrex::Array4<amrex::Real> const& Sx = face_areas[0]->array(mfi);
       amrex::Array4<amrex::Real> const& Sy = face_areas[1]->array(mfi);
       amrex::Array4<amrex::Real> const& Sz = face_areas[2]->array(mfi);

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
       const amrex::Dim3 lx_lo = amrex::lbound(lx);
       const amrex::Dim3 lx_hi = amrex::ubound(lx);
       const amrex::Dim3 lz_lo = amrex::lbound(lz);
       const amrex::Dim3 lz_hi = amrex::ubound(lz);
#endif

#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        amrex::ignore_unused(ly, Sx, Sz);
#elif defined(WARPX_DIM_1D_Z)
        amrex::ignore_unused(lx, ly, lz, Sx, Sy, Sz);
#endif

#else
       amrex::ignore_unused(edge_lengths, face_areas, field);
#endif

        amrex::ParallelFor (tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef AMREX_USE_EB
#ifdef WARPX_DIM_3D
                if((field=='E' and lx(i, j, k)<=0) or (field=='B' and Sx(i, j, k)<=0))  return;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                //In XZ and RZ Ex is associated with a x-edge, while Bx is associated with a z-edge
                if((field=='E' and lx(i, j, k)<=0) or (field=='B' and lz(i, j, k)<=0)) return;
#endif
#endif
                // Shift required in the x-, y-, or z- position
                // depending on the index type of the multifab
#if defined(WARPX_DIM_1D_Z)
                amrex::Real x = 0._rt;
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - x_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real z = j*dx_lev[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                amrex::Real fac_x = (1._rt - x_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - x_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#else
                amrex::Real fac_x = (1._rt - x_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                amrex::Real fac_y = (1._rt - x_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                amrex::Real fac_z = (1._rt - x_nodal_flag[2]) * dx_lev[2] * 0.5_rt;
                amrex::Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the x-component of the field.
                mfxfab(i,j,k) = xfield_parser(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef AMREX_USE_EB
#ifdef WARPX_DIM_3D
                if((field=='E' and ly(i, j, k)<=0) or (field=='B' and Sy(i, j, k)<=0))  return;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                //In XZ and RZ Ey is associated with a mesh node, so we need to check if  the mesh node is covered
                if((field=='E' and (lx(std::min(i  , lx_hi.x), std::min(j  , lx_hi.y), k)<=0
                                 || lx(std::max(i-1, lx_lo.x), std::min(j  , lx_hi.y), k)<=0
                                 || lz(std::min(i  , lz_hi.x), std::min(j  , lz_hi.y), k)<=0
                                 || lz(std::min(i  , lz_hi.x), std::max(j-1, lz_lo.y), k)<=0)) or
                   (field=='B' and Sy(i,j,k)<=0)) return;
#endif
#endif
#if defined(WARPX_DIM_1D_Z)
                amrex::Real x = 0._rt;
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - y_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real z = j*dx_lev[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                amrex::Real fac_x = (1._rt - y_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - y_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif defined(WARPX_DIM_3D)
                amrex::Real fac_x = (1._rt - y_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                amrex::Real fac_y = (1._rt - y_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                amrex::Real fac_z = (1._rt - y_nodal_flag[2]) * dx_lev[2] * 0.5_rt;
                amrex::Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the y-component of the field.
                mfyfab(i,j,k)  = yfield_parser(x,y,z);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
#ifdef AMREX_USE_EB
#ifdef WARPX_DIM_3D
                if((field=='E' and lz(i, j, k)<=0) or (field=='B' and Sz(i, j, k)<=0))  return;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                //In XZ and RZ Ez is associated with a z-edge, while Bz is associated with a x-edge
                if((field=='E' and lz(i, j, k)<=0) or (field=='B' and lx(i, j, k)<=0)) return;
#endif
#endif
#if defined(WARPX_DIM_1D_Z)
                amrex::Real x = 0._rt;
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - z_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real z = j*dx_lev[0] + real_box.lo(0) + fac_z;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                amrex::Real fac_x = (1._rt - z_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                amrex::Real y = 0._rt;
                amrex::Real fac_z = (1._rt - z_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real z = j*dx_lev[1] + real_box.lo(1) + fac_z;
#elif defined(WARPX_DIM_3D)
                amrex::Real fac_x = (1._rt - z_nodal_flag[0]) * dx_lev[0] * 0.5_rt;
                amrex::Real x = i*dx_lev[0] + real_box.lo(0) + fac_x;
                amrex::Real fac_y = (1._rt - z_nodal_flag[1]) * dx_lev[1] * 0.5_rt;
                amrex::Real y = j*dx_lev[1] + real_box.lo(1) + fac_y;
                amrex::Real fac_z = (1._rt - z_nodal_flag[2]) * dx_lev[2] * 0.5_rt;
                amrex::Real z = k*dx_lev[2] + real_box.lo(2) + fac_z;
#endif
                // Initialize the z-component of the field.
                mfzfab(i,j,k) = zfield_parser(x,y,z);
            }
        );
    }
}

void
WarpX::PerformanceHints ()
{
    // Check requested MPI ranks and available boxes
    amrex::Long total_nboxes = 0; // on all MPI ranks
    for (int ilev = 0; ilev <= finestLevel(); ++ilev) {
        total_nboxes += boxArray(ilev).size();
    }
    auto const nprocs = ParallelDescriptor::NProcs();

    // Check: are there more MPI ranks than Boxes?
    if (nprocs > total_nboxes) {
        std::stringstream warnMsg;
        warnMsg << "Too many resources / too little work!\n"
            << "  It looks like you requested more compute resources than "
            << "there are total number of boxes of cells available ("
            << total_nboxes << "). "
            << "You started with (" << nprocs
            << ") MPI ranks, so (" << nprocs - total_nboxes
            << ") rank(s) will have no work.\n"
#ifdef AMREX_USE_GPU
            << "  On GPUs, consider using 1-8 boxes per GPU that together fill "
            << "each GPU's memory sufficiently. If you do not rely on dynamic "
            << "load-balancing, then one large box per GPU is ideal.\n"
#endif
            << "Consider decreasing the amr.blocking_factor and "
            << "amr.max_grid_size parameters and/or using fewer MPI ranks.\n"
            << "  More information:\n"
            << "  https://warpx.readthedocs.io/en/latest/usage/workflows/parallelization.html\n";

        ablastr::warn_manager::WMRecordWarning(
          "Performance", warnMsg.str(), ablastr::warn_manager::WarnPriority::high);
    }

#ifdef AMREX_USE_GPU
    // Check: Are there more than 12 boxes per GPU?
    if (total_nboxes > nprocs * 12) {
        std::stringstream warnMsg;
        warnMsg << "Too many boxes per GPU!\n"
            << "  It looks like you split your simulation domain "
            << "in too many boxes (" << total_nboxes << "), which "
            << "results in an average number of ("
            << amrex::Long(total_nboxes/nprocs) << ") per GPU. "
            << "This causes severe overhead in the communication of "
            << "border/guard regions.\n"
            << "  On GPUs, consider using 1-8 boxes per GPU that together fill "
            << "each GPU's memory sufficiently. If you do not rely on dynamic "
            << "load-balancing, then one large box per GPU is ideal.\n"
            << "Consider increasing the amr.blocking_factor and "
            << "amr.max_grid_size parameters and/or using more MPI ranks.\n"
            << "  More information:\n"
            << "  https://warpx.readthedocs.io/en/latest/usage/workflows/parallelization.html\n";

        ablastr::warn_manager::WMRecordWarning(
          "Performance", warnMsg.str(), ablastr::warn_manager::WarnPriority::high);
    }
#endif

    // TODO: warn if some ranks have disproportionally more work than all others
    //       tricky: it can be ok to assign "vacuum" boxes to some ranks w/o slowing down
    //               all other ranks; we need to measure this with our load-balancing
    //               routines and issue a warning only of some ranks stall all other ranks
    // TODO: check MPI-rank to GPU ratio (should be 1:1)
    // TODO: check memory per MPI rank, especially if GPUs are underutilized
    // TODO: CPU tiling hints with OpenMP
}

void WarpX::CheckGuardCells()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int dim = 0; dim < 3; ++dim)
        {
            CheckGuardCells(*Efield_fp[lev][dim]);
            CheckGuardCells(*Bfield_fp[lev][dim]);
            CheckGuardCells(*current_fp[lev][dim]);

            if (WarpX::fft_do_time_averaging)
            {
                CheckGuardCells(*Efield_avg_fp[lev][dim]);
                CheckGuardCells(*Bfield_avg_fp[lev][dim]);
            }
        }

        if (rho_fp[lev])
        {
            CheckGuardCells(*rho_fp[lev]);
        }

        if (F_fp[lev])
        {
            CheckGuardCells(*F_fp[lev]);
        }

        if (G_fp[lev])
        {
            CheckGuardCells(*G_fp[lev]);
        }

        // MultiFabs on coarse patch
        if (lev > 0)
        {
            for (int dim = 0; dim < 3; ++dim)
            {
                CheckGuardCells(*Efield_cp[lev][dim]);
                CheckGuardCells(*Bfield_cp[lev][dim]);
                CheckGuardCells(*current_cp[lev][dim]);

                if (WarpX::fft_do_time_averaging)
                {
                    CheckGuardCells(*Efield_avg_cp[lev][dim]);
                    CheckGuardCells(*Bfield_avg_cp[lev][dim]);
                }
            }

            if (rho_cp[lev])
            {
                CheckGuardCells(*rho_cp[lev]);
            }

            if (F_cp[lev])
            {
                CheckGuardCells(*F_cp[lev]);
            }

            if (G_cp[lev])
            {
                CheckGuardCells(*G_cp[lev]);
            }
        }
    }
}

void WarpX::CheckGuardCells(amrex::MultiFab const& mf)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const amrex::IntVect vc = mfi.validbox().enclosedCells().size();
        const amrex::IntVect gc = mf.nGrowVect();

        std::stringstream ss_msg;
        ss_msg << "MultiFab " << mf.tags()[1].c_str() << ":" <<
            " the number of guard cells " << gc <<
            " is larger than or equal to the number of valid cells "
            << vc << ", please reduce the number of guard cells" <<
             " or increase the grid size by changing domain decomposition.";
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(vc.allGT(gc), ss_msg.str());
    }
}

void WarpX::InitializeEBGridData (int lev)
{
#ifdef AMREX_USE_EB
    if (lev == maxLevel()) {

        // Throw a warning if EB is on and particle_shape > 1
        bool flag_eb_on = not fieldEBFactory(lev).isAllRegular();

        if ((nox > 1 or noy > 1 or noz > 1) and flag_eb_on)
        {
            ablastr::warn_manager::WMRecordWarning("Particles",
              "when algo.particle_shape > 1, numerical artifacts will be present when\n"
              "particles are close to embedded boundaries");
        }

        if (WarpX::electromagnetic_solver_id != ElectromagneticSolverAlgo::PSATD ) {

            auto const eb_fact = fieldEBFactory(lev);

            ComputeEdgeLengths(m_edge_lengths[lev], eb_fact);
            ScaleEdges(m_edge_lengths[lev], CellSize(lev));
            ComputeFaceAreas(m_face_areas[lev], eb_fact);
            ScaleAreas(m_face_areas[lev], CellSize(lev));

            if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::ECT) {
                MarkCells();
                ComputeFaceExtensions();
            }
        }

        ComputeDistanceToEB();

    }
#else
    amrex::ignore_unused(lev);
#endif
}

void WarpX::CheckKnownIssues()
{
    if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD &&
        (std::any_of(do_pml_Lo[0].begin(),do_pml_Lo[0].end(),[](const auto& ee){return ee;}) ||
        std::any_of(do_pml_Hi[0].begin(),do_pml_Hi[0].end(),[](const auto& ee){return ee;})) )
        {
            ablastr::warn_manager::WMRecordWarning(
                "PML",
                "Using PSATD together with PML may lead to instabilities if the plasma touches the PML region. "
                "It is recommended to leave enough empty space between the plasma boundary and the PML region.",
                ablastr::warn_manager::WarnPriority::low);
        }
}
