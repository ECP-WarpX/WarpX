/* Copyright 2016-2020 Andrew Myers, Ann Almgren, Aurore Blelly
 * Axel Huebl, Burlen Loring, David Grote
 * Glenn Richardson, Jean-Luc Vay, Junmin Gu
 * Mathieu Lobet, Maxence Thevenet, Michael Rowan
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao, levinem
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "FieldSolver/WarpX_FDTD.H"
#ifdef WARPX_USE_PSATD
#include "FieldSolver/SpectralSolver/SpectralKSpace.H"
#endif
#include "Python/WarpXWrappers.h"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#ifdef BL_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif

#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <utility>

using namespace amrex;

Vector<Real> WarpX::E_external_grid(3, 0.0);
Vector<Real> WarpX::B_external_grid(3, 0.0);

std::string WarpX::authors = "";
std::string WarpX::B_ext_grid_s = "default";
std::string WarpX::E_ext_grid_s = "default";

// Parser for B_external on the grid
std::string WarpX::str_Bx_ext_grid_function;
std::string WarpX::str_By_ext_grid_function;
std::string WarpX::str_Bz_ext_grid_function;
// Parser for E_external on the grid
std::string WarpX::str_Ex_ext_grid_function;
std::string WarpX::str_Ey_ext_grid_function;
std::string WarpX::str_Ez_ext_grid_function;

int WarpX::do_moving_window = 0;
int WarpX::moving_window_dir = -1;
Real WarpX::moving_window_v = std::numeric_limits<amrex::Real>::max();

bool WarpX::fft_do_time_averaging = false;

Real WarpX::quantum_xi_c2 = PhysConst::xi_c2;
Real WarpX::gamma_boost = 1._rt;
Real WarpX::beta_boost = 0._rt;
Vector<int> WarpX::boost_direction = {0,0,0};
int WarpX::do_compute_max_step_from_zmax = 0;
Real WarpX::zmax_plasma_to_compute_max_step = 0._rt;

long WarpX::current_deposition_algo;
long WarpX::charge_deposition_algo;
long WarpX::field_gathering_algo;
long WarpX::particle_pusher_algo;
int WarpX::maxwell_solver_id;
long WarpX::load_balance_costs_update_algo;
int WarpX::do_dive_cleaning = 0;
int WarpX::em_solver_medium;
int WarpX::macroscopic_solver_algo;

int WarpX::n_rz_azimuthal_modes = 1;
int WarpX::ncomps = 1;

long WarpX::nox = 1;
long WarpX::noy = 1;
long WarpX::noz = 1;

// For momentum-conserving field gathering, order of interpolation from the
// staggered positions to the grid nodes
int WarpX::field_gathering_nox = 2;
int WarpX::field_gathering_noy = 2;
int WarpX::field_gathering_noz = 2;

bool WarpX::use_fdtd_nci_corr = false;
bool WarpX::galerkin_interpolation = true;

bool WarpX::use_filter        = false;
bool WarpX::use_kspace_filter       = false;
bool WarpX::use_filter_compensation = false;
bool WarpX::use_damp_fields_in_z_guard = false;

bool WarpX::serialize_ics     = false;
bool WarpX::refine_plasma     = false;

int WarpX::num_mirrors = 0;

IntervalsParser WarpX::sort_intervals;
amrex::IntVect WarpX::sort_bin_size(AMREX_D_DECL(1,1,1));

bool WarpX::do_back_transformed_diagnostics = false;
std::string WarpX::lab_data_directory = "lab_frame_data";
int  WarpX::num_snapshots_lab = std::numeric_limits<int>::lowest();
Real WarpX::dt_snapshots_lab  = std::numeric_limits<Real>::lowest();
bool WarpX::do_back_transformed_fields = true;
bool WarpX::do_back_transformed_particles = true;

int  WarpX::num_slice_snapshots_lab = 0;
Real WarpX::dt_slice_snapshots_lab;
Real WarpX::particle_slice_width_lab = 0.0_rt;

bool WarpX::do_dynamic_scheduling = true;

int WarpX::do_electrostatic;
Real WarpX::self_fields_required_precision = 1.e-11_rt;
int WarpX::self_fields_max_iters = 200;

int WarpX::do_subcycling = 0;
bool WarpX::safe_guard_cells = 0;

IntVect WarpX::filter_npass_each_dir(1);

int WarpX::n_field_gather_buffer = -1;
int WarpX::n_current_deposition_buffer = -1;

int WarpX::do_nodal = false;

#ifdef AMREX_USE_GPU
bool WarpX::do_device_synchronize_before_profile = true;
#else
bool WarpX::do_device_synchronize_before_profile = false;
#endif

WarpX* WarpX::m_instance = nullptr;

WarpX&
WarpX::GetInstance ()
{
    if (!m_instance) {
        m_instance = new WarpX();
    }
    return *m_instance;
}

void
WarpX::ResetInstance ()
{
    delete m_instance;
    m_instance = nullptr;
}

WarpX::WarpX ()
{
    m_instance = this;
    ReadParameters();

    BackwardCompatibility();

    InitEB();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    const int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
#if 0
    // no subcycling yet
    for (int lev = 1; lev < nlevs_max; ++lev) {
        nsubsteps[lev] = MaxRefRatio(lev-1);
    }
#endif

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, std::numeric_limits<Real>::lowest());
    dt.resize(nlevs_max, std::numeric_limits<Real>::max());

    // Particle Container
    mypc = std::make_unique<MultiParticleContainer>(this);
    warpx_do_continuous_injection = mypc->doContinuousInjection();
    if (warpx_do_continuous_injection){
        if (moving_window_v >= 0){
            // Inject particles continuously from the right end of the box
            current_injection_position = geom[0].ProbHi(moving_window_dir);
        } else {
            // Inject particles continuously from the left end of the box
            current_injection_position = geom[0].ProbLo(moving_window_dir);
        }
    }
    do_back_transformed_particles = mypc->doBackTransformedDiagnostics();

    // Diagnostics
    multi_diags = std::make_unique<MultiDiagnostics>();

    /** create object for reduced diagnostics */
    reduced_diags = new MultiReducedDiags();

    Efield_aux.resize(nlevs_max);
    Bfield_aux.resize(nlevs_max);

    Efield_avg_aux.resize(nlevs_max);
    Bfield_avg_aux.resize(nlevs_max);

    F_fp.resize(nlevs_max);
    rho_fp.resize(nlevs_max);
    phi_fp.resize(nlevs_max);
    current_fp.resize(nlevs_max);
    Efield_fp.resize(nlevs_max);
    Bfield_fp.resize(nlevs_max);
    Efield_avg_fp.resize(nlevs_max);
    Bfield_avg_fp.resize(nlevs_max);

    current_store.resize(nlevs_max);

    F_cp.resize(nlevs_max);
    rho_cp.resize(nlevs_max);
    current_cp.resize(nlevs_max);
    Efield_cp.resize(nlevs_max);
    Bfield_cp.resize(nlevs_max);
    Efield_avg_cp.resize(nlevs_max);
    Bfield_avg_cp.resize(nlevs_max);

    Efield_cax.resize(nlevs_max);
    Bfield_cax.resize(nlevs_max);
    current_buffer_masks.resize(nlevs_max);
    gather_buffer_masks.resize(nlevs_max);
    current_buf.resize(nlevs_max);
    charge_buf.resize(nlevs_max);

    pml.resize(nlevs_max);
    costs.resize(nlevs_max);
    load_balance_efficiency.resize(nlevs_max);

    m_field_factory.resize(nlevs_max);

    if (em_solver_medium == MediumForEM::Macroscopic) {
        // create object for macroscopic solver
        m_macroscopic_properties = std::make_unique<MacroscopicProperties>();
    }


    // Set default values for particle and cell weights for costs update;
    // Default values listed here for the case AMREX_USE_GPU are determined
    // from single-GPU tests on Summit.
    if (costs_heuristic_cells_wt<0. && costs_heuristic_particles_wt<0.
        && WarpX::load_balance_costs_update_algo==LoadBalanceCostsUpdateAlgo::Heuristic)
    {
#ifdef AMREX_USE_GPU
        if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            switch (WarpX::nox)
            {
                case 1:
                    costs_heuristic_cells_wt = 0.575_rt;
                    costs_heuristic_particles_wt = 0.425_rt;
                    break;
                case 2:
                    costs_heuristic_cells_wt = 0.405_rt;
                    costs_heuristic_particles_wt = 0.595_rt;
                    break;
                case 3:
                    costs_heuristic_cells_wt = 0.250_rt;
                    costs_heuristic_particles_wt = 0.750_rt;
                    break;
            }
        } else { // FDTD
            switch (WarpX::nox)
            {
                case 1:
                    costs_heuristic_cells_wt = 0.401_rt;
                    costs_heuristic_particles_wt = 0.599_rt;
                    break;
                case 2:
                    costs_heuristic_cells_wt = 0.268_rt;
                    costs_heuristic_particles_wt = 0.732_rt;
                    break;
                case 3:
                    costs_heuristic_cells_wt = 0.145_rt;
                    costs_heuristic_particles_wt = 0.855_rt;
                    break;
            }
        }
#else // CPU
        costs_heuristic_cells_wt = 0.1_rt;
        costs_heuristic_particles_wt = 0.9_rt;
#endif // AMREX_USE_GPU
    }

    // Allocate field solver objects
#ifdef WARPX_USE_PSATD
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        spectral_solver_fp.resize(nlevs_max);
        spectral_solver_cp.resize(nlevs_max);
    }
#endif
    if (WarpX::maxwell_solver_id != MaxwellSolverAlgo::PSATD) {
        m_fdtd_solver_fp.resize(nlevs_max);
        m_fdtd_solver_cp.resize(nlevs_max);
    }

    // NCI Godfrey filters can have different stencils
    // at different levels (the stencil depends on c*dt/dz)
    nci_godfrey_filter_exeybz.resize(nlevs_max);
    nci_godfrey_filter_bxbyez.resize(nlevs_max);

    // Sanity checks. Must be done after calling the MultiParticleContainer
    // constructor, as it reads additional parameters
    // (e.g., use_fdtd_nci_corr)
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        AMREX_ALWAYS_ASSERT(use_fdtd_nci_corr == 0);
        AMREX_ALWAYS_ASSERT(do_subcycling == 0);
    }
}

WarpX::~WarpX ()
{
    const int nlevs_max = maxLevel() +1;
    for (int lev = 0; lev < nlevs_max; ++lev) {
        ClearLevel(lev);
    }

    delete reduced_diags;
}

void
WarpX::ReadParameters ()
{
    {
        ParmParse pp;// Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        queryWithParser(pp, "stop_time", stop_time);
        pp.query("authors", authors);
    }

    {
        ParmParse pp("amr");// Traditionally, these have prefix, amr.

        pp.query("restart", restart_chkfile);
    }

    {
        ParmParse pp("algo");
        maxwell_solver_id = GetAlgorithmInteger(pp, "maxwell_solver");
    }

    {
        ParmParse pp("warpx");

        std::vector<int> numprocs_in;
        pp.queryarr("numprocs", numprocs_in);
        if (not numprocs_in.empty()) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE
                (numprocs_in.size() == AMREX_SPACEDIM,
                 "warpx.numprocs, if specified, must have AMREX_SPACEDIM numbers");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE
                (ParallelDescriptor::NProcs() == AMREX_D_TERM(numprocs_in[0],
                                                             *numprocs_in[1],
                                                             *numprocs_in[2]),
                 "warpx.numprocs, if specified, its product must be equal to the number of processes");
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                numprocs[idim] = numprocs_in[idim];
            }
        }

        // set random seed
        std::string random_seed = "default";
        pp.query("random_seed", random_seed);
        if ( random_seed != "default" ) {
            unsigned long myproc_1 = ParallelDescriptor::MyProc() + 1;
            if ( random_seed == "random" ) {
                std::random_device rd;
                std::uniform_int_distribution<int> dist(2, INT_MAX);
                unsigned long seed = myproc_1 * dist(rd);
                ResetRandomSeed(seed);
            } else if ( std::stoi(random_seed) > 0 ) {
                unsigned long seed = myproc_1 * std::stoul(random_seed);
                ResetRandomSeed(seed);
            } else {
                Abort("warpx.random_seed must be \"default\", \"random\" or an integer > 0.");
            }
        }

        queryWithParser(pp, "cfl", cfl);
        pp.query("verbose", verbose);
        pp.query("regrid_int", regrid_int);
        pp.query("do_subcycling", do_subcycling);
        pp.query("use_hybrid_QED", use_hybrid_QED);
        pp.query("safe_guard_cells", safe_guard_cells);
        std::vector<std::string> override_sync_intervals_string_vec = {"1"};
        pp.queryarr("override_sync_intervals", override_sync_intervals_string_vec);
        override_sync_intervals = IntervalsParser(override_sync_intervals_string_vec);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(do_subcycling != 1 || max_level <= 1,
                                         "Subcycling method 1 only works for 2 levels.");

        ReadBoostedFrameParameters(gamma_boost, beta_boost, boost_direction);

        pp.query("do_device_synchronize_before_profile", do_device_synchronize_before_profile);

        // pp.query returns 1 if argument zmax_plasma_to_compute_max_step is
        // specified by the user, 0 otherwise.
        do_compute_max_step_from_zmax =
            queryWithParser(pp, "zmax_plasma_to_compute_max_step",
                      zmax_plasma_to_compute_max_step);

        pp.query("do_moving_window", do_moving_window);
        if (do_moving_window)
        {
            std::string s;
            pp.get("moving_window_dir", s);
            if (s == "x" || s == "X") {
                moving_window_dir = 0;
            }
#if (AMREX_SPACEDIM == 3)
            else if (s == "y" || s == "Y") {
                moving_window_dir = 1;
            }
#endif
            else if (s == "z" || s == "Z") {
                moving_window_dir = AMREX_SPACEDIM-1;
            }
            else {
                const std::string msg = "Unknown moving_window_dir: "+s;
                amrex::Abort(msg.c_str());
            }

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Geom(0).isPeriodic(moving_window_dir) == 0,
                       "The problem must be non-periodic in the moving window direction");

            moving_window_x = geom[0].ProbLo(moving_window_dir);

            pp.get("moving_window_v", moving_window_v);
            moving_window_v *= PhysConst::c;
        }

        pp.query("do_back_transformed_diagnostics", do_back_transformed_diagnostics);
        if (do_back_transformed_diagnostics) {

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(gamma_boost > 1.0,
                   "gamma_boost must be > 1 to use the boosted frame diagnostic.");

            pp.query("lab_data_directory", lab_data_directory);

            std::string s;
            pp.get("boost_direction", s);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (s == "z" || s == "Z"),
                   "The boosted frame diagnostic currently only works if the boost is in the z direction.");

            pp.get("num_snapshots_lab", num_snapshots_lab);

            // Read either dz_snapshots_lab or dt_snapshots_lab
            bool snapshot_interval_is_specified = 0;
            Real dz_snapshots_lab = 0;
            snapshot_interval_is_specified += queryWithParser(pp, "dt_snapshots_lab", dt_snapshots_lab);
            if ( queryWithParser(pp, "dz_snapshots_lab", dz_snapshots_lab) ){
                dt_snapshots_lab = dz_snapshots_lab/PhysConst::c;
                snapshot_interval_is_specified = 1;
            }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                snapshot_interval_is_specified,
                "When using back-transformed diagnostics, user should specify either dz_snapshots_lab or dt_snapshots_lab.");

            getWithParser(pp, "gamma_boost", gamma_boost);

            pp.query("do_back_transformed_fields", do_back_transformed_fields);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(do_moving_window,
                   "The moving window should be on if using the boosted frame diagnostic.");

            pp.get("moving_window_dir", s);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (s == "z" || s == "Z"),
                   "The boosted frame diagnostic currently only works if the moving window is in the z direction.");
        }

        do_electrostatic = GetAlgorithmInteger(pp, "do_electrostatic");

        if (do_electrostatic == ElectrostaticSolverAlgo::LabFrame) {
            queryWithParser(pp, "self_fields_required_precision", self_fields_required_precision);
            pp.query("self_fields_max_iters", self_fields_max_iters);
            // Note that with the relativistic version, these parameters would be
            // input for each species.
        }

        pp.query("n_buffer", n_buffer);
        pp.query("const_dt", const_dt);

        // Read filter and fill IntVect filter_npass_each_dir with
        // proper size for AMREX_SPACEDIM
        pp.query("use_filter", use_filter);
        pp.query("use_filter_compensation", use_filter_compensation);
        Vector<int> parse_filter_npass_each_dir(AMREX_SPACEDIM,1);
        pp.queryarr("filter_npass_each_dir", parse_filter_npass_each_dir);
        filter_npass_each_dir[0] = parse_filter_npass_each_dir[0];
        filter_npass_each_dir[1] = parse_filter_npass_each_dir[1];
#if (AMREX_SPACEDIM == 3)
        filter_npass_each_dir[2] = parse_filter_npass_each_dir[2];
#endif

#ifdef WARPX_DIM_RZ
        if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            // With RZ spectral, only use k-space filtering
            use_kspace_filter = use_filter;
            use_filter = false;
        }
#endif

        pp.query("num_mirrors", num_mirrors);
        if (num_mirrors>0){
            mirror_z.resize(num_mirrors);
            pp.getarr("mirror_z", mirror_z, 0, num_mirrors);
            mirror_z_width.resize(num_mirrors);
            pp.getarr("mirror_z_width", mirror_z_width, 0, num_mirrors);
            mirror_z_npoints.resize(num_mirrors);
            pp.getarr("mirror_z_npoints", mirror_z_npoints, 0, num_mirrors);
        }

        pp.query("serialize_ics", serialize_ics);
        pp.query("refine_plasma", refine_plasma);
        pp.query("do_dive_cleaning", do_dive_cleaning);
        pp.query("n_field_gather_buffer", n_field_gather_buffer);
        pp.query("n_current_deposition_buffer", n_current_deposition_buffer);
#ifdef AMREX_USE_GPU
        std::vector<std::string>sort_intervals_string_vec = {"4"};
#else
        std::vector<std::string> sort_intervals_string_vec = {"-1"};
#endif
        pp.queryarr("sort_intervals", sort_intervals_string_vec);
        sort_intervals = IntervalsParser(sort_intervals_string_vec);

        Vector<int> vect_sort_bin_size(AMREX_SPACEDIM,1);
        bool sort_bin_size_is_specified = pp.queryarr("sort_bin_size", vect_sort_bin_size);
        if (sort_bin_size_is_specified){
            for (int i=0; i<AMREX_SPACEDIM; i++)
                sort_bin_size[i] = vect_sort_bin_size[i];
        }

        amrex::Real quantum_xi_tmp;
        int quantum_xi_is_specified = queryWithParser(pp, "quantum_xi", quantum_xi_tmp);
        if (quantum_xi_is_specified) {
            double const quantum_xi = quantum_xi_tmp;
            quantum_xi_c2 = static_cast<amrex::Real>(quantum_xi * PhysConst::c * PhysConst::c);
        }

        pp.query("do_pml", do_pml);
        pp.query("pml_ncell", pml_ncell);
        pp.query("pml_delta", pml_delta);
        pp.query("pml_has_particles", pml_has_particles);
        pp.query("do_pml_j_damping", do_pml_j_damping);
        pp.query("do_pml_in_domain", do_pml_in_domain);
#ifdef WARPX_DIM_RZ
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( do_pml==0,
            "PML are not implemented in RZ geometry ; please set `warpx.do_pml=0`");
#endif

        Vector<int> parse_do_pml_Lo(AMREX_SPACEDIM,1);
        pp.queryarr("do_pml_Lo", parse_do_pml_Lo);
        do_pml_Lo[0] = parse_do_pml_Lo[0];
        do_pml_Lo[1] = parse_do_pml_Lo[1];
#if (AMREX_SPACEDIM == 3)
        do_pml_Lo[2] = parse_do_pml_Lo[2];
#endif
        Vector<int> parse_do_pml_Hi(AMREX_SPACEDIM,1);
        pp.queryarr("do_pml_Hi", parse_do_pml_Hi);
        do_pml_Hi[0] = parse_do_pml_Hi[0];
        do_pml_Hi[1] = parse_do_pml_Hi[1];
#if (AMREX_SPACEDIM == 3)
        do_pml_Hi[2] = parse_do_pml_Hi[2];
#endif

        if ( (do_pml_j_damping==1)&&(do_pml_in_domain==0) ){
            amrex::Abort("J-damping can only be done when PML are inside simulation domain (do_pml_in_domain=1)");
        }

        {
            // Parameters below control all plotfile diagnostics
            bool plotfile_min_max = true;
            pp.query("plotfile_min_max", plotfile_min_max);
            if (plotfile_min_max) {
                plotfile_headerversion = amrex::VisMF::Header::Version_v1;
            } else {
                plotfile_headerversion = amrex::VisMF::Header::NoFabHeader_v1;
            }
            pp.query("usesingleread", use_single_read);
            pp.query("usesinglewrite", use_single_write);
            ParmParse ppv("vismf");
            ppv.add("usesingleread", use_single_read);
            ppv.add("usesinglewrite", use_single_write);
            pp.query("mffile_nstreams", mffile_nstreams);
            VisMF::SetMFFileInStreams(mffile_nstreams);
            pp.query("field_io_nfiles", field_io_nfiles);
            VisMF::SetNOutFiles(field_io_nfiles);
            pp.query("particle_io_nfiles", particle_io_nfiles);
            ParmParse ppp("particles");
            ppp.add("particles_nfiles", particle_io_nfiles);
        }

        if (maxLevel() > 0) {
            Vector<Real> lo, hi;
            pp.getarr("fine_tag_lo", lo);
            pp.getarr("fine_tag_hi", hi);
            fine_tag_lo = RealVect{lo};
            fine_tag_hi = RealVect{hi};
        }

        pp.query("do_dynamic_scheduling", do_dynamic_scheduling);

        pp.query("do_nodal", do_nodal);
        // Use same shape factors in all directions, for gathering
        if (do_nodal) galerkin_interpolation = false;

        // Only needs to be set with WARPX_DIM_RZ, otherwise defaults to 1
        pp.query("n_rz_azimuthal_modes", n_rz_azimuthal_modes);
    }

    {
        ParmParse pp("algo");
#ifdef WARPX_DIM_RZ
        if (maxwell_solver_id == MaxwellSolverAlgo::CKC) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
                "algo.maxwell_solver = ckc is not (yet) available for RZ geometry");
        }
#endif
#ifndef WARPX_USE_PSATD
        if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
                                              "algo.maxwell_solver = psatd is not supported because WarpX was built without spectral solvers");
        }
#endif

#ifdef WARPX_DIM_RZ
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Geom(0).isPeriodic(0) == 0,
            "The problem must not be periodic in the radial direction");

        if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            // Force do_nodal=true (that is, not staggered) and
            // use same shape factors in all directions, for gathering
            do_nodal = true;
            galerkin_interpolation = false;
        }
#endif

        // note: current_deposition must be set after maxwell_solver is already determined,
        //       because its default depends on the solver selection
        current_deposition_algo = GetAlgorithmInteger(pp, "current_deposition");
        charge_deposition_algo = GetAlgorithmInteger(pp, "charge_deposition");
        particle_pusher_algo = GetAlgorithmInteger(pp, "particle_pusher");

        field_gathering_algo = GetAlgorithmInteger(pp, "field_gathering");
        if (field_gathering_algo == GatheringAlgo::MomentumConserving) {
            // Use same shape factors in all directions, for gathering
            galerkin_interpolation = false;
        }

        em_solver_medium = GetAlgorithmInteger(pp, "em_solver_medium");
        if (em_solver_medium == MediumForEM::Macroscopic ) {
            macroscopic_solver_algo = GetAlgorithmInteger(pp,"macroscopic_sigma_method");
        }

        // Load balancing parameters
        std::vector<std::string> load_balance_intervals_string_vec = {"0"};
        pp.queryarr("load_balance_intervals", load_balance_intervals_string_vec);
        load_balance_intervals = IntervalsParser(load_balance_intervals_string_vec);
        pp.query("load_balance_with_sfc", load_balance_with_sfc);
        pp.query("load_balance_knapsack_factor", load_balance_knapsack_factor);
        queryWithParser(pp, "load_balance_efficiency_ratio_threshold",
                        load_balance_efficiency_ratio_threshold);
        load_balance_costs_update_algo = GetAlgorithmInteger(pp, "load_balance_costs_update");
        queryWithParser(pp, "costs_heuristic_cells_wt", costs_heuristic_cells_wt);
        queryWithParser(pp, "costs_heuristic_particles_wt", costs_heuristic_particles_wt);
    }
    {
        ParmParse pp("interpolation");
        pp.query("nox", nox);
        pp.query("noy", noy);
        pp.query("noz", noz);

#ifdef WARPX_USE_PSATD
        if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            // For momentum-conserving field gathering, read from input the order of
            // interpolation from the staggered positions to the grid nodes
            if (field_gathering_algo == GatheringAlgo::MomentumConserving) {
                pp.query("field_gathering_nox", field_gathering_nox);
                pp.query("field_gathering_noy", field_gathering_noy);
                pp.query("field_gathering_noz", field_gathering_noz);
            }

            if (maxLevel() > 0) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    field_gathering_nox == 2 && field_gathering_noy == 2 && field_gathering_noz == 2,
                    "High-order interpolation (order > 2) is not implemented with mesh refinement");
            }

            // Host vectors for Fornberg stencil coefficients
            amrex::Vector<amrex::Real> host_Fornberg_stencil_coeffs_x;
            amrex::Vector<amrex::Real> host_Fornberg_stencil_coeffs_y;
            amrex::Vector<amrex::Real> host_Fornberg_stencil_coeffs_z;

            host_Fornberg_stencil_coeffs_x = getFornbergStencilCoefficients(field_gathering_nox, false);
            host_Fornberg_stencil_coeffs_y = getFornbergStencilCoefficients(field_gathering_noy, false);
            host_Fornberg_stencil_coeffs_z = getFornbergStencilCoefficients(field_gathering_noz, false);

            // Host vectors for ordered Fornberg stencil coefficients used for finite-order centering
            host_centering_stencil_coeffs_x.resize(field_gathering_nox);
            host_centering_stencil_coeffs_y.resize(field_gathering_noy);
            host_centering_stencil_coeffs_z.resize(field_gathering_noz);

            // Re-order Fornberg stencil coefficients
            // example for order 6: (c_0,c_1,c_2) becomes (c_2,c_1,c_0,c_0,c_1,c_2)
            ReorderFornbergCoefficients(host_centering_stencil_coeffs_x,
                                        host_Fornberg_stencil_coeffs_x, field_gathering_nox);
            ReorderFornbergCoefficients(host_centering_stencil_coeffs_y,
                                        host_Fornberg_stencil_coeffs_y, field_gathering_noy);
            ReorderFornbergCoefficients(host_centering_stencil_coeffs_z,
                                        host_Fornberg_stencil_coeffs_z, field_gathering_noz);

            // Device vectors for ordered Fornberg stencil coefficients used for finite-order centering
            device_centering_stencil_coeffs_x.resize(host_centering_stencil_coeffs_x.size());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  host_centering_stencil_coeffs_x.begin(),
                                  host_centering_stencil_coeffs_x.end(),
                                  device_centering_stencil_coeffs_x.begin());

            device_centering_stencil_coeffs_y.resize(host_centering_stencil_coeffs_y.size());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  host_centering_stencil_coeffs_y.begin(),
                                  host_centering_stencil_coeffs_y.end(),
                                  device_centering_stencil_coeffs_y.begin());

            device_centering_stencil_coeffs_z.resize(host_centering_stencil_coeffs_z.size());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  host_centering_stencil_coeffs_z.begin(),
                                  host_centering_stencil_coeffs_z.end(),
                                  device_centering_stencil_coeffs_z.begin());

            amrex::Gpu::synchronize();
        }
#endif

        pp.query("galerkin_scheme",galerkin_interpolation);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( nox == noy and nox == noz ,
            "warpx.nox, noy and noz must be equal");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( nox >= 1, "warpx.nox must >= 1");

        if (maxLevel() > 0 and nox>1 and do_pml_j_damping==1) {
            amrex::Warning("WARNING: for nox>1, some numerical artifact will be present"
                           " at coarse-fine interface. nox=1 recommended to avoid this issue");
        }
    }

    if (maxwell_solver_id == MaxwellSolverAlgo::PSATD)
    {
        ParmParse pp("psatd");
        pp.query("periodic_single_box_fft", fft_periodic_single_box);
        pp.query("fftw_plan_measure", fftw_plan_measure);

        std::string nox_str;
        std::string noy_str;
        std::string noz_str;

        pp.query("nox", nox_str);
        pp.query("noy", noy_str);
        pp.query("noz", noz_str);

        if(nox_str == "inf"){
            nox_fft = -1;
        } else{
            pp.query("nox", nox_fft);
        }
        if(noy_str == "inf"){
            noy_fft = -1;
        } else{
            pp.query("noy", noy_fft);
        }
        if(noz_str == "inf"){
            noz_fft = -1;
        } else{
            pp.query("noz", noz_fft);
        }


        if (!fft_periodic_single_box) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nox_fft > 0, "PSATD order must be finite unless psatd.periodic_single_box_fft is used");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(noy_fft > 0, "PSATD order must be finite unless psatd.periodic_single_box_fft is used");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(noz_fft > 0, "PSATD order must be finite unless psatd.periodic_single_box_fft is used");
        }

        pp.query("current_correction", current_correction);
        pp.query("v_comoving", m_v_comoving);
        pp.query("do_time_averaging", fft_do_time_averaging);

        // Check whether the default Galilean velocity should be used
        bool use_default_v_galilean = false;
        pp.query("use_default_v_galilean", use_default_v_galilean);
        if (use_default_v_galilean) {
            m_v_galilean[2] = -std::sqrt(1._rt - 1._rt / (gamma_boost * gamma_boost));
        } else {
            pp.query("v_galilean", m_v_galilean);
        }
        // Scale the Galilean/comoving velocity by the speed of light
        for (int i=0; i<3; i++) m_v_galilean[i] *= PhysConst::c;
        for (int i=0; i<3; i++) m_v_comoving[i] *= PhysConst::c;

        // The comoving PSATD algorithm is not implemented nor tested with Esirkepov current deposition
        if (current_deposition_algo == CurrentDepositionAlgo::Esirkepov) {
            if (m_v_comoving[0] != 0. || m_v_comoving[1] != 0. || m_v_comoving[2] != 0.) {
                amrex::Abort("Esirkepov current deposition cannot be used with the comoving PSATD algorithm");
            }
        }

        if (m_v_galilean[0] == 0. && m_v_galilean[1] == 0. && m_v_galilean[2] == 0. &&
            m_v_comoving[0] == 0. && m_v_comoving[1] == 0. && m_v_comoving[2] == 0.) {
            update_with_rho = false; // standard PSATD
        }
        else {
            update_with_rho = true;  // Galilean PSATD or comoving PSATD
        }

        // Overwrite update_with_rho with value set in input file
        pp.query("update_with_rho", update_with_rho);

        if (m_v_comoving[0] != 0. || m_v_comoving[1] != 0. || m_v_comoving[2] != 0.) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(update_with_rho,
                "psatd.update_with_rho must be equal to 1 for comoving PSATD");
        }

#   ifdef WARPX_DIM_RZ
        if (!Geom(0).isPeriodic(1)) {
            use_damp_fields_in_z_guard = true;
        }
        pp.query("use_damp_fields_in_z_guard", use_damp_fields_in_z_guard);
#   endif

    }

    // for slice generation //
    {
       ParmParse pp("slice");
       amrex::Vector<Real> slice_lo(AMREX_SPACEDIM);
       amrex::Vector<Real> slice_hi(AMREX_SPACEDIM);
       Vector<int> slice_crse_ratio(AMREX_SPACEDIM);
       // set default slice_crse_ratio //
       for (int idim=0; idim < AMREX_SPACEDIM; ++idim )
       {
          slice_crse_ratio[idim] = 1;
       }
       pp.queryarr("dom_lo",slice_lo,0,AMREX_SPACEDIM);
       pp.queryarr("dom_hi",slice_hi,0,AMREX_SPACEDIM);
       pp.queryarr("coarsening_ratio",slice_crse_ratio,0,AMREX_SPACEDIM);
       pp.query("plot_int",slice_plot_int);
       slice_realbox.setLo(slice_lo);
       slice_realbox.setHi(slice_hi);
       slice_cr_ratio = IntVect(AMREX_D_DECL(1,1,1));
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
       {
          if (slice_crse_ratio[idim] > 1 ) {
             slice_cr_ratio[idim] = slice_crse_ratio[idim];
          }
       }

       if (do_back_transformed_diagnostics) {
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(gamma_boost > 1.0,
                 "gamma_boost must be > 1 to use the boost frame diagnostic");
          pp.query("num_slice_snapshots_lab", num_slice_snapshots_lab);
          if (num_slice_snapshots_lab > 0) {
             getWithParser(pp, "dt_slice_snapshots_lab", dt_slice_snapshots_lab );
             getWithParser(pp, "particle_slice_width_lab",particle_slice_width_lab);
          }
       }

    }
}

void
WarpX::BackwardCompatibility ()
{
    ParmParse ppa("amr");
    int backward_int;
    if (ppa.query("plot_int", backward_int)){
        amrex::Abort("amr.plot_int is not supported anymore. Please use the new syntax for diagnostics:\n"
            "diagnostics.diags_names = my_diag\n"
            "my_diag.intervals = 10\n"
            "for output every 10 iterations. See documentation for more information");
    }

    std::string backward_str;
    if (ppa.query("plot_file", backward_str)){
        amrex::Abort("amr.plot_file is not supported anymore. "
                     "Please use the new syntax for diagnostics, see documentation.");
    }

    ParmParse ppw("warpx");
    std::vector<std::string> backward_strings;
    if (ppw.queryarr("fields_to_plot", backward_strings)){
        amrex::Abort("warpx.fields_to_plot is not supported anymore. "
                     "Please use the new syntax for diagnostics, see documentation.");
    }
    if (ppw.query("plot_finepatch", backward_int)){
        amrex::Abort("warpx.plot_finepatch is not supported anymore. "
                     "Please use the new syntax for diagnostics, see documentation.");
    }
    if (ppw.query("plot_crsepatch", backward_int)){
        amrex::Abort("warpx.plot_crsepatch is not supported anymore. "
                     "Please use the new syntax for diagnostics, see documentation.");
    }
    if (ppw.queryarr("load_balance_int", backward_strings)){
        amrex::Abort("warpx.load_balance_int is no longer a valid option. "
                     "Please use the renamed option algo.load_balance_intervals instead.");
    }
    if (ppw.queryarr("load_balance_intervals", backward_strings)){
        amrex::Abort("warpx.load_balance_intervals is no longer a valid option. "
                     "Please use the renamed option algo.load_balance_intervals instead.");
    }

    amrex::Real backward_Real;
    if (ppw.query("load_balance_efficiency_ratio_threshold", backward_Real)){
        amrex::Abort("warpx.load_balance_efficiency_ratio_threshold is not supported anymore. "
                     "Please use the renamed option algo.load_balance_efficiency_ratio_threshold.");
    }
    if (ppw.query("load_balance_with_sfc", backward_int)){
        amrex::Abort("warpx.load_balance_with_sfc is not supported anymore. "
                     "Please use the renamed option algo.load_balance_with_sfc.");
    }
    if (ppw.query("load_balance_knapsack_factor", backward_Real)){
        amrex::Abort("warpx.load_balance_knapsack_factor is not supported anymore. "
                     "Please use the renamed option algo.load_balance_knapsack_factor.");
    }
    if (ppw.queryarr("override_sync_int", backward_strings)){
        amrex::Abort("warpx.override_sync_int is no longer a valid option. "
                     "Please use the renamed option warpx.override_sync_intervals instead.");
    }
    if (ppw.queryarr("sort_int", backward_strings)){
        amrex::Abort("warpx.sort_int is no longer a valid option. "
                     "Please use the renamed option warpx.sort_intervals instead.");
    }
    if (ppw.query("use_kspace_filter", backward_int)){
        amrex::Abort("warpx.use_kspace_filter is not supported anymore. "
                     "Please use the flag use_filter, see documentation.");
    }

    ParmParse ppalgo("algo");
    int backward_mw_solver;
    if (ppalgo.query("maxwell_fdtd_solver", backward_mw_solver)){
        amrex::Abort("algo.maxwell_fdtd_solver is not supported anymore. "
                     "Please use the renamed option algo.maxwell_solver");
    }

    ParmParse pparticles("particles");
    int nspecies;
    if (pparticles.query("nspecies", nspecies)){
        amrex::Print()<<"particles.nspecies is ignored. Just use particles.species_names please.\n";
    }
    ParmParse pcol("collisions");
    int ncollisions;
    if (pcol.query("ncollisions", ncollisions)){
        amrex::Print()<<"collisions.ncollisions is ignored. Just use particles.collision_names please.\n";
    }
    ParmParse plasers("lasers");
    int nlasers;
    if (plasers.query("nlasers", nlasers)){
        amrex::Print()<<"lasers.nlasers is ignored. Just use lasers.names please.\n";
    }
}

// This is a virtual function.
void
WarpX::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                const DistributionMapping& new_dmap)
{
    AllocLevelData(lev, new_grids, new_dmap);
    InitLevelData(lev, time);
}

void
WarpX::ClearLevel (int lev)
{
    for (int i = 0; i < 3; ++i) {
        Efield_aux[lev][i].reset();
        Bfield_aux[lev][i].reset();

        current_fp[lev][i].reset();
        Efield_fp [lev][i].reset();
        Bfield_fp [lev][i].reset();

        current_store[lev][i].reset();

        current_cp[lev][i].reset();
        Efield_cp [lev][i].reset();
        Bfield_cp [lev][i].reset();

        Efield_cax[lev][i].reset();
        Bfield_cax[lev][i].reset();
        current_buf[lev][i].reset();
    }

    charge_buf[lev].reset();

    current_buffer_masks[lev].reset();
    gather_buffer_masks[lev].reset();

    F_fp  [lev].reset();
    rho_fp[lev].reset();
    phi_fp[lev].reset();
    F_cp  [lev].reset();
    rho_cp[lev].reset();

    costs[lev].reset();
    load_balance_efficiency[lev] = -1;
}

void
WarpX::AllocLevelData (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
#ifdef AMREX_USE_EB
    m_field_factory[lev] = amrex::makeEBFabFactory(Geom(lev), ba, dm,
                                             {1,1,1}, // Not clear how many ghost cells we need yet
                                             amrex::EBSupport::full);
#else
    m_field_factory[lev] = std::make_unique<FArrayBoxFactory>();
#endif

    bool aux_is_nodal = (field_gathering_algo == GatheringAlgo::MomentumConserving);

#if   (AMREX_SPACEDIM == 2)
    amrex::RealVect dx = {WarpX::CellSize(lev)[0], WarpX::CellSize(lev)[2]};
#elif (AMREX_SPACEDIM == 3)
    amrex::RealVect dx = {WarpX::CellSize(lev)[0], WarpX::CellSize(lev)[1], WarpX::CellSize(lev)[2]};
#endif

    guard_cells.Init(
        dt[lev],
        dx,
        do_subcycling,
        WarpX::use_fdtd_nci_corr,
        do_nodal,
        do_moving_window,
        aux_is_nodal,
        moving_window_dir,
        WarpX::nox,
        nox_fft, noy_fft, noz_fft,
        NCIGodfreyFilter::m_stencil_width,
        maxwell_solver_id,
        maxLevel(),
        WarpX::m_v_galilean,
        WarpX::m_v_comoving,
        safe_guard_cells,
        WarpX::do_electrostatic);

    if (mypc->nSpeciesDepositOnMainGrid() && n_current_deposition_buffer == 0) {
        n_current_deposition_buffer = 1;
        // This forces the allocation of buffers and allows the code associated
        // with buffers to run. But the buffer size of `1` is in fact not used,
        // `deposit_on_main_grid` forces all particles (whether or not they
        // are in buffers) to deposition on the main grid.
    }

    if (n_current_deposition_buffer < 0) {
        n_current_deposition_buffer = guard_cells.ng_alloc_J.max();
    }
    if (n_field_gather_buffer < 0) {
        // Field gather buffer should be larger than current deposition buffers
        n_field_gather_buffer = n_current_deposition_buffer + 1;
    }

    AllocLevelMFs(lev, ba, dm, guard_cells.ng_alloc_EB, guard_cells.ng_alloc_J,
                  guard_cells.ng_alloc_Rho, guard_cells.ng_alloc_F,
                  guard_cells.ng_Extra, aux_is_nodal);
}

void
WarpX::AllocLevelMFs (int lev, const BoxArray& ba, const DistributionMapping& dm,
                      const IntVect& ngE, const IntVect& ngJ, const IntVect& ngRho,
                      const IntVect& ngF, const IntVect& ngextra, const bool aux_is_nodal)
{
    // Declare nodal flags
    IntVect Ex_nodal_flag, Ey_nodal_flag, Ez_nodal_flag;
    IntVect Bx_nodal_flag, By_nodal_flag, Bz_nodal_flag;
    IntVect jx_nodal_flag, jy_nodal_flag, jz_nodal_flag;
    IntVect rho_nodal_flag;
    IntVect phi_nodal_flag;

    // Set nodal flags
#if   (AMREX_SPACEDIM == 2)
    // AMReX convention: x = first dimension, y = missing dimension, z = second dimension
    Ex_nodal_flag = IntVect(0,1);
    Ey_nodal_flag = IntVect(1,1);
    Ez_nodal_flag = IntVect(1,0);
    Bx_nodal_flag = IntVect(1,0);
    By_nodal_flag = IntVect(0,0);
    Bz_nodal_flag = IntVect(0,1);
    jx_nodal_flag = IntVect(0,1);
    jy_nodal_flag = IntVect(1,1);
    jz_nodal_flag = IntVect(1,0);
#elif (AMREX_SPACEDIM == 3)
    Ex_nodal_flag = IntVect(0,1,1);
    Ey_nodal_flag = IntVect(1,0,1);
    Ez_nodal_flag = IntVect(1,1,0);
    Bx_nodal_flag = IntVect(1,0,0);
    By_nodal_flag = IntVect(0,1,0);
    Bz_nodal_flag = IntVect(0,0,1);
    jx_nodal_flag = IntVect(0,1,1);
    jy_nodal_flag = IntVect(1,0,1);
    jz_nodal_flag = IntVect(1,1,0);
#endif
    rho_nodal_flag = IntVect( AMREX_D_DECL(1,1,1) );
    phi_nodal_flag = IntVect::TheNodeVector();

    // Overwrite nodal flags if necessary
    if (do_nodal) {
        Ex_nodal_flag  = IntVect::TheNodeVector();
        Ey_nodal_flag  = IntVect::TheNodeVector();
        Ez_nodal_flag  = IntVect::TheNodeVector();
        Bx_nodal_flag  = IntVect::TheNodeVector();
        By_nodal_flag  = IntVect::TheNodeVector();
        Bz_nodal_flag  = IntVect::TheNodeVector();
        jx_nodal_flag  = IntVect::TheNodeVector();
        jy_nodal_flag  = IntVect::TheNodeVector();
        jz_nodal_flag  = IntVect::TheNodeVector();
        rho_nodal_flag = IntVect::TheNodeVector();
    }
#ifdef WARPX_DIM_RZ
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        // Force cell-centered IndexType in r and z
        Ex_nodal_flag  = IntVect::TheCellVector();
        Ey_nodal_flag  = IntVect::TheCellVector();
        Ez_nodal_flag  = IntVect::TheCellVector();
        Bx_nodal_flag  = IntVect::TheCellVector();
        By_nodal_flag  = IntVect::TheCellVector();
        Bz_nodal_flag  = IntVect::TheCellVector();
        jx_nodal_flag  = IntVect::TheCellVector();
        jy_nodal_flag  = IntVect::TheCellVector();
        jz_nodal_flag  = IntVect::TheCellVector();
        rho_nodal_flag = IntVect::TheCellVector();
    }

    // With RZ multimode, there is a real and imaginary component
    // for each mode, except mode 0 which is purely real
    // Component 0 is mode 0.
    // Odd components are the real parts.
    // Even components are the imaginary parts.
    ncomps = n_rz_azimuthal_modes*2 - 1;
#endif

    // set human-readable tag for each MultiFab
    auto const tag = [lev]( std::string tagname ) {
        tagname.append("[l=").append(std::to_string(lev)).append("]");
        return MFInfo().SetTag(std::move(tagname));
    };

    //
    // The fine patch
    //
    std::array<Real,3> dx = CellSize(lev);

    Bfield_fp[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Bx_nodal_flag),dm,ncomps,ngE+ngextra,tag("Bfield_fp[x]"));
    Bfield_fp[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,By_nodal_flag),dm,ncomps,ngE+ngextra,tag("Bfield_fp[y]"));
    Bfield_fp[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Bz_nodal_flag),dm,ncomps,ngE+ngextra,tag("Bfield_fp[z]"));

    Efield_fp[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Ex_nodal_flag),dm,ncomps,ngE+ngextra,tag("Efield_fp[x]"));
    Efield_fp[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,Ey_nodal_flag),dm,ncomps,ngE+ngextra,tag("Efield_fp[y]"));
    Efield_fp[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Ez_nodal_flag),dm,ncomps,ngE+ngextra,tag("Efield_fp[z]"));

    current_fp[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,jx_nodal_flag),dm,ncomps,ngJ,tag("current_fp[x]"));
    current_fp[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,jy_nodal_flag),dm,ncomps,ngJ,tag("current_fp[y]"));
    current_fp[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,jz_nodal_flag),dm,ncomps,ngJ,tag("current_fp[z]"));


    Bfield_avg_fp[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Bx_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_fp[x]"));
    Bfield_avg_fp[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,By_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_fp[y]"));
    Bfield_avg_fp[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Bz_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_fp[z]"));

    Efield_avg_fp[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Ex_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_fp[x]"));
    Efield_avg_fp[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,Ey_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_fp[y]"));
    Efield_avg_fp[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Ez_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_fp[z]"));

    bool deposit_charge = do_dive_cleaning || (plot_rho && do_back_transformed_diagnostics);
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
        deposit_charge = do_dive_cleaning || (plot_rho && do_back_transformed_diagnostics)
                         || update_with_rho || current_correction;
    }
    if (deposit_charge)
    {
        rho_fp[lev] = std::make_unique<MultiFab>(amrex::convert(ba,rho_nodal_flag),dm,2*ncomps,ngRho,tag("rho_fp"));
    }

    if (do_electrostatic == ElectrostaticSolverAlgo::LabFrame)
    {
        IntVect ngPhi = IntVect( AMREX_D_DECL(1,1,1) );
        phi_fp[lev] = std::make_unique<MultiFab>(amrex::convert(ba,phi_nodal_flag),dm,ncomps,ngPhi,tag("phi_fp"));
    }

    if (do_subcycling == 1 && lev == 0)
    {
        current_store[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,jx_nodal_flag),dm,ncomps,ngJ,tag("current_store[x]"));
        current_store[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,jy_nodal_flag),dm,ncomps,ngJ,tag("current_store[y]"));
        current_store[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,jz_nodal_flag),dm,ncomps,ngJ,tag("current_store[z]"));
    }

    if (do_dive_cleaning)
    {
        F_fp[lev] = std::make_unique<MultiFab>(amrex::convert(ba,IntVect::TheUnitVector()),dm,ncomps, ngF.max(),tag("F_fp"));
    }
    if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD)
    {
        // Allocate and initialize the spectral solver
#ifndef WARPX_USE_PSATD
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
            "WarpX::AllocLevelMFs: PSATD solver requires WarpX build with spectral solver support.");
#else
        if (!do_dive_cleaning)
            rho_fp[lev] = std::make_unique<MultiFab>(amrex::convert(ba,rho_nodal_flag),dm,2*ncomps,ngRho,tag("rho_fp"));

#   if (AMREX_SPACEDIM == 3)
        RealVect dx_vect(dx[0], dx[1], dx[2]);
#   elif (AMREX_SPACEDIM == 2)
        RealVect dx_vect(dx[0], dx[2]);
#   endif
        // Check whether the option periodic, single box is valid here
        if (fft_periodic_single_box) {
#   ifdef WARPX_DIM_RZ
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                geom[0].isPeriodic(1)          // domain is periodic in z
                && ba.size() == 1 && lev == 0, // domain is decomposed in a single box
                "The option `psatd.periodic_single_box_fft` can only be used for a periodic domain, decomposed in a single box");
#   else
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                geom[0].isAllPeriodic()        // domain is periodic in all directions
                && ba.size() == 1 && lev == 0, // domain is decomposed in a single box
                "The option `psatd.periodic_single_box_fft` can only be used for a periodic domain, decomposed in a single box");
#   endif
        }
        // Get the cell-centered box
        BoxArray realspace_ba = ba;  // Copy box
        realspace_ba.enclosedCells(); // Make it cell-centered
        // Define spectral solver
#   ifdef WARPX_DIM_RZ
        if ( fft_periodic_single_box == false ) {
            realspace_ba.grow(1, ngE[1]); // add guard cells only in z
        }
        spectral_solver_fp[lev] = std::make_unique<SpectralSolverRZ>( realspace_ba, dm,
            n_rz_azimuthal_modes, noz_fft, do_nodal, m_v_galilean, dx_vect, dt[lev], lev, update_with_rho );
        if (use_kspace_filter) {
            spectral_solver_fp[lev]->InitFilter(filter_npass_each_dir, use_filter_compensation);
        }
#   else
        if ( fft_periodic_single_box == false ) {
            realspace_ba.grow(ngE); // add guard cells
        }
        bool const pml_flag_false = false;
        spectral_solver_fp[lev] = std::make_unique<SpectralSolver>( realspace_ba, dm,
            nox_fft, noy_fft, noz_fft, do_nodal, m_v_galilean, m_v_comoving, dx_vect, dt[lev],
            pml_flag_false, fft_periodic_single_box, update_with_rho, fft_do_time_averaging );
#   endif
#endif
    } // MaxwellSolverAlgo::PSATD
    else {
        m_fdtd_solver_fp[lev] = std::make_unique<FiniteDifferenceSolver>(maxwell_solver_id, dx, do_nodal);
    }

    //
    // The Aux patch (i.e., the full solution)
    //
    if (aux_is_nodal and !do_nodal)
    {
        // Create aux multifabs on Nodal Box Array
        BoxArray const nba = amrex::convert(ba,IntVect::TheNodeVector());
        Bfield_aux[lev][0] = std::make_unique<MultiFab>(nba,dm,ncomps,ngE,tag("Bfield_aux[x]"));
        Bfield_aux[lev][1] = std::make_unique<MultiFab>(nba,dm,ncomps,ngE,tag("Bfield_aux[y]"));
        Bfield_aux[lev][2] = std::make_unique<MultiFab>(nba,dm,ncomps,ngE,tag("Bfield_aux[z]"));

        Efield_aux[lev][0] = std::make_unique<MultiFab>(nba,dm,ncomps,ngE,tag("Efield_aux[x]"));
        Efield_aux[lev][1] = std::make_unique<MultiFab>(nba,dm,ncomps,ngE,tag("Efield_aux[y]"));
        Efield_aux[lev][2] = std::make_unique<MultiFab>(nba,dm,ncomps,ngE,tag("Efield_aux[z]"));
    }
    else if (lev == 0)
    {
        for (int idir = 0; idir < 3; ++idir) {
            Efield_aux[lev][idir] = std::make_unique<MultiFab>(*Efield_fp[lev][idir], amrex::make_alias, 0, ncomps);
            Bfield_aux[lev][idir] = std::make_unique<MultiFab>(*Bfield_fp[lev][idir], amrex::make_alias, 0, ncomps);

            Efield_avg_aux[lev][idir] = std::make_unique<MultiFab>(*Efield_avg_fp[lev][idir], amrex::make_alias, 0, ncomps);
            Bfield_avg_aux[lev][idir] = std::make_unique<MultiFab>(*Bfield_avg_fp[lev][idir], amrex::make_alias, 0, ncomps);
        }
    }
    else
    {
        Bfield_aux[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Bx_nodal_flag),dm,ncomps,ngE,tag("Bfield_aux[x]"));
        Bfield_aux[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,By_nodal_flag),dm,ncomps,ngE,tag("Bfield_aux[y]"));
        Bfield_aux[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Bz_nodal_flag),dm,ncomps,ngE,tag("Bfield_aux[z]"));

        Efield_aux[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Ex_nodal_flag),dm,ncomps,ngE,tag("Efield_aux[x]"));
        Efield_aux[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,Ey_nodal_flag),dm,ncomps,ngE,tag("Efield_aux[y]"));
        Efield_aux[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Ez_nodal_flag),dm,ncomps,ngE,tag("Efield_aux[z]"));


        Bfield_avg_aux[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Bx_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_aux[x]"));
        Bfield_avg_aux[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,By_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_aux[y]"));
        Bfield_avg_aux[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Bz_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_aux[z]"));

        Efield_avg_aux[lev][0] = std::make_unique<MultiFab>(amrex::convert(ba,Ex_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_aux[x]"));
        Efield_avg_aux[lev][1] = std::make_unique<MultiFab>(amrex::convert(ba,Ey_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_aux[y]"));
        Efield_avg_aux[lev][2] = std::make_unique<MultiFab>(amrex::convert(ba,Ez_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_aux[z]"));

    }

    //
    // The coarse patch
    //
    if (lev > 0)
    {
        BoxArray cba = ba;
        cba.coarsen(refRatio(lev-1));
        std::array<Real,3> cdx = CellSize(lev-1);

        // Create the MultiFabs for B
        Bfield_cp[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,Bx_nodal_flag),dm,ncomps,ngE,tag("Bfield_cp[x]"));
        Bfield_cp[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,By_nodal_flag),dm,ncomps,ngE,tag("Bfield_cp[y]"));
        Bfield_cp[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,Bz_nodal_flag),dm,ncomps,ngE,tag("Bfield_cp[z]"));

        // Create the MultiFabs for E
        Efield_cp[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,Ex_nodal_flag),dm,ncomps,ngE,tag("Efield_cp[x]"));
        Efield_cp[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,Ey_nodal_flag),dm,ncomps,ngE,tag("Efield_cp[y]"));
        Efield_cp[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,Ez_nodal_flag),dm,ncomps,ngE,tag("Efield_cp[z]"));

        // Create the MultiFabs for B_avg
        Bfield_avg_cp[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,Bx_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_cp[x]"));
        Bfield_avg_cp[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,By_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_cp[y]"));
        Bfield_avg_cp[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,Bz_nodal_flag),dm,ncomps,ngE,tag("Bfield_avg_cp[z]"));

        // Create the MultiFabs for E_avg
        Efield_avg_cp[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,Ex_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_cp[x]"));
        Efield_avg_cp[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,Ey_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_cp[y]"));
        Efield_avg_cp[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,Ez_nodal_flag),dm,ncomps,ngE,tag("Efield_avg_cp[z]"));

        // Create the MultiFabs for the current
        current_cp[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,jx_nodal_flag),dm,ncomps,ngJ,tag("current_cp[x]"));
        current_cp[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,jy_nodal_flag),dm,ncomps,ngJ,tag("current_cp[y]"));
        current_cp[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,jz_nodal_flag),dm,ncomps,ngJ,tag("current_cp[z]"));

        if (do_dive_cleaning || (plot_rho && do_back_transformed_diagnostics)) {
            rho_cp[lev] = std::make_unique<MultiFab>(amrex::convert(cba,rho_nodal_flag),dm,2*ncomps,ngRho,tag("rho_cp"));
        }
        if (do_dive_cleaning)
        {
            F_cp[lev] = std::make_unique<MultiFab>(amrex::convert(cba,IntVect::TheUnitVector()),dm,ncomps, ngF.max(),tag("F_cp"));
        }
        if (WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD)
        {
            // Allocate and initialize the spectral solver
#ifndef WARPX_USE_PSATD
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false,
                "WarpX::AllocLevelMFs: PSATD solver requires WarpX build with spectral solver support.");
#else
            if (!do_dive_cleaning)
                rho_cp[lev] = std::make_unique<MultiFab>( amrex::convert(cba,rho_nodal_flag),dm,2*ncomps,ngRho,tag("rho_cp") );

#   if (AMREX_SPACEDIM == 3)
            RealVect cdx_vect(cdx[0], cdx[1], cdx[2]);
#   elif (AMREX_SPACEDIM == 2)
            RealVect cdx_vect(cdx[0], cdx[2]);
#   endif
            // Get the cell-centered box, with guard cells
            BoxArray c_realspace_ba = cba;// Copy box
            c_realspace_ba.enclosedCells(); // Make it cell-centered
            // Define spectral solver
#   ifdef WARPX_DIM_RZ
            c_realspace_ba.grow(1, ngE[1]); // add guard cells only in z
            spectral_solver_cp[lev] = std::make_unique<SpectralSolverRZ>( c_realspace_ba, dm,
                n_rz_azimuthal_modes, noz_fft, do_nodal, m_v_galilean, cdx_vect, dt[lev], lev, update_with_rho );
            if (use_kspace_filter) {
                spectral_solver_cp[lev]->InitFilter(filter_npass_each_dir, use_filter_compensation);
            }
#   else
            c_realspace_ba.grow(ngE); // add guard cells
            bool const pml_flag_false = false;
            spectral_solver_cp[lev] = std::make_unique<SpectralSolver>( c_realspace_ba, dm,
                nox_fft, noy_fft, noz_fft, do_nodal, m_v_galilean, m_v_comoving, cdx_vect, dt[lev],
                pml_flag_false, fft_periodic_single_box, update_with_rho, fft_do_time_averaging );
#   endif
#endif
        } // MaxwellSolverAlgo::PSATD
        else {
            m_fdtd_solver_cp[lev] = std::make_unique<FiniteDifferenceSolver>(maxwell_solver_id, cdx,
                                                                             do_nodal);
        }
    }

    //
    // Copy of the coarse aux
    //
    if (lev > 0 && (n_field_gather_buffer > 0 || n_current_deposition_buffer > 0 ||
                    mypc->nSpeciesGatherFromMainGrid() > 0))
    {
        BoxArray cba = ba;
        cba.coarsen(refRatio(lev-1));

        if (n_field_gather_buffer > 0 || mypc->nSpeciesGatherFromMainGrid() > 0) {
            if (aux_is_nodal) {
                BoxArray const& cnba = amrex::convert(cba,IntVect::TheNodeVector());
                Bfield_cax[lev][0] = std::make_unique<MultiFab>(cnba,dm,ncomps,ngE,tag("Bfield_cax[x]"));
                Bfield_cax[lev][1] = std::make_unique<MultiFab>(cnba,dm,ncomps,ngE,tag("Bfield_cax[y]"));
                Bfield_cax[lev][2] = std::make_unique<MultiFab>(cnba,dm,ncomps,ngE,tag("Bfield_cax[z]"));
                Efield_cax[lev][0] = std::make_unique<MultiFab>(cnba,dm,ncomps,ngE,tag("Efield_cax[x]"));
                Efield_cax[lev][1] = std::make_unique<MultiFab>(cnba,dm,ncomps,ngE,tag("Efield_cax[y]"));
                Efield_cax[lev][2] = std::make_unique<MultiFab>(cnba,dm,ncomps,ngE,tag("Efield_cax[z]"));
            } else {
                // Create the MultiFabs for B
                Bfield_cax[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,Bx_nodal_flag),dm,ncomps,ngE,tag("Bfield_cax[x]"));
                Bfield_cax[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,By_nodal_flag),dm,ncomps,ngE,tag("Bfield_cax[y]"));
                Bfield_cax[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,Bz_nodal_flag),dm,ncomps,ngE,tag("Bfield_cax[z]"));

                // Create the MultiFabs for E
                Efield_cax[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,Ex_nodal_flag),dm,ncomps,ngE,tag("Efield_cax[x]"));
                Efield_cax[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,Ey_nodal_flag),dm,ncomps,ngE,tag("Efield_cax[y]"));
                Efield_cax[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,Ez_nodal_flag),dm,ncomps,ngE,tag("Efield_cax[z]"));
            }

            gather_buffer_masks[lev] = std::make_unique<iMultiFab>(ba, dm, ncomps, 1 );
            // Gather buffer masks have 1 ghost cell, because of the fact
            // that particles may move by more than one cell when using subcycling.
        }

        if (n_current_deposition_buffer > 0) {
            current_buf[lev][0] = std::make_unique<MultiFab>(amrex::convert(cba,jx_nodal_flag),dm,ncomps,ngJ,tag("current_buf[x]"));
            current_buf[lev][1] = std::make_unique<MultiFab>(amrex::convert(cba,jy_nodal_flag),dm,ncomps,ngJ,tag("current_buf[y]"));
            current_buf[lev][2] = std::make_unique<MultiFab>(amrex::convert(cba,jz_nodal_flag),dm,ncomps,ngJ,tag("current_buf[z]"));
            if (rho_cp[lev]) {
                charge_buf[lev] = std::make_unique<MultiFab>(amrex::convert(cba,rho_nodal_flag),dm,2*ncomps,ngRho,tag("charge_buf"));
            }
            current_buffer_masks[lev] = std::make_unique<iMultiFab>(ba, dm, ncomps, 1);
            // Current buffer masks have 1 ghost cell, because of the fact
            // that particles may move by more than one cell when using subcycling.
        }
    }

    if (load_balance_intervals.isActivated())
    {
        costs[lev] = std::make_unique<LayoutData<Real>>(ba, dm);
        load_balance_efficiency[lev] = -1;
    }
}

std::array<Real,3>
WarpX::CellSize (int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const Real* dx = gm.CellSize();
#if (AMREX_SPACEDIM == 3)
    return { dx[0], dx[1], dx[2] };
#elif (AMREX_SPACEDIM == 2)
    return { dx[0], 1.0, dx[1] };
#else
    static_assert(AMREX_SPACEDIM != 1, "1D is not supported");
#endif
}

amrex::RealBox
WarpX::getRealBox(const Box& bx, int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const RealBox grid_box{bx, gm.CellSize(), gm.ProbLo()};
    return( grid_box );
}

std::array<Real,3>
WarpX::LowerCorner(const Box& bx, std::array<amrex::Real,3> galilean_shift, int lev)
{
    RealBox grid_box = getRealBox( bx, lev );

    const Real* xyzmin = grid_box.lo();

#if (AMREX_SPACEDIM == 3)
    return { xyzmin[0] + galilean_shift[0], xyzmin[1] + galilean_shift[1], xyzmin[2] + galilean_shift[2] };

#elif (AMREX_SPACEDIM == 2)
    return { xyzmin[0] + galilean_shift[0], std::numeric_limits<Real>::lowest(), xyzmin[1] + galilean_shift[2] };
#endif
}

std::array<Real,3>
WarpX::UpperCorner(const Box& bx, int lev)
{
    const RealBox grid_box = getRealBox( bx, lev );
    const Real* xyzmax = grid_box.hi();
#if (AMREX_SPACEDIM == 3)
    return { xyzmax[0], xyzmax[1], xyzmax[2] };
#elif (AMREX_SPACEDIM == 2)
    return { xyzmax[0], std::numeric_limits<Real>::max(), xyzmax[1] };
#endif
}

std::array<Real,3>
WarpX::LowerCornerWithGalilean (const Box& bx, const amrex::Array<amrex::Real,3>& v_galilean, int lev)
{
    amrex::Real cur_time = gett_new(lev);
    amrex::Real time_shift = (cur_time - time_of_last_gal_shift);
    amrex::Array<amrex::Real,3> galilean_shift = { v_galilean[0]*time_shift, v_galilean[1]*time_shift, v_galilean[2]*time_shift };
    return WarpX::LowerCorner(bx, galilean_shift, lev);
}

IntVect
WarpX::RefRatio (int lev)
{
    return GetInstance().refRatio(lev);
}

void
WarpX::ComputeDivB (amrex::MultiFab& divB, int const dcomp,
                    const std::array<const amrex::MultiFab* const, 3>& B,
                    const std::array<amrex::Real,3>& dx)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!do_nodal,
        "ComputeDivB not implemented with do_nodal."
        "Shouldn't be too hard to make it general with class FiniteDifferenceSolver");

    Real dxinv = 1._rt/dx[0], dyinv = 1._rt/dx[1], dzinv = 1._rt/dx[2];

#ifdef WARPX_DIM_RZ
    const Real rmin = GetInstance().Geom(0).ProbLo(0);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(divB, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& Bxfab = B[0]->array(mfi);
        auto const& Byfab = B[1]->array(mfi);
        auto const& Bzfab = B[2]->array(mfi);
        auto const& divBfab = divB.array(mfi);

        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            warpx_computedivb(i, j, k, dcomp, divBfab, Bxfab, Byfab, Bzfab, dxinv, dyinv, dzinv
#ifdef WARPX_DIM_RZ
                              ,rmin
#endif
                              );
        });
    }
}

void
WarpX::ComputeDivB (amrex::MultiFab& divB, int const dcomp,
                    const std::array<const amrex::MultiFab* const, 3>& B,
                    const std::array<amrex::Real,3>& dx, int const ngrow)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!do_nodal,
        "ComputeDivB not implemented with do_nodal."
        "Shouldn't be too hard to make it general with class FiniteDifferenceSolver");

    Real dxinv = 1._rt/dx[0], dyinv = 1._rt/dx[1], dzinv = 1._rt/dx[2];

#ifdef WARPX_DIM_RZ
    const Real rmin = GetInstance().Geom(0).ProbLo(0);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(divB, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(ngrow);
        auto const& Bxfab = B[0]->array(mfi);
        auto const& Byfab = B[1]->array(mfi);
        auto const& Bzfab = B[2]->array(mfi);
        auto const& divBfab = divB.array(mfi);

        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            warpx_computedivb(i, j, k, dcomp, divBfab, Bxfab, Byfab, Bzfab, dxinv, dyinv, dzinv
#ifdef WARPX_DIM_RZ
                              ,rmin
#endif
                              );
        });
    }
}

void
WarpX::ComputeDivE(amrex::MultiFab& divE, const int lev)
{
    if ( WarpX::maxwell_solver_id == MaxwellSolverAlgo::PSATD ) {
#ifdef WARPX_USE_PSATD
        spectral_solver_fp[lev]->ComputeSpectralDivE( Efield_aux[lev], divE );
#else
        amrex::Abort("ComputeDivE: PSATD requested but not compiled");
#endif
    } else {
        m_fdtd_solver_fp[lev]->ComputeDivE( Efield_aux[lev], divE );
    }
}

PML*
WarpX::GetPML (int lev)
{
    if (do_pml) {
        // This should check if pml was initialized
        return pml[lev].get();
    } else {
        return nullptr;
    }
}

std::vector< bool >
WarpX::getPMLdirections() const
{
    std::vector< bool > dirsWithPML( 6, false );
#if AMREX_SPACEDIM!=3
    dirsWithPML.resize( 4 );
#endif
    if( do_pml )
    {
        for( auto i = 0u; i < dirsWithPML.size() / 2u; ++i )
        {
            dirsWithPML.at( 2u*i      ) = bool(do_pml_Lo[i]);
            dirsWithPML.at( 2u*i + 1u ) = bool(do_pml_Hi[i]);
        }
    }
    return dirsWithPML;
}

amrex::LayoutData<amrex::Real>*
WarpX::getCosts (int lev)
{
    if (m_instance)
    {
        return m_instance->costs[lev].get();
    } else
    {
        return nullptr;
    }
}

void
WarpX::BuildBufferMasks ()
{
    for (int lev = 1; lev <= maxLevel(); ++lev)
    {
        for (int ipass = 0; ipass < 2; ++ipass)
        {
            int ngbuffer = (ipass == 0) ? n_current_deposition_buffer : n_field_gather_buffer;
            iMultiFab* bmasks = (ipass == 0) ? current_buffer_masks[lev].get() : gather_buffer_masks[lev].get();
            if (bmasks)
            {
                const int ngtmp = ngbuffer + bmasks->nGrow();
                iMultiFab tmp(bmasks->boxArray(), bmasks->DistributionMap(), 1, ngtmp);
                const int covered = 1;
                const int notcovered = 0;
                const int physbnd = 1;
                const int interior = 1;
                const Box& dom = Geom(lev).Domain();
                const auto& period = Geom(lev).periodicity();
                tmp.BuildMask(dom, period, covered, notcovered, physbnd, interior);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*bmasks, true); mfi.isValid(); ++mfi)
                {
                    const Box tbx = mfi.growntilebox();
                    BuildBufferMasksInBox( tbx, (*bmasks)[mfi], tmp[mfi], ngbuffer );
                }
            }
        }
    }
}

/**
 * \brief Build buffer mask within given FArrayBox
 *
 * \param tbx         Current FArrayBox
 * \param buffer_mask Buffer mask to be set
 * \param guard_mask  Guard mask used to set buffer_mask
 * \param ng          Number of guard cells
 */
void
WarpX::BuildBufferMasksInBox ( const amrex::Box tbx, amrex::IArrayBox &buffer_mask,
                               const amrex::IArrayBox &guard_mask, const int ng )
{
    bool setnull;
    const auto lo = amrex::lbound( tbx );
    const auto hi = amrex::ubound( tbx );
    Array4<int> msk = buffer_mask.array();
    Array4<int const> gmsk = guard_mask.array();
#if (AMREX_SPACEDIM == 2)
    int k = lo.z;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            setnull = false;
            // If gmsk=0 for any neighbor within ng cells, current cell is in the buffer region
            for     (int jj = j-ng; jj <= j+ng; ++jj) {
                for (int ii = i-ng; ii <= i+ng; ++ii) {
                    if ( gmsk(ii,jj,k) == 0 ) setnull = true;
                }
            }
            if ( setnull ) msk(i,j,k) = 0;
            else           msk(i,j,k) = 1;
        }
    }
#elif (AMREX_SPACEDIM == 3)
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                setnull = false;
                // If gmsk=0 for any neighbor within ng cells, current cell is in the buffer region
                for         (int kk = k-ng; kk <= k+ng; ++kk) {
                    for     (int jj = j-ng; jj <= j+ng; ++jj) {
                        for (int ii = i-ng; ii <= i+ng; ++ii) {
                            if ( gmsk(ii,jj,kk) == 0 ) setnull = true;
                        }
                    }
                }
                if ( setnull ) msk(i,j,k) = 0;
                else           msk(i,j,k) = 1;
            }
        }
    }
#endif
}

void WarpX::ReorderFornbergCoefficients(amrex::Vector<amrex::Real>& ordered_coeffs,
                                        amrex::Vector<amrex::Real>& unordered_coeffs,
                                        const int order)
{
    const int n = order / 2;
    for (int i = 0; i < n; i++) {
        ordered_coeffs[i] = unordered_coeffs[n-1-i];
    }
    for (int i = n; i < order; i++) {
        ordered_coeffs[i] = unordered_coeffs[i-n];
    }
}

const iMultiFab*
WarpX::CurrentBufferMasks (int lev)
{
    return GetInstance().getCurrentBufferMasks(lev);
}

const iMultiFab*
WarpX::GatherBufferMasks (int lev)
{
    return GetInstance().getGatherBufferMasks(lev);
}

void
WarpX::StoreCurrent (int lev)
{
    for (int idim = 0; idim < 3; ++idim) {
        if (current_store[lev][idim]) {
            MultiFab::Copy(*current_store[lev][idim], *current_fp[lev][idim],
                           0, 0, 1, current_store[lev][idim]->nGrowVect());
        }
    }
}

void
WarpX::RestoreCurrent (int lev)
{
    for (int idim = 0; idim < 3; ++idim) {
        if (current_store[lev][idim]) {
            std::swap(current_fp[lev][idim], current_store[lev][idim]);
        }
    }
}

std::string
WarpX::Version ()
{
#ifdef WARPX_GIT_VERSION
    return std::string(WARPX_GIT_VERSION);
#else
    return std::string("Unknown");
#endif
}

std::string
WarpX::PicsarVersion ()
{
#ifdef PICSAR_GIT_VERSION
    return std::string(PICSAR_GIT_VERSION);
#else
    return std::string("Unknown");
#endif
}
