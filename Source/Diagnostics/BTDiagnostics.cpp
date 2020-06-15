#include "BTDiagnostics.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
using namespace amrex::literals;

BTDiagnostics::BTDiagnostics (int i, std::string name)
    : Diagnostics(i, name)
{
    ReadParameters();
}


void BTDiagnostics::DerivedInitData ()
{
    // storing BTD related member variables
    auto & warpx = WarpX::GetInstance();
    m_gamma_boost = WarpX::gamma_boost;
    m_beta_boost = std::sqrt( 1._rt - 1._rt/( m_gamma_boost * m_gamma_boost) );
    m_moving_window_dir = warpx.moving_window_dir;
    // Currently, for BTD, all the data is averaged+coarsened to coarsest level
    // and then sliced+back-transformed+filled_to_buffer.
    // The number of levels to be output is nlev_output.
    nlev_output = 1;
    // temporary function related to customized output from previous BTD to verify accuracy
    TMP_writeMetaData();

    // allocate vector of m_t_lab with m_num_buffers;
    m_t_lab.resize(m_num_buffers);
    // allocate vector of RealBox of the diag domain
    m_buffer_domain_lab.resize(m_num_buffers);
    // define box correctly (one for all snapshots)
    m_buffer_box.resize(m_num_buffers);
    // allocate vector of m_current_z_lab
    m_current_z_lab.resize(m_num_buffers);
    // allocate vector of m_num_buffers
    m_current_z_boost.resize(m_num_buffers);
    // allocate vector of m_buff_counter
    m_buffer_counter.resize(m_num_buffers);
    // allocate vector of num_Cells in the lab-frame
    m_buffer_ncells_lab.resize(m_num_buffers);
    // allocate vector of file names for each buffer
    m_file_name.resize(m_num_buffers);
    // allocate vector of cell centered multifabs for nlevels
    m_cell_centered_data.resize(nmax_lev);
    // allocate vector of cell-center functors for nlevels
    m_cell_center_functors.resize(nmax_lev);

    for (int i = 0; i < m_num_buffers; ++i) {
        // TMP variable name for customized BTD output to verify accuracy
        m_file_name[i] = amrex::Concatenate(m_file_prefix +"/snapshots/snapshot",i,5);
    }
    for (int lev = 0; lev < nmax_lev; ++lev) {
        // Define cell-centered multifab over the whole domain with user-defined crse_ratio
        // for nlevels
        DefineCellCenteredMultiFab(lev);
    }

}

void
BTDiagnostics::ReadParameters ()
{
    BaseReadParameters();
    auto & warpx = WarpX::GetInstance();
    // Read list of back-transform diag parameters requested by the user //
    amrex::Print() << " in read parameters for BTD \n";

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ( warpx.do_back_transformed_diagnostics==true),
        "the do_back_transformed_diagnostics flag must be set to true for BTDiagnostics");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.gamma_boost > 1.0_rt,
        "gamma_boost must be > 1 to use the back-transformed diagnostics");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.boost_direction[2] == 1,
        "The back transformed diagnostics currently only works if the boost is in the z-direction");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.do_moving_window,
           "The moving window should be on if using the boosted frame diagnostic.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.moving_window_dir == 2,
           "The boosted frame diagnostic currently only works if the moving window is in the z direction.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "plotfile" || m_format == "openpmd",
        "<diag>.format must be plotfile or openpmd for back transformed diagnostics");
    amrex::ParmParse pp(m_diag_name);

    m_file_prefix = "diags/lab_frame_data/" + m_diag_name;
    pp.query("file_prefix", m_file_prefix);
    pp.query("do_back_transformed_fields", m_do_back_transformed_fields);
    pp.query("do_back_transformed_particles", m_do_back_transformed_particles);
    AMREX_ALWAYS_ASSERT(m_do_back_transformed_fields or m_do_back_transformed_particles);

    pp.get("num_snapshots_lab", m_num_snapshots_lab);
    m_num_buffers = m_num_snapshots_lab;

    // Read either dz_snapshots_lab or dt_snapshots_lab
    bool snapshot_interval_is_specified = false;
    amrex::Real m_dz_snapshots_lab = 0.0_rt;
    snapshot_interval_is_specified = pp.query("dt_snapshots_lab", m_dt_snapshots_lab);
    if ( pp.query("dz_snapshots_lab", m_dz_snapshots_lab) ) {
        m_dt_snapshots_lab = m_dz_snapshots_lab/PhysConst::c;
        snapshot_interval_is_specified = true;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(snapshot_interval_is_specified,
        "For back-transformed diagnostics, user should specify either dz_snapshots_lab or dt_snapshots_lab");

}

void
BTDiagnostics::TMP_writeMetaData ()
{
    // This function will have the same functionality as writeMetaData in
    // previously used BackTransformedDiagnostics class to write
    // back-transformed data in a customized format
}

bool
BTDiagnostics::DoDump (int step, int i_buffer, bool force_flush)
{
    // check if buffer is full (m_buffer_counter[i_buffer] == m_) or if force_flush == true
    if ( buffer_full(i_buffer) || force_flush) {
        return true;
    }
    return false;
}


bool
BTDiagnostics::DoComputeAndPack (int step, bool force_flush)
{
    // always set to true for BTDiagnostics since back-transform buffers are potentially
    // computed and packed every timstep.
    return true;
}

void
BTDiagnostics::InitializeFieldBufferData ( int i_buffer , int lev)
{
    auto & warpx = WarpX::GetInstance();
    // 1. Lab-frame time for the i^th snapshot
    // 2. Define domain in boosted frame at level, lev
    // 3. Initializing the m_buffer_box for the i^th snapshot.
    //    At initialization, the Box has the same index space as the boosted-frame
    //    As time-progresses, the z-dimension indices will be modified based on
    //    current_z_lab
    // 4. Define buffer_domain  in lab-frame for the i^th snapshot.
    //    Replace z-dimension with lab-frame co-ordinates.
    // 5. Initialize buffer counter and z-positions of the  i^th snapshot in
    //    boosted-frame and lab-frame
    // 6. Compute ncells_lab required for writing Header file and potentially to generate
    //    Back-Transform geometry to ensure compatibility with plotfiles //
    // 7. Call funtion to create directories for customized output format
    TMP_createLabFrameDirectories(i_buffer, lev);
}

void
BTDiagnostics::DefineCellCenteredMultiFab(int lev)
{
    // Creating MultiFab to store cell-centered data in boosted-frame for the entire-domain
    // This MultiFab will store all the user-requested fields in the boosted-frame
    auto & warpx = WarpX::GetInstance();
    // The BoxArray is coarsened based on the user-defined coarsening ratio
}

void
BTDiagnostics::InitializeFieldFunctors (int lev)
{
    auto & warpx = WarpX::GetInstance();
    // Clear any pre-existing vector to release stored data
    // This ensures that when domain is load-balanced, the functors point
    // to the correct field-data pointers
    m_all_field_functors[lev].clear();
    // For back-transformed data, all the components are cell-centered and stored in a single multifab. Therefore, size of functors at all levels is 1.
    int num_BT_functors = 1;
    m_all_field_functors[lev].resize(num_BT_functors);
    m_cell_center_functors[lev].clear();
    m_cell_center_functors[lev].resize( m_varnames.size() );
    // 1. create an object of class BackTransformFunctor

    // 2. Define all cell-centered functors required to compute cell-centere data
    //    Fill vector of cell-center functors for all components
    //    Only E,B, j, and rho are included in the cell-center functors for BackTransform Diags
}

// Temporary function only to debug the current implementation.
// Will be replaced with plotfile/OpenPMD functionality
void
BTDiagnostics::TMP_createLabFrameDirectories(int i_buffer, int lev)
{
    // This function will include relevant code from class BackTransformedDiagnostics
    // to create lab-frame directories. Note that here we will also add level, lev
    // if required.

}

void
BTDiagnostics::PrepareFieldDataForOutput ()
{
    // In this function, we will get cell-centered data for every level, lev,
    // using the cell-center functors and their respective opeators()
    // Call m_cell_center_functors->operator
}


amrex::Real
BTDiagnostics::dz_lab (amrex::Real dt, amrex::Real ref_ratio)
{
    return PhysConst::c * dt * 1._rt/m_beta_boost * 1._rt/m_gamma_boost * 1._rt/ref_ratio;
}
