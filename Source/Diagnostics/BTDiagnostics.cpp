#include "BTDiagnostics.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/BackTransformFunctor.H"
#include "Utils/CoarsenIO.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
using namespace amrex::literals;
using namespace amrex;

BTDiagnostics::BTDiagnostics (int i, std::string name)
    : Diagnostics(i, name)
{
    ReadParameters();
}


void BTDiagnostics::InitDerivedData ()
{
    // storing BTD related member variables
    auto & warpx = WarpX::GetInstance();
    m_gamma_boost = WarpX::gamma_boost;
    m_beta_boost = std::sqrt( 1._rt - 1._rt/( m_gamma_boost * m_gamma_boost) );
    m_moving_window_dir = warpx.moving_window_dir;

    // temporary function related to customized output from previous BTD to verify accuracy
    writeMetaData();

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
        // temporary variable name for customized BTD output to verify accuracy
        m_file_name[i] = amrex::Concatenate(m_file_prefix +"/snapshots/snapshot",i,5);
    }
    for (int lev = 0; lev < nmax_lev; ++lev) {
        // Define cell-centered multifab over the whole domain with
        // user-defined crse_ratio for nlevels
        DefineCellCenteredMultiFab(lev);
    }

}

void
BTDiagnostics::ReadParameters ()
{
    ReadBaseParameters();
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
    bool snapshot_interval_is_specified = 0;
    amrex::Real m_dz_snapshots_lab = 0.0_rt;
    snapshot_interval_is_specified = pp.query("dt_snapshots_lab", m_dt_snapshots_lab);
    if ( pp.query("dz_snapshots_lab", m_dz_snapshots_lab) ) {
        m_dt_snapshots_lab = m_dz_snapshots_lab/PhysConst::c;
        snapshot_interval_is_specified = 1;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(snapshot_interval_is_specified,
        "For back-transformed diagnostics, user should specify either dz_snapshots_lab or dt_snapshots_lab");

}

void
BTDiagnostics::writeMetaData ()
{
    // This function will have the same functionality as writeMetaData in
    // previously used BackTransformedDiagnostics class to write
    // back-transformed data in a customized format

    if (ParallelDescriptor::IOProcessor()) {
        const std::string fullpath = m_file_prefix + "/snapshots";
        if (!UtilCreateDirectory(fullpath, 0755)) CreateDirectoryFailed(fullpath);

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(m_file_prefix + "/snapshots/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if (!HeaderFile.good()) FileOpenFailed( HeaderFileName );

        HeaderFile.precision(17);
        HeaderFile << m_num_snapshots_lab << "\n";
        HeaderFile << m_dt_snapshots_lab << "\n";
        HeaderFile << m_gamma_boost << "\n";
        HeaderFile << m_beta_boost << "\n";

    }

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
    m_t_lab[i_buffer] = i_buffer * m_dt_snapshots_lab;
    // 2. Define domain in boosted frame at level, lev
    amrex::RealBox diag_dom;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim ) {
        diag_dom.setLo(idim, max(m_lo[idim],warpx.Geom(lev).ProbLo(idim)) );
        diag_dom.setHi(idim, min(m_hi[idim],warpx.Geom(lev).ProbHi(idim)) );
    }
    // 3. Initializing the m_buffer_box for the i^th snapshot.
    //    At initialization, the Box has the same index space as the boosted-frame
    //    As time-progresses, the z-dimension indices will be modified based on
    //    current_z_lab
    IntVect lo(0);
    IntVect hi(1);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        // lo index with same cell-size as simulation at level, lev.
        lo[idim] = max( static_cast<int>( floor (
                      ( diag_dom.lo(idim) - warpx.Geom(lev).ProbLo(idim)) /
                        warpx.Geom(lev).CellSize(idim)) ), 0 );
        // hi index with same cell-size as simulation at level, lev.
        hi[idim] = max( static_cast<int> ( ceil (
                      ( diag_dom.hi(idim) - warpx.Geom(lev).ProbLo(idim)) /
                        warpx.Geom(lev).CellSize(idim) ) ), 0) - 1 ;
        // if hi<=lo, then hi = lo + 1, to ensure one cell in that dimension
        if ( hi[idim] <= lo[idim] ) {
             hi[idim]  = lo[idim] + 1;
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_crse_ratio[idim]==1, "coarsening ratio in reduced dimension must be 1."
             );
        }
    }
    Box diag_box( lo, hi );
    // The box is not coarsened yet. Should this be coarsened here or when
    // the buffer multifab is initialized?
    m_buffer_box[i_buffer] = diag_box;

    // 4. Define buffer_domain  in lab-frame for the i^th snapshot.
    //    Replace z-dimension with lab-frame co-ordinates.
    amrex::Real zmin_lab = diag_dom.lo(m_boost_direction)
                           / ( (1.0_rt + m_beta_boost) * m_gamma_boost);
    amrex::Real zmax_lab = diag_dom.hi(m_boost_direction)
                           / ( (1.0_rt + m_beta_boost) * m_gamma_boost);

    m_buffer_domain_lab[i_buffer] = warpx.Geom(lev).ProbDomain();
    m_buffer_domain_lab[i_buffer].setLo(m_boost_direction, zmin_lab + warpx.moving_window_v * m_t_lab[i_buffer]);
    m_buffer_domain_lab[i_buffer].setHi(m_boost_direction, zmin_lab + warpx.moving_window_v * m_t_lab[i_buffer]);

    // 5. Initialize buffer counter and z-positions of the  i^th snapshot in
    //    boosted-frame and lab-frame
    m_buffer_counter[i_buffer] = 0;
    m_current_z_lab[i_buffer] = 0.0;
    m_current_z_boost[i_buffer] = 0.0;
    // Now Update Current Z Positions
    m_current_z_boost[i_buffer] = UpdateCurrentZBoostCoordinate(m_t_lab[i_buffer],
                                                              warpx.gett_new(lev) );
    m_current_z_lab[i_buffer] = UpdateCurrentZLabCoordinate(m_t_lab[i_buffer],
                                                              warpx.gett_new(lev) );


    // 6. Compute ncells_lab required for writing Header file and potentially to generate
    //    Back-Transform geometry to ensure compatibility with plotfiles //
    IntVect ref_ratio = WarpX::RefRatio(lev);
    int Nz_lab = max( static_cast<int>( floor ( ( zmax_lab - zmin_lab)
                      / dz_lab(warpx.getdt(lev), ref_ratio[AMREX_SPACEDIM-1]) ) ), 0);
    int Nx_lab = max( static_cast<int>( floor( (diag_dom.hi(0) - diag_dom.lo(0) )
                      / warpx.Geom(lev).CellSize(0) )), 0);
#if (AMREX_SPACEDIM == 3)
    int Ny_lab = max( static_cast<int>( floor( (diag_dom.hi(1) - diag_dom.lo(1) )
                      / warpx.Geom(lev).CellSize(1) )), 0);
    m_buffer_ncells_lab[i_buffer] = {Nx_lab, Ny_lab, Nz_lab};
#else
    int Ny_lab = 0;
    m_buffer_ncells_lab[i_buffer] = {Nx_lab, Nz_lab};
#endif

    // 7. Call funtion to create directories for customized output format
    createLabFrameDirectories(i_buffer, lev);
}

void
BTDiagnostics::DefineCellCenteredMultiFab(int lev)
{
    // Creating MultiFab to store cell-centered data in boosted-frame for the entire-domain
    // This MultiFab will store all the user-requested fields in the boosted-frame
    auto & warpx = WarpX::GetInstance();
    // The BoxArray is coarsened based on the user-defined coarsening ratio
    BoxArray ba = warpx.boxArray(lev);
    ba.coarsen(m_crse_ratio);
    DistributionMapping dmap = warpx.DistributionMap(lev);
    int ngrow = 1;
    m_cell_centered_data[lev].reset( new amrex::MultiFab(ba, dmap, m_varnames.size(), ngrow) );

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
    for (int i = 0; i < num_BT_functors; ++i)
    {
        // coarsening ratio is not provided since the source MultiFab, m_cell_centered_data
        // is coarsened based on the user-defined m_crse_ratio
        m_all_field_functors[lev][i] = std::make_unique<BackTransformFunctor>(
                  m_cell_centered_data[lev].get(), lev, m_varnames.size() );
    }

    // 2. Define all cell-centered functors required to compute cell-centere data
    //    Fill vector of cell-center functors for all components
    //    Only E,B, j, and rho are included in the cell-center functors for BackTransform Diags
    for (int comp=0, n=m_cell_center_functors[lev].size(); comp<n; comp++){
        if        ( m_varnames[comp] == "Ex" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 0), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Ey" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 1), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Ez" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 2), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Bx" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 0), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "By" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 1), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "Bz" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 2), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "jx" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 0), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "jy" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 1), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "jz" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 2), lev, m_crse_ratio);
        } else if ( m_varnames[comp] == "rho" ){
            // rho_new is stored in component 1 of rho_fp when using PSATD
#ifdef WARPX_USE_PSATD
            MultiFab* rho_new = new MultiFab(*warpx.get_pointer_rho_fp(lev), amrex::make_alias, 1, 1);
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(rho_new, lev, m_crse_ratio);
#else
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_rho_fp(lev), lev, m_crse_ratio);
#endif
        }
    }

}

// Temporary function only to debug the current implementation.
// Will be replaced with plotfile/OpenPMD functionality
void
BTDiagnostics::createLabFrameDirectories(int i_buffer, int lev)
{
    // This function will include relevant code from class BackTransformedDiagnostics
    // to create lab-frame directories. Note that here we will also add level, lev
    // if required.

}

void
BTDiagnostics::PrepareFieldDataForOutput ()
{
    auto & warpx = WarpX::GetInstance();
    // In this function, we will get cell-centered data for every level, lev,
    // using the cell-center functors and their respective opeators()
    // Call m_cell_center_functors->operator
    for (int lev = 0; lev < nmax_lev; ++lev) {
        int icomp_dst = 0;
        for (int icomp = 0, n=m_cell_center_functors[0].size(); icomp<n; ++icomp) {
            // Call all the cell-center functors in m_cell_center_functors.
            // Each of them computes cell-centered data for a field and
            // stores it in cell-centered MultiFab, m_cell_centered_data[lev].
            m_cell_center_functors[lev][icomp]->operator()(*m_cell_centered_data[lev], icomp_dst);
            icomp_dst += m_cell_center_functors[lev][icomp]->nComp();
//            m_cell_centered_data[lev]->FillBoundary(warpx.Geom(lev).periodocity() );
        }
        // Check that the proper number of user-requested components are cell-centered
        AMREX_ALWAYS_ASSERT( icomp_dst == m_varnames.size() );
        // fill boundary call is required to average_down (flatten) data to 
        // the coarsest level.
        m_cell_centered_data[lev]->FillBoundary(warpx.Geom(lev).periodicity() );
    }
    // Flattening out MF over levels -- this will be moved to a general place
    // so that it can be called for FullDiagnostics as well.
    for (int lev = nmax_lev; lev > 0; --lev) {
        CoarsenIO::Coarsen( *m_cell_centered_data[lev-1], *m_cell_centered_data[lev], 0, 0, m_varnames.size(), 0, WarpX::RefRatio(lev-1) );
    }
}

