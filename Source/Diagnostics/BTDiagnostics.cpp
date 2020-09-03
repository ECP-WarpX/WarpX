#include "BTDiagnostics.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "ComputeDiagFunctors/BackTransformFunctor.H"
#include "ComputeDiagFunctors/RhoFunctor.H"
#include "Utils/CoarsenIO.H"
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
    // allocate vector of RealBost of the simulation domain in lab-frame
    m_prob_domain_lab.resize(m_num_buffers);
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
        // Define cell-centered multifab over the whole domain with
        // user-defined crse_ratio for nlevels
        DefineCellCenteredMultiFab(lev);
    }


}

void
BTDiagnostics::ReadParameters ()
{
    BaseReadParameters();
    auto & warpx = WarpX::GetInstance();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( ( warpx.do_back_transformed_diagnostics==true),
        "the do_back_transformed_diagnostics flag must be set to true for BTDiagnostics");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.gamma_boost > 1.0_rt,
        "gamma_boost must be > 1 to use the back-transformed diagnostics");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.boost_direction[2] == 1,
        "The back transformed diagnostics currently only works if the boost is in the z-direction");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.do_moving_window,
           "The moving window should be on if using the boosted frame diagnostic.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( warpx.moving_window_dir == AMREX_SPACEDIM-1,
           "The boosted frame diagnostic currently only works if the moving window is in the z direction.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "plotfile" || m_format == "openpmd",
        "<diag>.format must be plotfile or openpmd for back transformed diagnostics");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_crse_ratio == amrex::IntVect(1),
        "Only support for coarsening ratio of 1 in all directions is included for BTD\n"
        );

    // Read list of back-transform diag parameters requested by the user //
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
    snapshot_interval_is_specified = pp.query("dt_snapshots_lab", m_dt_snapshots_lab);
    if ( pp.query("dz_snapshots_lab", m_dz_snapshots_lab) ) {
        m_dt_snapshots_lab = m_dz_snapshots_lab/PhysConst::c;
        snapshot_interval_is_specified = true;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(snapshot_interval_is_specified,
        "For back-transformed diagnostics, user should specify either dz_snapshots_lab or dt_snapshots_lab");
    // For BTD, we always need rho to perform Lorentz Transform of current-density
    if (WarpXUtilStr::is_in(m_cellcenter_varnames, "rho")) warpx.setplot_rho(true);
}

void
BTDiagnostics::TMP_writeMetaData ()
{
    // This function will have the same functionality as writeMetaData in
    // previously used BackTransformedDiagnostics class to write
    // back-transformed data in a customized format

    if (amrex::ParallelDescriptor::IOProcessor()) {
        const std::string fullpath = m_file_prefix + "/snapshots";
        if (!amrex::UtilCreateDirectory(fullpath, 0755)) amrex::CreateDirectoryFailed(fullpath);

        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(m_file_prefix + "/snapshots/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if (!HeaderFile.good()) amrex::FileOpenFailed( HeaderFileName );

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
    // Return true if buffer is full or if force_flush == true
    // Return false if timestep < 0, i.e., at initialization when step == -1.
    if (step < 0 ) return false;
    else if ( buffer_full(i_buffer) || force_flush) {
        return true;
    }
    return false;
}


bool
BTDiagnostics::DoComputeAndPack (int step, bool /*force_flush*/)
{
    // always set to true for BTDiagnostics since back-transform buffers are potentially
    // computed and packed every timstep, except at initialization when step == -1.
    return (step>=0);
}

void
BTDiagnostics::InitializeFieldBufferData ( int i_buffer , int lev)
{
    auto & warpx = WarpX::GetInstance();
    // Lab-frame time for the i^th snapshot
    m_t_lab.at(i_buffer) = i_buffer * m_dt_snapshots_lab;


    // Compute lab-frame co-ordinates that correspond to the simulation domain
    // at level, lev, and time, m_t_lab[i_buffer] for each ith buffer.
    m_prob_domain_lab[i_buffer] = warpx.Geom(lev).ProbDomain();
    amrex::Real zmin_prob_domain_lab = m_prob_domain_lab[i_buffer].lo(m_moving_window_dir)
                                      / ( (1.0_rt + m_beta_boost) * m_gamma_boost);
    amrex::Real zmax_prob_domain_lab = m_prob_domain_lab[i_buffer].hi(m_moving_window_dir)
                                      / ( (1.0_rt + m_beta_boost) * m_gamma_boost);
    m_prob_domain_lab[i_buffer].setLo(m_moving_window_dir, zmin_prob_domain_lab +
                                               warpx.moving_window_v * m_t_lab[i_buffer] );
    m_prob_domain_lab[i_buffer].setHi(m_moving_window_dir, zmax_prob_domain_lab +
                                               warpx.moving_window_v * m_t_lab[i_buffer] );


    // Define buffer domain in boosted frame at level, lev, with user-defined lo and hi
    amrex::RealBox diag_dom;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim ) {
        // Setting lo-coordinate for the diag domain by taking the max of user-defined
        // lo-cordinate and lo-coordinat of the simulation domain at level, lev
        diag_dom.setLo(idim, std::max(m_lo[idim],warpx.Geom(lev).ProbLo(idim)) );
        // Setting hi-coordinate for the diag domain by taking the max of user-defined
        // hi-cordinate and hi-coordinate of the simulation domain at level, lev
        diag_dom.setHi(idim, std::min(m_hi[idim],warpx.Geom(lev).ProbHi(idim)) );
    }
    // Initializing the m_buffer_box for the i^th snapshot.
    // At initialization, the Box has the same index space as the boosted-frame
    // As time-progresses, the z-dimension indices will be modified based on
    // current_z_lab
    amrex::IntVect lo(0);
    amrex::IntVect hi(-1);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        // lo index with same cell-size as simulation at level, lev.
        const int lo_index = static_cast<int>( floor(
                ( diag_dom.lo(idim) - warpx.Geom(lev).ProbLo(idim) ) /
                  warpx.Geom(lev).CellSize(idim) ) );
        // Taking max of (0,lo_index) because lo_index must always be >=0
        lo[idim] = std::max( 0, lo_index );
        // hi index with same cell-size as simulation at level, lev.
        const int hi_index =  static_cast<int>( ceil(
                ( diag_dom.hi(idim) - warpx.Geom(lev).ProbLo(idim) ) /
                  warpx.Geom(lev).CellSize(idim) ) );
        // Taking max of (0,hi_index) because hi_index must always be >=0
        // Subtracting by 1 because lo,hi indices are set to cell-centered staggering.
        hi[idim] = std::max( 0, hi_index) - 1;
        // if hi<=lo, then hi = lo + 1, to ensure one cell in that dimension
        if ( hi[idim] <= lo[idim] ) {
             hi[idim]  = lo[idim] + 1;
             AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_crse_ratio[idim]==1, "coarsening ratio in reduced dimension must be 1."
             );
        }
    }
    amrex::Box diag_box( lo, hi );
    m_buffer_box[i_buffer] = diag_box;

    // Define buffer_domain in lab-frame for the i^th snapshot.
    // Replace z-dimension with lab-frame co-ordinates.
    amrex::Real zmin_buffer_lab = diag_dom.lo(m_moving_window_dir)
                                 / ( (1.0_rt + m_beta_boost) * m_gamma_boost);
    amrex::Real zmax_buffer_lab = diag_dom.hi(m_moving_window_dir)
                                 / ( (1.0_rt + m_beta_boost) * m_gamma_boost);


    m_buffer_domain_lab[i_buffer] = diag_dom;
    m_buffer_domain_lab[i_buffer].setLo(m_moving_window_dir,
                                  zmin_buffer_lab + warpx.moving_window_v * m_t_lab[i_buffer]);
    m_buffer_domain_lab[i_buffer].setHi(m_moving_window_dir,
                                  zmax_buffer_lab + warpx.moving_window_v * m_t_lab[i_buffer]);

    // Initialize buffer counter and z-positions of the  i^th snapshot in
    // boosted-frame and lab-frame
    m_buffer_counter[i_buffer] = 0;
    m_current_z_lab[i_buffer] = 0.0_rt;
    m_current_z_boost[i_buffer] = 0.0_rt;
    // Now Update Current Z Positions
    m_current_z_boost[i_buffer] = UpdateCurrentZBoostCoordinate(m_t_lab[i_buffer],
                                                              warpx.gett_new(lev) );
    m_current_z_lab[i_buffer] = UpdateCurrentZLabCoordinate(m_t_lab[i_buffer],
                                                              warpx.gett_new(lev) );


    // Compute number of cells in lab-frame required for writing Header file
    // and potentially to generate Back-Transform geometry to ensure
    // compatibility with plotfiles.
    // For the z-dimension, number of cells in the lab-frame is
    // computed using the coarsened cell-size in the lab-frame obtained using
    // the ref_ratio at level, lev-1.
    amrex::IntVect ref_ratio = amrex::IntVect(1);
    if (lev > 0 ) ref_ratio = WarpX::RefRatio(lev-1);
    // Number of lab-frame cells in z-direction at level, lev
    const int num_zcells_lab = static_cast<int>( floor (
                                   ( zmax_buffer_lab - zmin_buffer_lab)
                                   / dz_lab(warpx.getdt(lev), ref_ratio[m_moving_window_dir])                               ) );
    // Take the max of 0 and num_zcells_lab
    int Nz_lab = std::max( 0, num_zcells_lab );
    // Number of lab-frame cells in x-direction at level, lev
    const int num_xcells_lab = static_cast<int>( floor (
                                  ( diag_dom.hi(0) - diag_dom.lo(0) )
                                  / warpx.Geom(lev).CellSize(0)
                              ) );
    // Take the max of 0 and num_ycells_lab
    int Nx_lab = std::max( 0, num_xcells_lab);
#if (AMREX_SPACEDIM == 3)
    // Number of lab-frame cells in the y-direction at level, lev
    const int num_ycells_lab = static_cast<int>( floor (
                                   ( diag_dom.hi(1) - diag_dom.lo(1) )
                                   / warpx.Geom(lev).CellSize(1)
                               ) );
    // Take the max of 0 and num_xcells_lab
    int Ny_lab = std::max( 0, num_ycells_lab );
    m_buffer_ncells_lab[i_buffer] = {Nx_lab, Ny_lab, Nz_lab};
#else
    m_buffer_ncells_lab[i_buffer] = {Nx_lab, Nz_lab};
#endif
    // Call funtion to create directories for customized output format
    TMP_createLabFrameDirectories(i_buffer, lev);
}

void
BTDiagnostics::DefineCellCenteredMultiFab(int lev)
{
    // Creating MultiFab to store cell-centered data in boosted-frame for the entire-domain
    // This MultiFab will store all the user-requested fields in the boosted-frame
    auto & warpx = WarpX::GetInstance();
    // The BoxArray is coarsened based on the user-defined coarsening ratio
    amrex::BoxArray ba = warpx.boxArray(lev);
    ba.coarsen(m_crse_ratio);
    amrex::DistributionMapping dmap = warpx.DistributionMap(lev);
    int ngrow = 1;
    m_cell_centered_data[lev].reset( new amrex::MultiFab(ba, dmap,
                                     m_cellcenter_varnames.size(), ngrow) );

}

void
BTDiagnostics::InitializeFieldFunctors (int lev)
{
    auto & warpx = WarpX::GetInstance();
    // Clear any pre-existing vector to release stored data
    // This ensures that when domain is load-balanced, the functors point
    // to the correct field-data pointers
    m_all_field_functors[lev].clear();
    // For back-transformed data, all the components are cell-centered and stored
    // in a single multifab, m_cell_centered_data.
    // Therefore, size of functors at all levels is 1.
    int num_BT_functors = 1;
    m_all_field_functors[lev].resize(num_BT_functors);
    m_cell_center_functors[lev].clear();
    m_cell_center_functors[lev].resize( m_cellcenter_varnames.size() );
    // Create an object of class BackTransformFunctor
    for (int i = 0; i < num_BT_functors; ++i)
    {
        // coarsening ratio is not provided since the source MultiFab, m_cell_centered_data
        // is coarsened based on the user-defined m_crse_ratio
        m_all_field_functors[lev][i] = std::make_unique<BackTransformFunctor>(
                  m_cell_centered_data[lev].get(), lev,
                  m_varnames.size(), m_num_buffers, m_varnames);
    }

    // Define all cell-centered functors required to compute cell-centere data
    // Fill vector of cell-center functors for all field-components, namely,
    // Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, and rho are included in the
    // cell-center functors for BackTransform Diags
    for (int comp=0, n=m_cell_center_functors[lev].size(); comp<n; comp++){
        if        ( m_cellcenter_varnames[comp] == "Ex" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 0), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "Ey" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 1), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "Ez" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev, 2), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "Bx" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 0), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "By" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 1), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "Bz" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Bfield_aux(lev, 2), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "jx" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 0), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "jy" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 1), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "jz" ){
            m_cell_center_functors[lev][comp] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_current_fp(lev, 2), lev, m_crse_ratio);
        } else if ( m_cellcenter_varnames[comp] == "rho" ){
            m_cell_center_functors[lev][comp] = std::make_unique<RhoFunctor>(lev, m_crse_ratio);
        }
    }

}

// Temporary function only to debug the current implementation.
// Will be replaced with plotfile/OpenPMD functionality
void
BTDiagnostics::TMP_createLabFrameDirectories(int i_buffer, int lev)
{
    // This is a creates lab-frame directories for writing out customized BTD output.
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if ( !amrex::UtilCreateDirectory (m_file_name[i_buffer], 0755) )
            amrex::CreateDirectoryFailed(m_file_name[i_buffer]);

        const std::string &fullpath = amrex::LevelFullPath(lev, m_file_name[i_buffer]);
        if ( !amrex::UtilCreateDirectory(fullpath, 0755) )
            amrex::CreateDirectoryFailed(fullpath);
    }
    amrex::ParallelDescriptor::Barrier();
    TMP_writeLabFrameHeader(i_buffer);

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
        }
        // Check that the proper number of user-requested components are cell-centered
        AMREX_ALWAYS_ASSERT( icomp_dst == m_cellcenter_varnames.size() );
        // fill boundary call is required to average_down (flatten) data to
        // the coarsest level.
        m_cell_centered_data[lev]->FillBoundary(warpx.Geom(lev).periodicity() );
    }
    // Flattening out MF over levels

    for (int lev = warpx.finestLevel(); lev > 0; --lev) {
        CoarsenIO::Coarsen( *m_cell_centered_data[lev-1], *m_cell_centered_data[lev], 0, 0,
                             m_cellcenter_varnames.size(), 0, WarpX::RefRatio(lev-1) );
    }

    int num_BT_functors = 1;

    for (int lev = 0; lev < nlev_output; ++lev)
    {
        for (int i = 0; i < num_BT_functors; ++i)
        {
            for (int i_buffer = 0; i_buffer < m_num_buffers; ++i_buffer )
            {
                // Update z-boost and z-lab positions
                m_current_z_boost[i_buffer] = UpdateCurrentZBoostCoordinate(m_t_lab[i_buffer],
                                                                      warpx.gett_new(lev) );
                m_current_z_lab[i_buffer] = UpdateCurrentZLabCoordinate(m_t_lab[i_buffer],
                                                                      warpx.gett_new(lev) );
                bool ZSliceInDomain = GetZSliceInDomainFlag (i_buffer, lev);
                // Initialize and define field buffer multifab if buffer is empty
                if (ZSliceInDomain) {
                    if ( buffer_empty(i_buffer) ) DefineFieldBufferMultiFab(i_buffer, lev);
                }
                m_all_field_functors[lev][i]->PrepareFunctorData (
                                             i_buffer, ZSliceInDomain,
                                             m_current_z_boost[i_buffer],
                                             m_buffer_box[i_buffer],
                                             k_index_zlab(i_buffer, lev) );

                if (ZSliceInDomain) ++m_buffer_counter[i_buffer];
            }
        }
    }

}


amrex::Real
BTDiagnostics::dz_lab (amrex::Real dt, amrex::Real ref_ratio){
    return PhysConst::c * dt * 1._rt/m_beta_boost * 1._rt/m_gamma_boost * 1._rt/ref_ratio;
}


int
BTDiagnostics::k_index_zlab (int i_buffer, int lev)
{
    auto & warpx = WarpX::GetInstance();
    amrex::Real prob_domain_zmin_lab = m_prob_domain_lab[i_buffer].lo( m_moving_window_dir );
    amrex::IntVect ref_ratio = amrex::IntVect(1);
    if (lev > 0 ) ref_ratio = WarpX::RefRatio(lev-1);
    int k_lab = static_cast<unsigned>( (
                          ( m_current_z_lab[i_buffer] - prob_domain_zmin_lab )
                          / dz_lab( warpx.getdt(lev), ref_ratio[m_moving_window_dir] )
                      ) );
    return k_lab;
}



void
BTDiagnostics::DefineFieldBufferMultiFab (const int i_buffer, const int lev)
{
    if ( m_do_back_transformed_fields ) {

        const int k_lab = k_index_zlab (i_buffer, lev);
        m_buffer_box[i_buffer].setSmall( m_moving_window_dir, k_lab - m_buffer_size + 1);
        m_buffer_box[i_buffer].setBig( m_moving_window_dir, k_lab);

        amrex::BoxArray buffer_ba( m_buffer_box[i_buffer] );
        buffer_ba.maxSize(m_max_box_size);

        // Generate a new distribution map for the back-transformed buffer multifab
        amrex::DistributionMapping buffer_dmap(buffer_ba);
        // Number of guard cells for the output buffer is zero.
        // Unlike FullDiagnostics, "m_format == sensei" option is not included here.
        int ngrow = 0;
        m_mf_output[i_buffer][lev] = amrex::MultiFab ( buffer_ba, buffer_dmap,
                                                  m_varnames.size(), ngrow ) ;
    }
}

bool
BTDiagnostics::GetZSliceInDomainFlag (const int i_buffer, const int lev)
{
    auto & warpx = WarpX::GetInstance();
    const amrex::RealBox& boost_domain = warpx.Geom(lev).ProbDomain();

    amrex::Real buffer_zmin_lab = m_buffer_domain_lab[i_buffer].lo( m_moving_window_dir );
    amrex::Real buffer_zmax_lab = m_buffer_domain_lab[i_buffer].hi( m_moving_window_dir );

    if ( ( m_current_z_boost[i_buffer] < boost_domain.lo(m_moving_window_dir) ) or
         ( m_current_z_boost[i_buffer] > boost_domain.hi(m_moving_window_dir) ) or
         ( m_current_z_lab[i_buffer] < buffer_zmin_lab ) or
         ( m_current_z_lab[i_buffer] > buffer_zmax_lab ) )
    {
        // the slice is not in the boosted domain or lab-frame domain
        return false;
    }

    return true;
}

void
BTDiagnostics::TMP_writeLabFrameHeader (int i_buffer)
{
    if (amrex::ParallelDescriptor::IOProcessor() )
    {
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(m_file_name[i_buffer] + "/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out    |
                                                std::ofstream::trunc  |
                                                std::ofstream::binary );

        if ( !HeaderFile.good() ) amrex::FileOpenFailed( HeaderFileName );

        HeaderFile.precision(17);

        HeaderFile << m_t_lab[i_buffer] << "\n";
        // Write buffer number of cells
        HeaderFile << m_buffer_ncells_lab[i_buffer][0] << ' '
#if ( AMREX_SPACEDIM==3 )
                   << m_buffer_ncells_lab[i_buffer][1] << ' '
#endif
                   << m_buffer_ncells_lab[i_buffer][m_moving_window_dir] << "\n";
        // Write physical boundary of the buffer
        // Lower bound
        HeaderFile << m_buffer_domain_lab[i_buffer].lo(0) << ' '
#if ( AMREX_SPACEDIM == 3 )
                   << m_buffer_domain_lab[i_buffer].lo(1) << ' '
#endif
                   << m_buffer_domain_lab[i_buffer].lo(m_moving_window_dir) << "\n";
        // Higher bound
        HeaderFile << m_buffer_domain_lab[i_buffer].hi(0) << ' '
#if ( AMREX_SPACEDIM == 3 )
                   << m_buffer_domain_lab[i_buffer].hi(1) << ' '
#endif
                   << m_buffer_domain_lab[i_buffer].hi(m_moving_window_dir) << "\n";
        // List of fields to flush to file
        for (int i = 0; i < m_varnames.size(); ++i)
        {
            HeaderFile << m_varnames[i] << ' ';
        }
        HeaderFile << "\n";
    }
}

void
BTDiagnostics::Flush (int i_buffer)
{
    TMP_FlushLabFrameData (i_buffer);
    // Reset the buffer counter to zero after flushing out data stored in the buffer.
    ResetBufferCounter(i_buffer);
}

void
BTDiagnostics::TMP_FlushLabFrameData ( int i_buffer )
{
    // customized output format for writing data stored in buffers to the disk.
    amrex::VisMF::Header::Version current_version = amrex::VisMF::GetHeaderVersion();
    amrex::VisMF::SetHeaderVersion(amrex::VisMF::Header::NoFabHeader_v1);

    if ( ! buffer_empty(i_buffer) ) {
        for ( int lev = 0; lev < nlev_output; ++lev) {
            const int k_lab = k_index_zlab (i_buffer, lev);
            const amrex::BoxArray& ba = m_mf_output[i_buffer][lev].boxArray();
            const int hi = ba[0].bigEnd(m_moving_window_dir);
            const int lo = hi - m_buffer_counter[i_buffer] + 1;

            amrex::Box buffer_box = m_buffer_box[i_buffer];
            buffer_box.setSmall(m_moving_window_dir, lo);
            buffer_box.setBig(m_moving_window_dir, hi);
            amrex::BoxArray buffer_ba(buffer_box);
            buffer_ba.maxSize(m_max_box_size);
            amrex::DistributionMapping buffer_dm(buffer_ba);

            amrex::MultiFab tmp (buffer_ba, buffer_dm, m_varnames.size(), 0);
            tmp.copy(m_mf_output[i_buffer][lev], 0, 0, m_varnames.size());

            std::stringstream ss;
            ss << m_file_name[i_buffer] << "/Level_0/"
               << amrex::Concatenate("buffer", k_lab, 5);
            amrex::VisMF::Write(tmp, ss.str());

       }
    }
    amrex::VisMF::SetHeaderVersion(current_version);
}
