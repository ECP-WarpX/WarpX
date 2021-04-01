#include "BTDiagnostics.H"
#include "BTD_Plotfile_Header_Impl.H"
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

#include <memory>

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

    // allocate vector of m_t_lab with m_num_buffers;
    m_t_lab.resize(m_num_buffers);
    // allocate vector of RealBost of the simulation domain in lab-frame
    m_prob_domain_lab.resize(m_num_buffers);
    // allocate vector of RealBox of the diag domain
    m_snapshot_domain_lab.resize(m_num_buffers);
    // allocate vector of RealBox of the buffers that fill the snapshot
    m_buffer_domain_lab.resize(m_num_buffers);
    // define box correctly (one for all snapshots)
    m_snapshot_box.resize(m_num_buffers);
    // define box for each buffer that fills the snapshot
    m_buffer_box.resize(m_num_buffers);
    // allocate vector of m_current_z_lab
    m_current_z_lab.resize(m_num_buffers);
    // allocate vector of m_num_buffers
    m_current_z_boost.resize(m_num_buffers);
    // allocate vector of m_buff_counter to counter number of slices filled in the buffer
    m_buffer_counter.resize(m_num_buffers);
    // allocate vector of num_Cells in the lab-frame
    m_snapshot_ncells_lab.resize(m_num_buffers);
    // allocate vector of cell centered multifabs for nlevels
    m_cell_centered_data.resize(nmax_lev);
    // allocate vector of cell-center functors for nlevels
    m_cell_center_functors.resize(nmax_lev);
    // allocate vector to estimate maximum number of buffer multifabs needed to
    // obtain the lab-frame snapshot.
    m_max_buffer_multifabs.resize(m_num_buffers);
    // allocate vector to count number of times the buffer multifab
    // has been flushed and refilled
    m_buffer_flush_counter.resize(m_num_buffers);
    // allocate vector of geometry objects corresponding to each snapshot
    m_geom_snapshot.resize( m_num_buffers );
    for (int i = 0; i < m_num_buffers; ++i) {
        m_geom_snapshot[i].resize(nmax_lev);
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
    amrex::ParmParse pp_diag_name(m_diag_name);

    m_file_prefix = "diags/" + m_diag_name;
    pp_diag_name.query("file_prefix", m_file_prefix);
    pp_diag_name.query("do_back_transformed_fields", m_do_back_transformed_fields);
    pp_diag_name.query("do_back_transformed_particles", m_do_back_transformed_particles);
    AMREX_ALWAYS_ASSERT(m_do_back_transformed_fields or m_do_back_transformed_particles);

    pp_diag_name.get("num_snapshots_lab", m_num_snapshots_lab);
    m_num_buffers = m_num_snapshots_lab;

    // Read either dz_snapshots_lab or dt_snapshots_lab
    bool snapshot_interval_is_specified = false;
    snapshot_interval_is_specified = queryWithParser(pp_diag_name, "dt_snapshots_lab", m_dt_snapshots_lab);
    if ( queryWithParser(pp_diag_name, "dz_snapshots_lab", m_dz_snapshots_lab) ) {
        m_dt_snapshots_lab = m_dz_snapshots_lab/PhysConst::c;
        snapshot_interval_is_specified = true;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(snapshot_interval_is_specified,
        "For back-transformed diagnostics, user should specify either dz_snapshots_lab or dt_snapshots_lab");
    // For BTD, we always need rho to perform Lorentz Transform of current-density
    if (WarpXUtilStr::is_in(m_cellcenter_varnames, "rho")) warpx.setplot_rho(true);

    if (pp_diag_name.query("buffer_size", m_buffer_size)) {
        if(m_max_box_size < m_buffer_size) m_max_box_size = m_buffer_size;
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
    amrex::IntVect hi(1);
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
    m_snapshot_box[i_buffer] = diag_box;
    // Define box array
    amrex::BoxArray diag_ba(diag_box);
    diag_ba.maxSize( warpx.maxGridSize( lev ) );
    // Update the physical co-ordinates m_lo and m_hi using the final index values
    // from the coarsenable, cell-centered BoxArray, ba.
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        diag_dom.setLo( idim, warpx.Geom(lev).ProbLo(idim) +
            diag_ba.getCellCenteredBox(0).smallEnd(idim) * warpx.Geom(lev).CellSize(idim));
        diag_dom.setHi( idim, warpx.Geom(lev).ProbLo(idim) +
            (diag_ba.getCellCenteredBox( diag_ba.size()-1 ).bigEnd(idim) + 1) * warpx.Geom(lev).CellSize(idim));
    }

    // Define buffer_domain in lab-frame for the i^th snapshot.
    // Replace z-dimension with lab-frame co-ordinates.
    amrex::Real zmin_buffer_lab = diag_dom.lo(m_moving_window_dir)
                                 / ( (1.0_rt + m_beta_boost) * m_gamma_boost);
    amrex::Real zmax_buffer_lab = diag_dom.hi(m_moving_window_dir)
                                 / ( (1.0_rt + m_beta_boost) * m_gamma_boost);


    m_snapshot_domain_lab[i_buffer] = diag_dom;
    m_snapshot_domain_lab[i_buffer].setLo(m_moving_window_dir,
                                  zmin_buffer_lab + warpx.moving_window_v * m_t_lab[i_buffer]);
    m_snapshot_domain_lab[i_buffer].setHi(m_moving_window_dir,
                                  zmax_buffer_lab + warpx.moving_window_v * m_t_lab[i_buffer]);

    // Initialize buffer counter and z-positions of the  i^th snapshot in
    // boosted-frame and lab-frame
    m_buffer_flush_counter[i_buffer] = 0;
    m_buffer_counter[i_buffer] = 0;
    m_current_z_lab[i_buffer] = 0._rt;
    m_current_z_boost[i_buffer] = 0._rt;
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
    m_snapshot_ncells_lab[i_buffer] = {Nx_lab, Ny_lab, Nz_lab};
#else
    m_snapshot_ncells_lab[i_buffer] = {Nx_lab, Nz_lab};
#endif
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
    m_cell_centered_data[lev] = std::make_unique<amrex::MultiFab>(ba, dmap,
                                     m_cellcenter_varnames.size(), ngrow);

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
                // Check if the zslice is in domain
                bool ZSliceInDomain = GetZSliceInDomainFlag (i_buffer, lev);
                // Initialize and define field buffer multifab if buffer is empty
                if (ZSliceInDomain) {
                    if ( buffer_empty(i_buffer) ) {
                        if ( m_buffer_flush_counter[i_buffer] == 0) {
                            // Compute the geometry, snapshot lab-domain extent
                            // and box-indices
                            DefineSnapshotGeometry(i_buffer, lev);
                        }
                        DefineFieldBufferMultiFab(i_buffer, lev);
                    }
                }
                m_all_field_functors[lev][i]->PrepareFunctorData (
                                             i_buffer, ZSliceInDomain,
                                             m_current_z_boost[i_buffer],
                                             m_buffer_box[i_buffer],
                                             k_index_zlab(i_buffer, lev), m_max_box_size );

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
                          ( m_current_z_lab[i_buffer]
                            - (prob_domain_zmin_lab + 0.5*dz_lab(warpx.getdt(lev), ref_ratio[m_moving_window_dir]) ) )
                          / dz_lab( warpx.getdt(lev), ref_ratio[m_moving_window_dir] )
                      ) );
    return k_lab;
}



void
BTDiagnostics::DefineFieldBufferMultiFab (const int i_buffer, const int lev)
{
    if ( m_do_back_transformed_fields ) {
        auto & warpx = WarpX::GetInstance();

        const int k_lab = k_index_zlab (i_buffer, lev);
        m_buffer_box[i_buffer].setSmall( m_moving_window_dir, k_lab - m_buffer_size+1);
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
        m_mf_output[i_buffer][lev].setVal(0.);

        amrex::IntVect ref_ratio = amrex::IntVect(1);
        if (lev > 0 ) ref_ratio = WarpX::RefRatio(lev-1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            amrex::Real cellsize;
            if (idim < AMREX_SPACEDIM-1) {
                cellsize = warpx.Geom(lev).CellSize(idim);
            } else {
                cellsize = dz_lab(warpx.getdt(lev), ref_ratio[m_moving_window_dir]);
            }
            amrex::Real buffer_lo = m_prob_domain_lab[i_buffer].lo(idim)
                                    + (buffer_ba.getCellCenteredBox(0).smallEnd(idim) ) * cellsize;
            amrex::Real buffer_hi = m_prob_domain_lab[i_buffer].lo(idim) +
                                          (buffer_ba.getCellCenteredBox( buffer_ba.size()-1 ).bigEnd(idim) + 1) * cellsize;
            m_buffer_domain_lab[i_buffer].setLo(idim, buffer_lo);
            m_buffer_domain_lab[i_buffer].setHi(idim, buffer_hi);
        }

        // Define the geometry object at level, lev, for the ith buffer.
        if (lev == 0) {
            // The extent of the physical domain covered by the ith buffer mf, m_mf_output
            // Default non-periodic geometry for diags
            amrex::Vector<int> BTdiag_periodicity(AMREX_SPACEDIM, 0);
            // Box covering the extent of the user-defined diag in the back-transformed frame
            amrex::Box domain = buffer_ba.minimalBox();
            // define the geometry object for the ith buffer using Physical co-oridnates
            // of m_buffer_domain_lab[i_buffer].
            m_geom_output[i_buffer][lev].define( domain, &m_buffer_domain_lab[i_buffer],
                                                 amrex::CoordSys::cartesian,
                                                 BTdiag_periodicity.data() );
        } else if (lev > 0 ) {
            // Refine the geometry object defined at the previous level, lev-1
            m_geom_output[i_buffer][lev] = amrex::refine( m_geom_output[i_buffer][lev-1],
                                                          warpx.RefRatio(lev-1) );
        }
    }
}


void
BTDiagnostics::DefineSnapshotGeometry (const int i_buffer, const int lev)
{
    if ( m_do_back_transformed_fields ) {
        auto & warpx = WarpX::GetInstance();
        const int k_lab = k_index_zlab (i_buffer, lev);
        // Box covering the extent of the user-defined diag in the back-transformed frame
        // for the ith snapshot
        // estimating the maximum number of buffer multifabs needed to obtain the
        // full lab-frame snapshot
        m_max_buffer_multifabs[i_buffer] = static_cast<int>( ceil (
            amrex::Real(m_snapshot_ncells_lab[i_buffer][m_moving_window_dir]) /
            amrex::Real(m_buffer_size) ) );
        // number of cells in z is modified since each buffer multifab always
        // contains a minimum m_buffer_size=256 cells
        int num_z_cells_in_snapshot = m_max_buffer_multifabs[i_buffer] * m_buffer_size;
        // Modify the domain indices according to the buffers that are flushed out
        m_snapshot_box[i_buffer].setSmall( m_moving_window_dir,
                                           k_lab - (num_z_cells_in_snapshot-1) );
        m_snapshot_box[i_buffer].setBig( m_moving_window_dir, k_lab);

        // Modifying the physical coordinates of the lab-frame snapshot to be
        // consistent with the above modified domain-indices in m_snapshot_box.
        amrex::IntVect ref_ratio = amrex::IntVect(1);
        amrex::Real new_lo = m_snapshot_domain_lab[i_buffer].hi(m_moving_window_dir) -
                             num_z_cells_in_snapshot *
                             dz_lab(warpx.getdt(lev), ref_ratio[m_moving_window_dir]);
        m_snapshot_domain_lab[i_buffer].setLo(m_moving_window_dir, new_lo);
        if (lev == 0) {
            // The extent of the physical domain covered by the ith snapshot
            // Default non-periodic geometry for diags
            amrex::Vector<int> BTdiag_periodicity(AMREX_SPACEDIM, 0);
            // define the geometry object for the ith snapshot using Physical co-oridnates
            // of m_snapshot_domain_lab[i_buffer], that corresponds to the full snapshot
            // in the back-transformed frame
            m_geom_snapshot[i_buffer][lev].define( m_snapshot_box[i_buffer],
                                                   &m_snapshot_domain_lab[i_buffer],
                                                   amrex::CoordSys::cartesian,
                                                   BTdiag_periodicity.data() );

        } else if (lev > 0) {
            // Refine the geometry object defined at the previous level, lev-1
            m_geom_snapshot[i_buffer][lev] = amrex::refine( m_geom_snapshot[i_buffer][lev-1],
                                                            warpx.RefRatio(lev-1) );
        }
    }
}

bool
BTDiagnostics::GetZSliceInDomainFlag (const int i_buffer, const int lev)
{
    auto & warpx = WarpX::GetInstance();
    const amrex::RealBox& boost_domain = warpx.Geom(lev).ProbDomain();

    amrex::Real buffer_zmin_lab = m_snapshot_domain_lab[i_buffer].lo( m_moving_window_dir );
    amrex::Real buffer_zmax_lab = m_snapshot_domain_lab[i_buffer].hi( m_moving_window_dir );
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
BTDiagnostics::Flush (int i_buffer)
{
    auto & warpx = WarpX::GetInstance();
    std::string file_name = m_file_prefix;
    if (m_format=="plotfile") {
        file_name = amrex::Concatenate(m_file_prefix,i_buffer,5);
        file_name = file_name+"/buffer";
    }
    bool isLastBTDFlush = ( ( m_max_buffer_multifabs[i_buffer]
                               - m_buffer_flush_counter[i_buffer]) == 1) ? true : false;
    bool const isBTD = true;
    double const labtime = m_t_lab[i_buffer];
    m_flush_format->WriteToFile(
        m_varnames, m_mf_output[i_buffer], m_geom_output[i_buffer], warpx.getistep(),
        labtime, m_output_species, nlev_output, file_name,
        m_plot_raw_fields, m_plot_raw_fields_guards, m_plot_raw_rho, m_plot_raw_F,
        isBTD, i_buffer, m_geom_snapshot[i_buffer][0], isLastBTDFlush);

    if (m_format == "plotfile") {
        MergeBuffersForPlotfile(i_buffer);
    }

    // Reset the buffer counter to zero after flushing out data stored in the buffer.
    ResetBufferCounter(i_buffer);
    IncrementBufferFlushCounter(i_buffer);
}

void BTDiagnostics::TMP_ClearSpeciesDataForBTD ()
{
    m_output_species.clear();
    m_output_species_names.clear();
}

void BTDiagnostics::MergeBuffersForPlotfile (int i_snapshot)
{
    auto & warpx = WarpX::GetInstance();
    const amrex::Vector<int> iteration = warpx.getistep();
    if (amrex::ParallelContext::IOProcessorSub()) {
        // Path to final snapshot plotfiles
        std::string snapshot_path = amrex::Concatenate(m_file_prefix,i_snapshot,5);
        // BTD plotfile have only one level, Level0.
        std::string snapshot_Level0_path = snapshot_path + "/Level_0";
        std::string snapshot_Header_filename = snapshot_path + "/Header";
        // Create directory only when the first buffer is flushed out.
        if (m_buffer_flush_counter[i_snapshot] == 0 ) {
            // Create Level_0 directory to store all Cell_D and Cell_H files
            if (!amrex::UtilCreateDirectory(snapshot_Level0_path, 0755) )
                amrex::CreateDirectoryFailed(snapshot_Level0_path);
        }

        // Path of the buffer recently flushed
        std::string BufferPath_prefix = snapshot_path + "/buffer";
        const std::string recent_Buffer_filepath = amrex::Concatenate(BufferPath_prefix,iteration[0]);
        // Header file of the recently flushed buffer
        std::string recent_Header_filename = recent_Buffer_filepath+"/Header";
        std::string recent_Buffer_Level0_path = recent_Buffer_filepath + "/Level_0";
        std::string recent_Buffer_FabHeaderFilename = recent_Buffer_Level0_path + "/Cell_H";
        // Read the header file to get the fab on disk string
        BTDMultiFabHeaderImpl Buffer_FabHeader(recent_Buffer_FabHeaderFilename);
        Buffer_FabHeader.ReadMultiFabHeader();
        if (Buffer_FabHeader.ba_size() > 1) amrex::Abort("BTD Buffer has more than one fabs.");
        // Every buffer that is flushed only has a single fab.
        std::string recent_Buffer_FabFilename = recent_Buffer_Level0_path + "/"
                                              + Buffer_FabHeader.FabName(0);
        // Existing snapshot Fab Header Filename
        std::string snapshot_FabHeaderFilename = snapshot_Level0_path + "/Cell_H";
        std::string snapshot_FabFilename = amrex::Concatenate(snapshot_Level0_path+"/Cell_D_",m_buffer_flush_counter[i_snapshot], 5);
        // Name of the newly appended fab in the snapshot
        std::string new_snapshotFabFilename = amrex::Concatenate("Cell_D_",m_buffer_flush_counter[i_snapshot],5);

        if ( m_buffer_flush_counter[i_snapshot] == 0) {
            std::rename(recent_Header_filename.c_str(), snapshot_Header_filename.c_str());
        Buffer_FabHeader.SetFabName(0, Buffer_FabHeader.fodPrefix(0),
                            new_snapshotFabFilename,
                    Buffer_FabHeader.FabHead(0));
        Buffer_FabHeader.WriteMultiFabHeader();
            std::rename(recent_Buffer_FabHeaderFilename.c_str(),
                        snapshot_FabHeaderFilename.c_str());
            std::rename(recent_Buffer_FabFilename.c_str(),
                        snapshot_FabFilename.c_str());
        } else {
            // Interleave Header file
            InterleaveBufferAndSnapshotHeader(recent_Header_filename,
                                              snapshot_Header_filename);
            InterleaveFabArrayHeader(recent_Buffer_FabHeaderFilename,
                                     snapshot_FabHeaderFilename,
                                     new_snapshotFabFilename);
            std::rename(recent_Buffer_FabFilename.c_str(),
                        snapshot_FabFilename.c_str());
        }
        // Destroying the recently flushed buffer directory since it is already merged.
        amrex::FileSystem::RemoveAll(recent_Buffer_filepath);

    } // ParallelContext if ends
    amrex::ParallelDescriptor::Barrier();
}

void
BTDiagnostics::InterleaveBufferAndSnapshotHeader ( std::string buffer_Header_path,
                                                   std::string snapshot_Header_path)
{
    BTDPlotfileHeaderImpl snapshot_HeaderImpl(snapshot_Header_path);
    snapshot_HeaderImpl.ReadHeaderData();

    BTDPlotfileHeaderImpl buffer_HeaderImpl(buffer_Header_path);
    buffer_HeaderImpl.ReadHeaderData();

    // Update timestamp of snapshot with recently flushed buffer
    snapshot_HeaderImpl.set_time( buffer_HeaderImpl.time() );
    snapshot_HeaderImpl.set_timestep( buffer_HeaderImpl.timestep() );

    amrex::Box snapshot_Box = snapshot_HeaderImpl.probDomain();
    amrex::Box buffer_Box = buffer_HeaderImpl.probDomain();
    amrex::IntVect box_lo(0);
    amrex::IntVect box_hi(1);
    // Update prob_lo with min of buffer and snapshot
    for (int idim = 0; idim < snapshot_HeaderImpl.spaceDim(); ++idim) {
        amrex::Real min_prob_lo = amrex::min(buffer_HeaderImpl.problo(idim),
                                             snapshot_HeaderImpl.problo(idim));
        amrex::Real max_prob_hi = amrex::max(buffer_HeaderImpl.probhi(idim),
                                             snapshot_HeaderImpl.probhi(idim));
        snapshot_HeaderImpl.set_problo(idim, min_prob_lo);
        snapshot_HeaderImpl.set_probhi(idim, max_prob_hi);
        // Update prob_hi with max of buffer and snapshot
        box_lo[idim] = amrex::min(buffer_Box.smallEnd(idim),
                                  snapshot_Box.smallEnd(idim));
        box_hi[idim] = amrex::max(buffer_Box.bigEnd(idim),
                                  snapshot_Box.bigEnd(idim));
    }
    amrex::Box domain_box(box_lo, box_hi);
    snapshot_HeaderImpl.set_probDomain(domain_box);

    // Increment numFabs
    snapshot_HeaderImpl.IncrementNumFabs();
    // The number of fabs in the recently written buffer is always 1.
    snapshot_HeaderImpl.AppendNewFabLo( buffer_HeaderImpl.FabLo(0));
    snapshot_HeaderImpl.AppendNewFabHi( buffer_HeaderImpl.FabHi(0));

    snapshot_HeaderImpl.WriteHeader();
}


void
BTDiagnostics::InterleaveFabArrayHeader(std::string Buffer_FabHeader_path,
                                        std::string snapshot_FabHeader_path,
                                        std::string newsnapshot_FabFilename)
{
    BTDMultiFabHeaderImpl snapshot_FabHeader(snapshot_FabHeader_path);
    snapshot_FabHeader.ReadMultiFabHeader();

    BTDMultiFabHeaderImpl Buffer_FabHeader(Buffer_FabHeader_path);
    Buffer_FabHeader.ReadMultiFabHeader();

    // Increment existing fabs in snapshot with the number of fabs in the buffer
    snapshot_FabHeader.IncreaseMultiFabSize( Buffer_FabHeader.ba_size() );
    snapshot_FabHeader.ResizeFabData();

    for (int ifab = 0; ifab < Buffer_FabHeader.ba_size(); ++ifab) {
        int new_ifab = snapshot_FabHeader.ba_size() - 1 + ifab;
        snapshot_FabHeader.SetBox(new_ifab, Buffer_FabHeader.ba_box(ifab) );
        // Set Name of the new fab using newsnapshot_FabFilename.
        snapshot_FabHeader.SetFabName(new_ifab, Buffer_FabHeader.fodPrefix(ifab),
                                                newsnapshot_FabFilename,
                                                Buffer_FabHeader.FabHead(ifab) );
        snapshot_FabHeader.SetMinVal(new_ifab, Buffer_FabHeader.minval(ifab));
        snapshot_FabHeader.SetMaxVal(new_ifab, Buffer_FabHeader.maxval(ifab));
    }

    snapshot_FabHeader.WriteMultiFabHeader();

}
