#include "BTDiagnostics.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
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

void
BTDiagnostics::InitData ()
{
    InitBaseData();
    amrex::Print() << " init data for BTD \n";
    // BTD related variables
    auto & warpx = WarpX::GetInstance();
    m_gamma_boost = warpx.gamma_boost;
    m_moving_window_dir = warpx.moving_window_dir;
    m_dt_boost = warpx.getdt(0);

    m_inv_gamma_boost = 1._rt / m_gamma_boost;
    m_beta_boost = std::sqrt( 1._rt - m_inv_gamma_boost * m_inv_gamma_boost);
    m_inv_beta_boost = 1._rt / m_beta_boost;
    m_dz_lab = PhysConst::c * m_dt_boost * m_inv_beta_boost * m_inv_gamma_boost;
    m_inv_dz_lab = 1._rt / m_dz_lab;    


    writeMetaData();


    // The ncells are required for writing the Header file //
    // Note : the cell size must be different for each level //
    amrex::Real zmin_lab = m_lo[m_moving_window_dir]
                           / ( (1.0_rt + m_beta_boost) * m_gamma_boost);
    amrex::Real zmax_lab = m_hi[m_moving_window_dir]
                           / ( (1.0_rt + m_beta_boost) * m_gamma_boost);

    int Nz_lab = static_cast<int>( (zmax_lab - zmin_lab) * m_inv_dz_lab );
    int Nx_lab = static_cast<int>(m_hi[0] - m_lo[0])/warpx.Geom(0).CellSize(0);
#if (AMREX_SPACEDIM == 3)
    // compute ny_lab m_hi - m_lo 
    int Ny_lab = static_cast<int>(m_hi[1] - m_lo[1])/warpx.Geom(0).CellSize(1);
    amrex::IntVect prob_ncells_lab = {Nx_lab, Ny_lab, Nz_lab};
#else
    amrex::IntVect prob_ncells_lab = {Nx_lab, Nz_lab};
#endif    


   


    // allocate vector of m_t_lab (m_num_snapshots_lab) ;
    m_t_lab.resize(m_num_snapshots_lab);
    // allocate vector of RealBox of the diag domain
    m_domain_lab.resize(m_num_snapshots_lab);
    // define box correctly (one for all snapshots)
    m_buffer_box.resize(m_num_snapshots_lab);
    // allocate vector of m_current_z_lab, m_current_z_boost
    m_current_z_lab.resize(m_num_snapshots_lab);
    m_current_z_boost.resize(m_num_snapshots_lab);
    // allocate vector of m_buff_counter
    m_buffer_counter.resize(m_num_snapshots_lab);

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

    m_map_actual_fields_to_dump.resize(m_varnames.size());
    for (int i=0; i < m_varnames.size(); ++i) {
        m_map_actual_fields_to_dump[i] = m_possible_fields_to_dump[ m_varnames[i] ];
    }

}

void
BTDiagnostics::writeMetaData ()
{
    WARPX_PROFILE("BTDiagnostic::WriteMetaData");
    
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
BTDiagnostics::DoDump (int step, bool force_flush)
{
    return true;
}
