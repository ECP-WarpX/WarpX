#ifndef WARPX_BTDIAGNOSTICS_H_
#define WARPX_BTDIAGNOSTICS_H_

#include "Diagnostics.H"

class
BTDiagnostics final : public Diagnostics
{
public:
    BTDiagnostics (int i, std::string name);
    void InitData () override;
    

private:
    void ReadParameters ();
    /** Flush m_mf_output and particles to file. */
    void Flush () override {} ; 
    /** Flush raw data */
    void FlushRaw () override {} ;
    bool DoDump (int step, bool force_flush=false) override;
    void InitializeFieldFunctors (int lev){};
    
    bool m_do_back_transformed_fields = true;
    bool m_do_back_transformed_particles = true;
    // currently not used yet //
    std::string m_lab_data_directory = "lab_frame_data";


    // Map field names and component number in cell_centered_data
    std::map<std::string, int> m_possible_fields_to_dump = {
        {"Ex" , 0},
        {"Ey" , 1},
        {"Ez" , 2},
        {"Bx" , 3},
        {"By" , 4},
        {"Bz" , 5},
        {"jx" , 6},
        {"jy" , 7},
        {"jz" , 8},
        {"rho", 9} };

    // maps field index in data_buffer_[i] -> cell_centered_data for
    // snapshots i. By default, all fields in cell_centered_data are dumped.
    // Needs to be amrex::Vector because used in a ParallelFor kernel.
    amrex::Gpu::ManagedDeviceVector<int> m_map_actual_fields_to_dump;
    // Name of fields to dump. By default, all fields in cell_centered_data.
    // Needed for file headers only.
    std::vector<std::string> m_mesh_field_names = {"Ex", "Ey", "Ez",
                                                 "Bx", "By", "Bz",
                                                 "jx", "jy", "jz", "rho"};
    int m_ncomp_to_dump = 10;


    // BTD related variables
    amrex::Real m_gamma_boost;
    amrex::Real m_inv_gamma_boost;
    amrex::Real m_beta_boost;
    amrex::Real m_inv_beta_boost;
    amrex::Real m_dz_lab;
    amrex::Real m_inv_dz_lab;
    amrex::Real m_dt_boost;
    int m_boost_direction;    
    int m_moving_window_dir;
    
    int m_num_snapshots_lab = std::numeric_limits<int>::lowest();
    amrex::Real m_dt_snapshots_lab = std::numeric_limits<amrex::Real>::lowest();
    amrex::Real m_dz_snapshots_lab;

    int m_num_buffer = 256;
    int m_max_box_size = 256;
   
    void writeMetaData ();
//    // common with full diags
//    std::string m_diag_name;
//    std::string m_file_prefix;
//    amrex::Vector< std::string > m_varnames;
};

#endif // WARPX_BTDIAGNOSTICS_H_
