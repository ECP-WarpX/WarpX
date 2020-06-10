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
    /** Read relevant parameters for BTD */
    void ReadParameters ();
    /** Flush m_mf_output and particles to file. */
    void Flush () override {} ; 
    /** Flush raw data */
    void FlushRaw () override {};
    /** whether to flush at this time step
      * \param[in] step current time step
      * \param[in] force_flush if true, return true for any step
      *  The return bool is always set to true since data is back-transformed
      *  to lab-frame and added to the output buffer every timestep.
      *  The data is flushed when the buffer is full and/or when the simulation ends.
      */
    bool DoDump (int step, bool force_flush=false) override;
    /** Initialize functors that store pointers to the fields requested by the user.
      * \param[in] lev level on which the vector of unique_ptrs to field functors is initialized.
      */    
    void InitializeFieldFunctors (int lev){};
    /** Whether to compute back-tranformed values for field-data */ 
    bool m_do_back_transformed_fields = true;
    /** Whether to compute back-tranformed values for particle-data */ 
    bool m_do_back_transformed_particles = true;

    /** Map field names and component number in cell-centered_data */
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

    /** Map of field index in m_mf_output[i] that contained cell-centered BTD for
      * snapshots i. By default, all fields are dumped.
      */
    amrex::Gpu::ManagedDeviceVector<int> m_map_actual_fields_to_dump;
    
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
    /** Number of back-transformed snapshots in the lab-frame requested by the user */
    int m_num_snapshots_lab = std::numeric_limits<int>::lowest();
    /** Time interval in lab-frame between the back-transformed snapshots */
    amrex::Real m_dt_snapshots_lab = std::numeric_limits<amrex::Real>::lowest();
    /** Distance between the back-transformed snapshots in lab-frame
      * m_dz_snapshots_lab = m_dt_snapshots_lab/c
      */
    amrex::Real m_dz_snapshots_lab;
    
   
    int m_buffer_size = 256;
    int m_max_box_size = 256;

    /**  Time in lab-frame of the back-transformed snapshot */
    amrex::Vector<amrex::Real> m_t_lab;
    /**  Physical domain with lab-frame co-ordinates of the back-transformed snapshot */
    amrex::Vector<amrex::RealBox> m_domain_lab;
    /**  Number of cells in the lab-frame for all back-transformed snapshots */
    amrex::IntVect m_ncells_lab;
    /**  Box-dimension in boosted-frame index space for each back-transformed snapshot */
    amrex::Vector<amrex::Box> m_buffer_box;
    /**  Lab-frame z co-ordinate of the back-transformed snapshot at the current timestep */
    amrex::Vector<amrex::Real> m_current_z_lab;
    /**  Boosted-frame z co-ordinate corresponding to the back-transformed 
         snapshot at the current timestep */
    amrex::Vector<amrex::Real> m_current_z_boost;
    /**  Number of back-transformed z-slices in the output buffer multifab */
    amrex::Vector<int> m_buffer_counter;
   
    /** Customized meta-data for BTD */
    void writeMetaData ();
};

#endif // WARPX_BTDIAGNOSTICS_H_
