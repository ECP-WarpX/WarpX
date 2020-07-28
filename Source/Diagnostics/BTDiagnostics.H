#ifndef WARPX_BTDIAGNOSTICS_H_
#define WARPX_BTDIAGNOSTICS_H_

#include "Diagnostics.H"

class
BTDiagnostics final : public Diagnostics
{
public:

    BTDiagnostics (int i, std::string name);

private:
    /** Read relevant parameters for BTD */
    void ReadParameters ();
    /** Flush m_mf_output and particles to file. */
    // Currently this is an empty function. When OpenPMD/Plotfile is supported
    // this function can be moved to the base class Diagnostics.
    void Flush (int i_buffer) override {} ;
    /** whether to write output files at this time step
     *  The data is flushed when the buffer is full and/or
     *   when the simulation ends or when forced.
     * \param[in] step current time step
     * \param[in] force_flush if true, return true for any step
     *  The return bool is used to determine when to write buffer data out
     */
    bool DoDump (int step, int i_buffer, bool force_flush=false) override;
    bool DoComputeAndPack (int step, bool force_flush=false) override;
    /** Initialize Data required to compute back-transformed diagnostics */
    void DerivedInitData ();
    /** Initialize functors that store pointers to the fields requested by the user.
     * \param[in] lev level on which the vector of unique_ptrs to field functors is initialized.
     */
    void InitializeFieldFunctors (int lev);
    /** This function initialized and allocates particle buffers for all the snapshots.
     *  The particle containers required for this must be added to populate this function.
     */
    void InitializeParticleBuffer () override {};
    /** Compute the cell-centered data for all user-requested fields and store in
     *  the multi-level cell-centered multifab, m_mf_cc. This MultiFab extends
     *  over the entire domain and is coarsened using the user-defined crse_ratio
     */
    void PrepareFieldDataForOutput ();
    /** Initialize buffer domain, buffer box and lab-frame parameters such as
     *  m_t_lab, and z-positions for the i^th snapshot, i_buffer, and level, lev.
     * \param[in] i_buffer i^th snapshot or buffer
     * \param[in] lev level
     */
    void InitializeFieldBufferData ( int i_buffer , int lev) override;
    /** Whether to compute back-tranformed values for field-data */
    bool m_do_back_transformed_fields = true;
    /** Whether to compute back-tranformed values for particle-data */
    bool m_do_back_transformed_particles = true;

    // boost parameters
    amrex::Real m_gamma_boost;
    amrex::Real m_beta_boost;
    int m_moving_window_dir;

    /** Number of back-transformed snapshots in the lab-frame requested by the user */
    int m_num_snapshots_lab = std::numeric_limits<int>::lowest();
    /** Time interval in lab-frame between the back-transformed snapshots */
    amrex::Real m_dt_snapshots_lab = std::numeric_limits<amrex::Real>::lowest();
    /** Distance between the back-transformed snapshots in lab-frame
      * m_dz_snapshots_lab = m_dt_snapshots_lab/c
      */
    amrex::Real m_dz_snapshots_lab;

    /** Number of z-slices in each buffer of the snapshot */
    int m_buffer_size = 256;
    /** max grid size used to generate BoxArray to define output MultiFabs */
    int m_max_box_size = 256;

    /**  Time in lab-frame of the back-transformed snapshot */
    amrex::Vector<amrex::Real> m_t_lab;
    /**  Physical domain with lab-frame co-ordinates of the back-transformed snapshot */
    amrex::Vector<amrex::RealBox> m_buffer_domain_lab;
    /**  Number of cells in the lab-frame for all back-transformed snapshots */
    amrex::Vector<amrex::IntVect> m_buffer_ncells_lab;
    /**  Box-dimension in boosted-frame index space for each back-transformed snapshot */
    amrex::Vector<amrex::Box> m_buffer_box;
    /**  Lab-frame z co-ordinate of the back-transformed snapshot at the current timestep */
    amrex::Vector<amrex::Real> m_current_z_lab;
    /**  Boosted-frame z co-ordinate corresponding to the back-transformed
         snapshot at the current timestep */
    amrex::Vector<amrex::Real> m_current_z_boost;
    /**  Number of back-transformed z-slices in the output buffer multifab */
    amrex::Vector<int> m_buffer_counter;
    /** Multi-level cell-centered multifab with all field-data components
        which will be sliced and transformed. The same multifab will be used
        for all snapshots.
      */
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_cell_centered_data;
    std::vector<std::string> m_file_name;
    /** Vector of pointers to compute cell-centered data, per level, per component
     *  using the coarsening-ratio provided by the user.
     */
    amrex::Vector< amrex::Vector <std::unique_ptr<ComputeDiagFunctor const> > > m_cell_center_functors;

    /** Customized function to write back-transformed data */
    void TMP_writeMetaData ();
    /** Customized functions to write back-transformed data */
    void TMP_createLabFrameDirectories (int i_buffer, int lev);
    /** Define the cell-centered multi-component MultiFab at level, lev.
     *  \param[in] lev level at which the cell-centered MultiFab is defined.
     */
    void DefineCellCenteredMultiFab(int lev);
    /** Returns z-position in the boosted-frame
      * \param[in] lab-frame time of the snapshot
      * \param[in] boosted-frame time at level, lev
      */
    amrex::Real UpdateCurrentZBoostCoordinate(amrex::Real t_lab, amrex::Real t_boost)
    {
        amrex::Real current_z_boost = (t_lab / m_gamma_boost - t_boost) * PhysConst::c / m_beta_boost;
        return current_z_boost;
    }
    /** Returns z-position in the lab-frame
      * \param[in] lab-frame time of the snapshot
      * \param[in] boosted-frame time at level, lev
      */
    amrex::Real UpdateCurrentZLabCoordinate(amrex::Real t_lab, amrex::Real t_boost)
    {
        amrex::Real current_z_lab = (t_lab - t_boost / m_gamma_boost ) * PhysConst::c / m_beta_boost;
        return current_z_lab;
    }
    /** Returns cell-size of the grid in z-direction in the lab-frame at level, lev
     *  \param[in] dt in boosted-frame at level, lev
     *  \param[in] ref_ratio, refinement ratio at level, lev
     */
    amrex::Real dz_lab (amrex::Real dt, amrex::Real ref_ratio);
    /** whether field buffer is full
     * \param[in] i_buffer buffer id for which the buffer size is checked.
     * returns bool = true is buffer is full
     */
    bool buffer_full (int i_buffer) {
        return ( m_buffer_counter[i_buffer] == m_buffer_size );
    }
};

#endif // WARPX_BTDIAGNOSTICS_H_
