#ifndef WARPX_BACKTRANSFORMFUNCTOR_H_
#define WARPX_BACKTRANSFORMFUNCTOR_H_

#include "ComputeDiagFunctor.H"

/**
 * \brief Functor to back-transform cell-centered data and store result in mf_out
 *
 * The cell-centered data is a ten-component multifab with field-data averaged-down
 * from the finest to coarsest level, and stored as single-level data.
 * For every i^th buffer, a z-slice corresponding to the z-boost location of the
 * slice at the current timestep is extracted. This slice containing field-data
 * in the boosted-frame is Lorentz-transformed to the lab-frame. The user-requested
 * lab-frame field data is then stored in mf_dst.
 */

class
BackTransformFunctor final : public ComputeDiagFunctor
{
public:
    /** Constructor description
     *
     * \param[in] mf_src cell-centered multifab containing all user-requested
                         fields in boosted-frame
     * \param[in] lev level of multifab.
     * \param[in] ncom number of components of mf_src to Lorentz-Transform
                       and store in destination multifab.
     */
    BackTransformFunctor ( const amrex::MultiFab * const mf_src, const int lev,
                           const int ncomp, const int num_buffers,
                           amrex::Vector< std::string > varnames,
                           const amrex::IntVect crse_ratio= amrex::IntVect(1));

    /** \brief Lorentz-transform mf_src for the ith buffer and write the result in mf_dst.
     *
     * The source multifab, is a ten-component cell-centered multifab storing
     * field-data in the boosted-frame. An z-slice is generated
     * at the z-boost location for the ith buffer, stored in m_current_z_boost[i_buffer].
     * The data is then lorentz-transformed in-place using LorenzTransformZ ().
     * The user-requested fields are then copied to mf_dst.
     *
     * \param[out] mf_dst output MuliFab where the back-transformed data is written
     * \param[in] dcomp first component of mf_dst in which the back-transformed
     *            lab-frame data for the user-request fields is written.
     * \param[in] i_buffer buffer index for which the data is transformed.
     */
    void operator ()(amrex::MultiFab& mf_dst, int dcomp, const int i_buffer) const override;

    /** \brief Prepare data required to back-transform fields for lab-frame snapshot, i_buffer
     *
     * \param[in] i_buffer, index of the snapshot
     * \param[in] ZSliceInDomain if the z-slice at current_z_boost is within
     *            the boosted-frame and lab-frame domain.
     *            The fields are sliced and back-transformed only if this value is true.
     * \param[in] buffer_box, Box with index-space in lab-frame for the ith buffer
     * \param[in] k_index_zlab, k-index in the lab-frame corresponding to the
     *                          current z co-ordinate in the lab-frame for the ith buffer.
     */
    void PrepareFunctorData ( int i_buffer, bool ZSliceInDomain,
                              amrex::Real current_z_boost,
                              amrex::Box buffer_box, const int k_index_zlab ) override;
    /** Allocate and initialize member variables and arrays required to back-transform
     *  field-data from boosted-frame to lab-frame.
     */
    void InitData () override;
    /** \brief In-place Lorentz-transform of MultiFab, data, from boosted-frame to the
     *  lab-frame for all fields, Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, and rho.
     *
     *  \param[in] data z-slice field-data MultiFab to be back-transformed in-place
     *                 from boosted-frame to lab-frame
     *  \param[in] gamma_boost The Lorentz factor of the boosted-frame in which
                       the simulation is run
     *  \param[in] beta_boost The ratio of boost velocity to the speed of light
     */
    void LorentzTransformZ (amrex::MultiFab& data, amrex::Real gamma_boost,
                           amrex::Real beta_boost) const;
private:
    /** pointer to source multifab (cell-centered multi-component multifab) */
    amrex::MultiFab const * const m_mf_src = nullptr;
    /** level at which m_mf_src is defined */
    int const m_lev;
    /** Number of buffers or snapshots */
    int const m_num_buffers;
    /** Vector of amrex::Box with index-space in the lab-frame */
    amrex::Vector<amrex::Box> m_buffer_box;
    /** Vector of current z co-ordinate in the boosted-frame for each buffer */
    amrex::Vector<amrex::Real> m_current_z_boost;
    /** Vector of 0s and 1s stored to check if back-transformation is to be performed
     *  for the ith buffer. The value is set 0 (false) or 1 (true) using the boolean
     *  ZSliceInDomain in PrepareFunctorData().
     */
    amrex::Vector<int> m_perform_backtransform;
    /** Vector of k-index correspoding to the current lab-frame z co-ordinate for each buffer */
    amrex::Vector<int> m_k_index_zlab;
    /** Vector of user-defined field names to be stored in the output multifab */
    amrex::Vector< std::string > m_varnames;

    /** max grid size used to generate BoxArray to define output MultiFabs */
    int m_max_box_size = 256;
    /** Indices that map user-defined fields to plot to the fields
     *  stored in the cell-centered MultiFab, m_mf_src.
     *  The cell-centered MultiFab stores Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, and rho.
     */
    amrex::Vector<int> m_map_varnames;
};

#endif
