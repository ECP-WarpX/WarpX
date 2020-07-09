#ifndef WARPX_BACKTRANSFORMFUNCTOR_H_
#define WARPX_BACKTRANSFORMFUNCTOR_H_

#include "ComputeDiagFunctor.H"

/**
 * \brief Functor to back-transform cell-centered data and store result in mf_out
 */

class
BackTransformFunctor final : public ComputeDiagFunctor
{
public:
    // This class will handle the slicing of the cell-centered field-data (stored as mf_src)
    // It will Lorentz-Transform the sliced data from boosted-frame at z=zboost
    // to lab-frame at z=zlab and add the lab-frame slice to the destination
    // multifab that will be provided in the operator () function, for all ncomponents
    // starting from icomp_dst.
    // The operator() function will also need a few boosted-frame parameters, such as,
    // current_z_boost of the corresponding snapshot, to generate appropriate slices
    // in boosted-frame and lab-frame.i
    // Unlike cell-center functor, coarsening ratio is not provided to this functor,
    // since the data is coarsened when preparing fields in m_mf_cc (mf_src).

    /** Constructor description
     * \param[in] mf_src cell-centered multifab containing all user-requested
                         fields in boosted-frame
     * \param[in] lev level of multifab.
     * \param[in] ncom number of components of mf_src to Lorentz-Transform
                       and store in destination multifab.
     */

    // initializing crse_ratio with 1, to accurately initialize ComputeDiagFunctor

    BackTransformFunctor ( const amrex::MultiFab * const mf_src, const int lev,
                           const int ncomp, const int num_buffers,
                           amrex::Vector< std::string > varnames,
                           const amrex::IntVect crse_ratio= amrex::IntVect(1));

    void operator ()(amrex::MultiFab& mf_dst, int dcomp, const int i_buffer) const override;
 
    /** Prepare data for back-transform functor for each lab-frame snapshot, i_buffer
     */
    void PrepareFunctorData ( int i_buffer, bool ZSliceInDomain,
                              amrex::Real current_z_boost,
                              amrex::Box buffer_box, const int k_index_zlab ) override;
    void InitData () override;
    /** In-place conversion of MultiFab data from boosted-frame to the 
     *  back-transformed lab-frame/
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
    int const m_num_buffers;
    amrex::Vector<amrex::Box> m_buffer_box;
    amrex::Vector<amrex::Real> m_current_z_boost;
    amrex::Vector<int> m_perform_backtransform;
    amrex::Vector<int> m_k_index_zlab;
    amrex::Vector< std::string > m_varnames;

    /** max grid size used to generate BoxArray to define output MultiFabs */
    int m_max_box_size = 256;
    /** Indices that map user-defined fields to plot to the fields
     *  stored in the cell-centered MultiFab, m_mf_src.
     */
    amrex::Gpu::ManagedDeviceVector<int> m_map_varnames;
};

#endif
