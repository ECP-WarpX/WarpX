#ifndef WARPX_COMPUTEDIAGFUNCTOR_H_
#define WARPX_COMPUTEDIAGFUNCTOR_H_

#include "ComputeDiagFunctor_fwd.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>

/**
 * \brief Functor to compute a diagnostic and store the result in existing
 * MultiFab
 */
class
ComputeDiagFunctor
{
public:
    ComputeDiagFunctor( int ncomp, amrex::IntVect crse_ratio) :
                        m_ncomp(ncomp), m_crse_ratio(crse_ratio) {}
    //** Virtual Destructor to handle clean destruction of derived classes */
    virtual ~ComputeDiagFunctor() = default;
    /** Compute a field and store the result in mf_dst
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which the result is written
     * \param[in] ncomp number of components of mf_dst to update
     * \param[in] crse_ratio for interpolating field values from simulation MultiFab, mf_src,
                  to the output diagnostic MultiFab, mf_dst.
     */
    virtual void operator() (amrex::MultiFab& mf_dst, int dcomp, const int i_buffer = 0) const = 0;
    /** Number of component from the input multifab to write to the output
     * multifab */
    int nComp () const { return m_ncomp; }

    /** \brief Prepare data required to process fields in the operator()
     *         Note that this function has parameters that are specific to
     *         back-transformed diagnostics, that are unused for regular diagnostics.
     *
     * \param[in] i_buffer       index of the back-transform snapshot
     * \param[in] z_slice_in_domain if the z-slice at current_z_boost is within
     *            the boosted-frame and lab-frame domain.
     *            The fields are sliced and back-transformed only if this value is true.
     * \param[in] buffer_box     Box with index-space in lab-frame for the ith buffer
     * \param[in] k_index_zlab   k-index in the lab-frame corresponding to the
     *                           current z co-ordinate in the lab-frame for the ith buffer.
     * \param[in] max_box_size   maximum box size for the multifab to generate box arrays
     * \param[in] snapshot_full   if the current snapshot, with index, i_buffer, is full (1)
                                 or not (0). If it is full, then Lorentz-transform
                                 is not performed by setting m_perform_backtransform to 0;
     */
    virtual void PrepareFunctorData ( int i_buffer, bool z_slice_in_domain,
                                      amrex::Real current_z_boost,
                                      amrex::Box buffer_box, const int k_index_zlab,
                                      const int max_box_size, const int snapshot_full) {
                                          amrex::ignore_unused(i_buffer,
                                          z_slice_in_domain,
                                          current_z_boost, buffer_box,
                                          k_index_zlab, max_box_size, snapshot_full);
                                      }
    virtual void InitData() {}
private:
    /** Number of components of mf_dst that this functor updates. */
    int m_ncomp;
protected:
    /** Coarsening ratio used to interpolate fields from simulation MultiFabs to output MultiFab. */
    amrex::IntVect m_crse_ratio;
};

#endif // WARPX_COMPUTEDIAGFUNCTOR_H_
