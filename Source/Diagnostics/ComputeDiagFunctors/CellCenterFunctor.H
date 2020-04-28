#ifndef WARPX_CELLCENTERFUNCTOR_H_
#define WARPX_CELLCENTERFUNCTOR_H_

#include "ComputeDiagFunctor.H"

struct FieldTypes { enum TypeEnum { E=0, B, J, rho, F, divE, divB }; };

/**
 * \brief Functor to cell-center MF and store result in mf_out.
 */
template <FieldTypes::TypeEnum FIELDTYPE>
class
CellCenterFunctor : public ComputeDiagFunctor
{
public:

    /** Constructor.
     * \param[in] lev level of multifab. Used for averaging in rz.
     * \param[in] dir direction of field to use
     * \param[in] convert_RZmodes2cartesian (in cylindrical) whether to
     *            sum all modes in mf_src before cell-centering into dst multifab.
     * \param[in] ncomp Number of component of mf_src to cell-center in dst multifab.
     */
    CellCenterFunctor(const int lev, const int dir=0,
                      const bool convertRZmodes2cartesian=true, const int ncomp=1);

    /** \brief Fetch the field, and cell-center it, and write the result in mf_dst.
     *
     * In cylindrical geometry, by default this functor average all components
     * of a MultiFab and writes into one single component.
     *
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which cell-centered
     *            data is stored
     */
     void operator()(amrex::MultiFab& mf_dst, int dcomp, const amrex::IntVect crse_ratio) const;


protected:

    int m_lev; /**< level on which mf_src is defined (used in cylindrical) */
    /**< (for cylindrical) whether to average all modes into 1 comp */
    int m_dir;
    bool m_convertRZmodes2cartesian;

};

#endif // WARPX_CELLCENTERFUNCTOR_H_
