/* Copyright 2020 Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_RHOMAXIMUM_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_RHOMAXIMUM_H_

#include "ReducedDiags.H"
#include "Diagnostics/ComputeDiagFunctors/ComputeDiagFunctor.H"

/**
 *  This class mainly contains a function that computes the extrema of the total charge density
 *  and of the charge density of each charged species.
 */
class RhoMaximum : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names
     */
    RhoMaximum(std::string rd_name);

    /** This function computes the maximum and minimum values of rho (summed over all species) and
     * the maximum absolute value of rho for each species.
     */
    virtual void ComputeDiags(int step) override final;

private:
    /** Vector of (pointers to) functors to compute rho, per level, per species. We reuse here the
     * same functors as those used for regular diagnostics.
     *  */
    amrex::Vector< amrex::Vector <std::unique_ptr<ComputeDiagFunctor > > > m_rho_functors;

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_RHOMAXIMUM_H_
