/* Copyright 2023 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_CHARGEINSIDEBOUNDARY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_CHARGEINSIDEBOUNDARY_H_

#include "ReducedDiags.H"

#include <string>

/**
 *  This class mainly contains a function that
 *  computes the charge inside a given boundary, by
 *  using Gauss's theorem: $Q_{int} = \epsilon_0 \iint dS \cdot E$
 *  where the integral is performed over that boundary
 */
class ChargeInsideBoundary : public ReducedDiags
{
public:

    /**
     * constructor
     * @param[in] rd_name reduced diags names
     */
    ChargeInsideBoundary(std::string rd_name);

    /**
     * This function computes the charge inside a given boundary:
     * $Q_{int} = \epsilon_0 \iint dS \cdot E$
     *  where the integral is performed over that boundary
     *
     * @param[in] step current time step
     */
    virtual void ComputeDiags(int step) override final;

};

#endif
