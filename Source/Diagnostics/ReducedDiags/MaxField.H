/* Copyright 2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_MAXFIELD_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_MAXFIELD_H_

#include "ReducedDiags.H"

/**
 *  This class mainly contains a function that computes the maximum value of each component
 * of the EM field and of the vector norm of the E and B fields.
 */
class MaxField : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    MaxField(std::string rd_name);

    /** This function computes the maximum value of Ex, Ey, Ez, |E|, Bx, By, Bz and |B| */
    virtual void ComputeDiags(int step) override final;

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_MAXFIELD_H_
