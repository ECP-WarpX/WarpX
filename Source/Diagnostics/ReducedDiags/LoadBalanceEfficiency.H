/* Copyright 2020-2021 Michael Rowan, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_LOADBALANCEEFFICIENCY_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_LOADBALANCEEFFICIENCY_H_

#include "WarpX.H"
#include "ReducedDiags.H"


/**
 *  This class mainly contains a function that gets the
 *  current load balance efficiency for writing to output.
 */
class LoadBalanceEfficiency : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    LoadBalanceEfficiency(std::string rd_name);

    /** This function gets the current load balance efficiency  */
    virtual void ComputeDiags(int step) override final;
};

#endif
